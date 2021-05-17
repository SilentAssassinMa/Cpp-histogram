// This project is CppHistogram, and licenced under the terms of the
// DO WHAT THE FUCK YOU WANT TO PUBLIC LICENCE, version 3.1,
// as published by dtf on July 2019. See the LICENCSE.txt file or
// https://ph.dtf.wtf/w/wtfpl/#version-3-1 for more details.

// Author : SilentAssassinMa
// Contact : silentassassin@mail.ustc.edu.cn
// May be modified later
// Comments on more features are welcome
//
// Basic Usage :
// 1D histogram with fixed bin width : CppHistogram myhist(Axis::Regular(100, 0, 100));
// 2D histogram with varied bin width : CppHistogram myhist(Axis::Variable({0, 2, 5, 6}), Axis::Variable({0, 2, 5, 9}));
// in principle you can define histogram with any dimension.
// Now Fill(), GetBinContent(), GetBinError(), GetBin(), GetAxis<>() are supported, similar usage as THXD in ROOT.
// Also a stats counter of type CppCounter<dim> is supported, one can get the statistics of the filled histogram.
// CppCounter provides a series of functions like GetMean(), GetStdVar() ..., see CppCounter.h for more details.
#ifndef CPPHISTOGRAM_HEADER
#define CPPHISTOGRAM_HEADER
#include <iostream>
#include <tuple>
#include <vector>
#include <algorithm>
#include <functional>
#include "CppCounter.h"

template <
    size_t Index = 0,                                                      // start iteration at 0 index
    typename TTuple,                                                       // the tuple type
    size_t Size = std::tuple_size<std::remove_reference_t<TTuple>>::value, // tuple size
    typename TCallable,                                                    // the callable to bo invoked for each tuple item
    typename... TArgs                                                      // other arguments to be passed to the callable
    >
inline void tuple_for_each(TTuple &&tuple, TCallable &&callable, TArgs &&...args)
{
    if constexpr (Index < Size)
    {
        std::invoke(callable, args..., std::get<Index>(tuple));

        if constexpr (Index + 1 < Size)
            tuple_for_each<Index + 1>(
                std::forward<TTuple>(tuple),
                std::forward<TCallable>(callable),
                std::forward<TArgs>(args)...);
    }
}

enum class OverflowOption
{
    WithOverflows,
    WithoutOverflows
};

template <typename Axis, typename... MoreAxis>
class CppHistogram
{
    using axis_type = std::tuple<std::decay_t<Axis>, std::decay_t<MoreAxis>...>;
    using counter_type = CppCounter<sizeof...(MoreAxis) + 1>;

public:
    CppHistogram(Axis &&axis, MoreAxis &&...moreaxis)
        : myaxis(std::forward<Axis>(axis), std::forward<MoreAxis>(moreaxis)...),
          dimension(sizeof...(MoreAxis) + 1)
    {
        int contentsize = 1;
        tuple_for_each(myaxis, [&](auto &item) { contentsize *= item.Nbins(OverflowOption::WithOverflows); });
        content = std::vector<double>(contentsize);
        error = std::vector<double>(contentsize);
    }

    template <typename... FillPar>
    int Fill(FillPar &&...par)
    {
        constexpr std::size_t naxis = sizeof...(MoreAxis) + 1;
        constexpr std::size_t npar = sizeof...(FillPar);
        static_assert(npar == naxis || npar == naxis + 1, "must be match");
        int binnumber = 0;
        bool outrange = false;
        double weight = 1;
        if constexpr (naxis == npar)
            FindBin(std::forward_as_tuple(par...), binnumber, outrange);
        if constexpr (naxis + 1 == npar)
            FindBinWithWeight(std::forward_as_tuple(par...), binnumber, weight, outrange);
        content[binnumber] += weight;
        error[binnumber] += weight * weight;
        if (!outrange || StatsOption == OverflowOption::WithOverflows)
            StatsCounter.Add(par...);
        return binnumber;
    }

    template <typename... BinIndex>
    double GetBinContent(BinIndex &&...indices)
    {
        static_assert(sizeof...(BinIndex) == std::tuple_size<std::remove_reference_t<axis_type>>::value, "must be match");
        return content[GetBin_impl(std::forward_as_tuple(indices...))];
    }

    template <typename... BinIndex>
    double GetBinError(BinIndex &&...indices)
    {
        static_assert(sizeof...(BinIndex) == std::tuple_size<std::remove_reference_t<axis_type>>::value, "must be match");
        return sqrt(fabs(error[GetBin_impl(std::forward_as_tuple(indices...))]));
    }

    template <typename... BinIndex>
    int GetBin(BinIndex &&...indices)
    {
        static_assert(sizeof...(BinIndex) == std::tuple_size<std::remove_reference_t<axis_type>>::value, "must be match");
        return GetBin_impl(std::forward_as_tuple(indices...));
    }

    template <size_t Index>
    // decltype(auto) GetAxis()
    auto GetAxis() -> const typename std::tuple_element<Index, axis_type>::type &
    {
        static_assert(Index < std::tuple_size<std::remove_reference_t<axis_type>>::value, "index must be lower than dimension");
        return std::cref(std::get<Index>(myaxis));
    }

    void SetStatOverflows(OverflowOption option = OverflowOption::WithOverflows)
    {
        this->StatsOption = option;
    }

    const counter_type &GetStats() { return StatsCounter; }

private:
    template <
        typename TTuple,                                                            // the tuple type
        size_t Index = std::tuple_size<std::remove_reference_t<TTuple>>::value - 1, // start iteration at Size - 1 index
        size_t Size = std::tuple_size<std::remove_reference_t<TTuple>>::value       // tuple size
        >
    inline void FindBin(TTuple &&value, int &binnumber, bool &outrange)
    {
        if constexpr (Index >= 0)
        {
            if constexpr (Index == Size - 1)
            {
                int temp = std::get<Index>(myaxis).FindBin(std::get<Index>(value));
                binnumber += temp;
                if (temp == 0 || temp > std::get<Index>(myaxis).Nbins())
                    outrange = true;
            }

            if constexpr (Index < Size - 1)
            {
                binnumber *= std::get<Index>(myaxis).Nbins(OverflowOption::WithOverflows);
                int temp = std::get<Index>(myaxis).FindBin(std::get<Index>(value));
                binnumber += temp;
                if (temp == 0 || temp > std::get<Index>(myaxis).Nbins())
                    outrange = true;
            }

            if constexpr (Index >= 1)
                FindBin<TTuple, Index - 1>(std::forward<TTuple>(value), std::ref(binnumber), outrange);
        }
    }

    template <
        typename TTuple,                                                            // the tuple type
        size_t Index = std::tuple_size<std::remove_reference_t<TTuple>>::value - 2, // start iteration at Size - 2 index
        size_t Size = std::tuple_size<std::remove_reference_t<TTuple>>::value       // tuple size
        >
    inline void FindBinWithWeight(TTuple &&value, int &binnumber, double &weight, bool &outrange)
    {
        // the last of the value is the weight
        weight = std::get<Size - 1>(value);
        if constexpr (Index >= 0)
            FindBin<TTuple, Index, Size - 1>(std::forward<TTuple>(value), std::ref(binnumber), outrange);
    }

    template <
        typename TTuple,                                                      // the tuple type
        size_t Index = 0,                                                     // start iteration at Size - 1 index
        size_t Size = std::tuple_size<std::remove_reference_t<TTuple>>::value // tuple size
        >
    inline int GetBin_impl(TTuple &&value)
    {
        if constexpr (Index == Size - 1)
            return std::get<Index>(value);
        else
            return std::get<Index>(myaxis).Nbins(OverflowOption::WithOverflows) * GetBin_impl<TTuple, Index + 1>(std::forward<TTuple>(value)) + std::get<Index>(value);
    }

    axis_type myaxis;
    std::size_t dimension;
    std::vector<double> content;
    std::vector<double> error;
    OverflowOption StatsOption = OverflowOption::WithoutOverflows;

    counter_type StatsCounter;
};

namespace Axis
{
    class Regular
    {
    public:
        Regular(unsigned int _nbins_, double _min_, double _max_)
            : nbins(_nbins_ + 2), min(_min_), max(_max_)
        {
            width = (max - min) / _nbins_;
        }

        int FindBin(double val) const
        {
            int __binnum = val < min ? 0 : (val - min) / width + 1;
            return __binnum > nbins - 2 ? nbins - 1 : __binnum;
        }

        int Nbins(OverflowOption option = OverflowOption::WithoutOverflows) const
        {
            if (option == OverflowOption::WithoutOverflows)
                return nbins - 2;
            else
                return nbins;
        }

    private:
        int nbins;
        double min;
        double max;
        double width;
    };

    class Variable
    {
    public:
        Variable(const std::vector<double> &_bindiv_)
            : bindiv(_bindiv_) {}

        int FindBin(double val) const
        {
            return std::upper_bound(bindiv.begin(), bindiv.end(), val) - bindiv.begin();
        }

        int Nbins(OverflowOption option = OverflowOption::WithoutOverflows) const
        {
            if (option == OverflowOption::WithoutOverflows)
                return (int)bindiv.size() - 1;
            else
                return bindiv.size() + 1;
        }

    private:
        std::vector<double> bindiv;
    };
}

#endif
