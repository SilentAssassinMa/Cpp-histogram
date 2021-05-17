// Multi dimensional counter
#ifndef CPPCOUNTER_HEADER
#define CPPCOUNTER_HEADER
#include <iostream>
#include <tuple>
#include <vector>
#include <algorithm>
#include <functional>

template <size_t Dimension>
class CppCounter
{
public:
    CppCounter() : SumX(Dimension, 0), SumX2(Dimension, std::vector<double>(Dimension, 0)){};

    template <typename... FillPar>
    inline void Add(FillPar &&...par)
    {
        constexpr std::size_t npar = sizeof...(FillPar);
        static_assert(npar == Dimension || npar == Dimension + 1, "must be match");
        if constexpr (npar == Dimension)
        {
            ++SumW;
            ++SumW2;
            Add_impl(std::forward_as_tuple(par...));
        }
        if constexpr (npar == Dimension + 1)
            AddWithWeight_impl(std::forward_as_tuple(par...));
    }

    // other functions are easy to write
    inline double GetSum(int index) const
    {
        if (index >= 0 && index < SumX.size())
            return SumX[index];
        else
            return 0;
    }

    inline double GetSum2(int index) const
    {
        if (index >= 0 && index < SumX.size())
            return SumX2[index][index];
        else
            return 0;
    }

    inline double GetSumW() const { return SumW; }
    inline double GetSumW2() const { return SumW2; }

    inline double GetMean(int index) const
    {
        if (index >= 0 && index < SumX.size())
            return SumW != 0 ? SumX[index] / SumW : 0;
        else
            return 0;
    }

    inline double GetSumStatUnc(int index) const
    {
        if (index >= 0 && index < SumX.size())
            return GetMean() * sqrt(SumW2);
        else
            return 0;
    }

    inline double GetMeanStatUnc(int index) const
    {
        return SumW != 0 ? GetStdVar(index) * sqrt(SumW2) / SumW : 0;
    }

    inline double GetVar(int index) const
    {
        if (index >= 0 && index < SumX.size())
            return SumW != 0 ? (SumX2[index][index] / SumW - (SumX[index] / SumW) * (SumX[index] / SumW)) : 0;
        else
            return 0;
    }

    inline double GetStdVar(int index) const
    {
        return sqrt(fabs(GetVar()));
    }

    inline double GetCovariance(int i1, int i2) const
    {
        if (i1 >= 0 && i1 < SumX.size() && i2 >= 0 && i2 < SumX.size())
            return SumW != 0 ? (SumX2[i1][i2] / SumW - (SumX[i1] / SumW) * (SumX[i2] / SumW)) : 0;
        else
            return 0;
    }

    inline double GetCorrelation(int i1, int i2) const
    {
        double stdvar1 = GetStdVar(i1);
        double stdvar2 = GetStdVar(i2);
        if (stdvar1 != 0 && stdvar2 != 0)
            return GetCovariance(i1, i2) / (stdvar1 * stdvar2);
        else
            return 0;
    }

private:
    template <
        typename TTuple,
        size_t Index = 0,
        size_t Size = std::tuple_size<std::remove_reference_t<TTuple>>::value>
    inline void Add_impl(TTuple &&value)
    {
        if constexpr (Index < Size)
        {
            SumX[Index] += std::get<Index>(value);
            SumX2[Index][Index] += std::get<Index>(value) * std::get<Index>(value);
            if constexpr (Index > 0)
                AddX2_impl<TTuple, Index, Index - 1>(std::forward<TTuple>(value));
        }

        if constexpr (Index < Size - 1)
            Add_impl<TTuple, Index + 1, Size>(std::forward<TTuple>(value));
    }

    template <
        typename TTuple,
        size_t Index = 0,
        size_t Size = std::tuple_size<std::remove_reference_t<TTuple>>::value>
    inline void AddWithWeight_impl(TTuple &&value)
    {
        SumW += std::get<Size - 1>(value);
        SumW2 += std::get<Size - 1>(value) * std::get<Size - 1>(value);
        Add_impl<TTuple, Index, Size - 1>(std::forward<TTuple>(value));
    }

    template <typename TTuple, size_t Index1, size_t Index2>
    inline void AddX2_impl(TTuple &&value)
    {
        if constexpr (Index2 >= 0)
            SumX2[Index2][Index1] += std::get<Index1>(value) * std::get<Index2>(value);
        if constexpr (Index2 > 0)
            AddX2_impl<TTuple, Index1, Index2 - 1>(std::forward<TTuple>(value));
    }

    double SumW = 0;
    double SumW2 = 0;
    std::vector<double> SumX;
    std::vector<std::vector<double>> SumX2;
};

#endif
