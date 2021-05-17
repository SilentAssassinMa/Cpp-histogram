# Cpp-histogram
This is a light-weight histogram class on c++ platform which provides similar usage as THX in CERN ROOT. This project aims to provide a alternative of THXD in ROOT due to its troublesome in concurrency coding. Basic usage : 

```c++
CppHistogram histogram(Axis::Regular(10, 1, 10), Axis::Variable({0, 1, 4, 9, 16}));
histogram.Fill(3, 14);
histogram.GetBinContent(1, 2);
histogram.GetAxis<0>().FindBin(10);
histogram.GetStats().GetMean(0);
```

`CppHistogram` uses `CppCounter` as a counter for calculating mean and standard variation like CERN ROOT does.

This project is inspired by boost::histogram but it's much more lighter without any outside dependence. However compiler with c++ 17 support is needed since I used `if constexpr`(I know some people are still using gcc 4.8, but...)

# LICENSE

This project follows a DO-WHAT-THE-FUCK-YOU-WANT-TO license. Feed-backs are welcome and I will add more methods in `CppHistogram` in the future, like converting to THXD in ROOT and integrate over some axes.