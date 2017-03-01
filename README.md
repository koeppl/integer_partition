Computational Aspects of Ordered Integer Partition with Upper Bounds
====================================================================

This is a C++14 implementation of the paper 
"Computational Aspects of Ordered Integer Partition with Upper Bounds"
by Roland Glück and Dominik Köppl, SEA 2014 (http://dx.doi.org/10.1007/978-3-642-38527-8_9).
It is an enhanced version supporting parallel execution.

# Goal

The goal is an enhanced version of the classical urn problem: 

> In how many ways can we distribute `z` indistinguishable balls into `n` distinguishable urns?

We restrict each urn to have a specific size.
In other terms:

> The goal is to compute all possible partitions of a given integer `z` as a sum of an ordered sequence of `n` integers,
> with the restriction that each integer of the sequence has an individual upper bound.

# Dependencies

- Command line tools
  - cmake
  - make
  - a C++14 compiler like gcc or clang 
  - The GNU MP Bignum Library (https://gmplib.org)
  - glog
  - gtest
  - gflags
  - celero (optional)

# Installation

This package ships as a library with a test program.
Invoke `cmake` and `make` to compile, `make test` to test the compilation.

# Demo

After compilation, a test program is located at `demo/integer_partition_demo`.
This program can be used to compute the number of distributions of n balls into m urns with constrained capacities `i_1,...,i_m`.
So a call of `./demo/integer_partition_demo 10 100 100` will output the cardinatily of solutions for n = 10 and the capacities to {100,100}.
