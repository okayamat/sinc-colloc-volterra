# sinc-colloc-volterra

Numerical solvers for Volterra integral equations of the second kind using
Sinc-Nyström and Sinc-collocation methods.

## Overview

The programs reproduce the numerical experiments for the following three
examples considered in [2]:

1. a problem with a smooth solution from [3];
2. a problem with an endpoint singularity from [5];
3. a problem with an oscillatory kernel from [6].

Each example is solved using the following five methods:

- SE-Sinc-Nyström method [1]
- DE-Sinc-Nyström method [1]
- SE-Sinc-collocation method of Rashidinia and Zarebnia [3]
- SE-Sinc-collocation method of Stenger [4]
- DE-Sinc-collocation method proposed in [2]

There are therefore 15 main programs. Their names identify the method and
example number:

| Prefix | Method |
| --- | --- |
| `SE_nyst` | SE-Sinc-Nyström |
| `DE_nyst` | DE-Sinc-Nyström |
| `RZ_coll` | SE-Sinc-collocation of Rashidinia and Zarebnia |
| `SE_coll` | SE-Sinc-collocation of Stenger |
| `DE_coll` | DE-Sinc-collocation proposed in [2] |

For example, `SE_nyst_ex1.c` implements the SE-Sinc-Nyström method for
Example 1, and `RZ_coll_ex3.c` implements the Rashidinia--Zarebnia method for
Example 3.

## Requirements

The supplied makefiles are configured for the following environment:

- macOS
- C compiler (`gcc` or Apple clang)
- `make`
- LAPACK from Apple's Accelerate framework
- Gnuplot for regenerating the figures

The Cephes Math Library is used to evaluate the sine integral. The required
files (`mconf.h`, `polevl.c`, and `sici.c`) are included in this repository.

To use another LAPACK implementation, modify the compiler and linker settings
in the `.mak` files according to the local installation. The LAPACK interface
used in `matrixvector.h` may also need to be adjusted on non-macOS systems.

## Building and running

Specify the makefile with the `-f` option. For example, to build and run the
SE-Sinc-collocation method of Stenger for Example 1, use

```sh
make -f SE_coll_ex1.mak
./SE_coll_ex1
```

To save the output in the format used by the files in `data/`, use

```sh
./SE_coll_ex1 > data/SE_coll_ex1.dat
```

To remove the executable and object files, use

```sh
make -f SE_coll_ex1.mak clean
```

Replace `SE_coll_ex1` with the desired method and example name. Each program
increases the discretization parameter as `N = 5, 10, 15, 20, ...` and writes
three columns to standard output:

1. `N`;
2. the maximum error evaluated at 2047 equally spaced interior test points;
3. the computation time in seconds.

The reported time includes construction of the coefficient matrix and
right-hand side, solution of the linear system, and evaluation of the
approximate solution at the test points.

## Numerical data and figures

The `data/` directory contains the output data (`.dat`), Gnuplot scripts
(`.plt`), and generated figures (`.eps`). For example,

```sh
cd data
gnuplot example1_N.plt
gnuplot example1_t.plt
```

generates the two figures for Example 1. The files `example2_*.plt` and
`example3_*.plt` generate the corresponding figures for Examples 2 and 3.

The three columns of each `.dat` file are the values of `N`, the maximum
error, and the computation time. Tables 1--3 of [2] list selected rows from
these files.

## Reproducing the frequency comparison

The main results for Example 3 use `omega = 20.0`. To reproduce the frequency
comparison in Table 4 of [2], edit

```c
const double omega = 20.0;
```

in each `*_ex3.c` program, replacing `20.0` successively with `5.0`, `10.0`,
`20.0`, and `40.0`. Recompile and rerun each program after every change. The
row with `N = 100` gives the value reported in Table 4. For example,

```sh
make -f DE_coll_ex3.mak clean
make -f DE_coll_ex3.mak
./DE_coll_ex3
```

## Computation environment used in [2]

- OS: macOS Big Sur
- CPU: 1.7 GHz Intel Core i7
- Memory: 8 GB 1600 MHz DDR3
- Compiler: Apple clang version 13.0.0
- Linear algebra library: LAPACK from Apple's Accelerate framework
- Floating-point arithmetic: double precision

## References

1. M. Muhammad, A. Nurmuhammad, M. Mori, and M. Sugihara, "Numerical
   solution of integral equations by means of the Sinc collocation method
   based on the double exponential transformation," *J. Comput. Appl. Math.*,
   177 (2005), 269--286.
2. T. Okayama, "Relation between two Sinc-collocation methods for Volterra
   integral equations of the second kind and further improvement,"
   arXiv:2502.20221 [math.NA].
3. J. Rashidinia and M. Zarebnia, "Solution of a Volterra integral equation by
   the Sinc-collocation method," *J. Comput. Appl. Math.*, 206 (2007),
   801--813.
4. F. Stenger, *Numerical Methods Based on Sinc and Analytic Functions*,
   Springer-Verlag, New York, 1993.
5. A. D. Polyanin and A. V. Manzhirov, *Handbook of Integral Equations*,
   2nd ed., Chapman & Hall/CRC, Boca Raton, FL, 2008.
6. L. Fermo and C. van der Mee, "Volterra integral equations with highly
   oscillatory kernels: a new numerical method with applications,"
   *Electron. Trans. Numer. Anal.*, 54 (2021), 333--354.
