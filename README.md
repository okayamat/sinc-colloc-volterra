# sinc-colloc-volterra
Numerical solvers for Volterra integral equations of the second kind by Sinc-collocation methods

## Overview
These programs solve two examples of Volterra integral equations of the
second kind, conducted in [2].

Those problems are solved by means of the following 5 methods:
* SE-Sinc-Nyström method [1]
* DE-Sinc-Nyström method [1]
* SE-Sinc-collocation method (by Rashidinia and Zarebnia) [3]
* Another SE-Sinc-collocation method (by Stenger) [4]
* DE-Sinc-collocation method [2]

The name of the program denotes the method and example number. For
example, SE_nyst_ex1.c denotes the SE-Sinc-Nyström method for Example 1,
and RZ_coll_ex2.c denotes the SE-Sinc-collocation method by Rashidinia
and Zarebnia for Example 2. LAPACK in Apple's Accelerate framework is
used for computation of the system of linear equations. If you want to
use another LAPACK library, modify make files according to your
installation.

Each program solves those problems increasing N as N = 5, 10, 15, 20, ...,
and outputs N, maximum error over the target interval, and computation
time.

## Results
Outputs by those programs are stored in data/ directory, with .dat extension.
Gnuplot programs for creating graphs are also stored in the directory.

computation environment:

OS: Mac OS Big Sur  
CPU: 1.7 GHz Intel Core i7  
Memory: 8 GB 1600 MHz DDR3  
Compiler: Apple clang version 13.0.0  
Library: LAPACK (Apple's Accelerate framework)

## References
1. M. Muhammad, A. Nurmuhammad, M. Mori, and M. Sugihara: Numerical solution
 of integral equations by means of the Sinc collocation method based on the
 double exponential transformation, J. Comput. Appl. Math., Vol. 177 (2005),
 pp. 269--286.
2. T. Okayama: On the relation between two Sinc-collocation methods for
 Volterra integral equations of the second kind and further improvement,
 arXiv.
3. J. Rashidinia and M. Zarebnia: Solution of a Volterra integral equation by
 the sinc-collocation method, J. Comput. Appl. Math., Vol. 206 (2007),
 pp. 801--813.
4. F. Stenger: Numerical Methods Based on Sinc and Analytic Functions,
 Springer-Verlag, New York, 1993.
