#include <time.h>
#include "sigma.h"
#include "matrixvector.h"
#include "DE_basis_func.h"

/* kDE(t, a, b, tau) = k(t, DE_trans(a, b, tau)) */
double kDE(double t, double a, double b, double tau)
{
  double s = DE_trans(a, b, tau);
  return 6 * (sqrt(t) - sqrt(s));
}

double g(double t)
{
  return 1 + sqrt(t) * (1 - 2*t) - t*t;
}

double u(double t)
{
  return 1 + sqrt(t);
}

double vDEn(double a, double b, double t, double h, int N, double* c_N, int n)
{
  int j;
  double x = DE_trans_inv(a, b, t);
  double ans = 0;

  for (j = N; j > 0; j--) {
    ans += c_N[ j+N] * S( j, h, x);
    ans += c_N[-j+N] * S(-j, h, x);
  }
    ans += c_N[ 0+N] * S( 0, h, x);

    ans += c_N[n]*wa(a, b, t) + c_N[n+1]*wb(a, b, t);

  return ans;
}

double* substitute_tN(double a, double b, double h, int N, int n)
{
  int i;
  double* t_N = AllocVec(n);

  for (i = -N; i <= N; i++) {
    t_N[i+N] = DE_trans(a, b, i*h);
  }

  return t_N;
}

double* substitute_AN(double a, double b, double h, int N, int n, double* t_N)
{
  int i, j;
  double* A_N = AllocVec(n*n); /* Column-major */

  for (i = -N; i <= N; i++) {
    A_N[i+N + (i+N)*n] = 1.0;
  }

  for (j = -N; j <= N; j++) {
    for (i = -N; i <= N; i++) {
      A_N[i+N + (j+N)*n]
        -= kDE(t_N[i+N], a, b, j*h)*DE_trans_div(a, b, j*h)*h*(0.5+sigma[i-j]);
    }
  }

  return A_N;
}

double* substitute_fN(int N, int n, double* t_N)
{
  int i;
  double* f_N = AllocVec(n);

  for (i = -N; i <= N; i++) {
    f_N[i+N] = g(t_N[i+N]);
  }

  return f_N;
}

double* calculate_cN(double h, int N, int n, double* f_N)
{
  int j;
  double* c_N = AllocVec(n+2);

  for (j = -N; j <= N; j++) {
    c_N[j+N] = f_N[j+N] - f_N[0]*waDE(j*h) - f_N[n-1]*wbDE(j*h);
  }

  c_N[n]  = f_N[0];
  c_N[n+1]= f_N[n-1];

  return c_N;
}

int main()
{
  double a = 0.0;
  double b = 1.0;
  double d = 1.57;
  double alpha = 0.5;
  double *A_N, *f_N, *t_N, *c_N;
  int i, n, N, info;
  double h, err, maxerr, t;
  int SAMPLE = 2048;
  double hh = (b - a)/SAMPLE;
  clock_t start, end;
  double time;

  for (N = 5; N <= 80; N += 5) {
    start = clock();

    n = 2*N+1;
    h = log(2 * d * N / alpha) / N;

    t_N = substitute_tN(a, b, h, N, n);
    A_N = substitute_AN(a, b, h, N, n, t_N);
    f_N = substitute_fN(N, n, t_N);

    info = lapack_linsolve(A_N, f_N, n);

    if ( info == 0 ) {

      maxerr = 0;
      c_N = calculate_cN(h, N, n, f_N);

      for (i = 1; i < SAMPLE; i++) {
        t = a + i*hh;

        err = fabs(u(t) - vDEn(a, b, t, h, N, c_N, n));

        maxerr = fmax(err, maxerr);
      }

      printf("%d\t%e\t", N, maxerr);
    } else fprintf(stderr, "error in lapack_linsolve!\n");

    FreeVec(f_N);
    FreeVec(A_N);
    FreeVec(t_N);

    end = clock();
    time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("%e\n", time);
  }

  return EXIT_SUCCESS;
}
