#include <time.h>
#include "sigma.h"
#include "matrixvector.h"
#include "SE_basis_func.h"

/* kSE(t, a, b, tau) = k(t, SE_trans(a, b, tau)) */
double kSE(double t, double a, double b, double tau)
{
  double s = SE_trans(a, b, tau);
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

double vRZn(double a, double b, double t, double h, int N, double* f_N, int n)
{
  int j;
  double x = SE_trans_inv(a, b, t);
  double ans = 0;

  for (j = N-1; j > 0; j--) {
    ans += f_N[ j+N] * S( j, h, x);
    ans += f_N[-j+N] * S(-j, h, x);
  }
    ans += f_N[ 0+N] * S( 0, h, x);

    ans += f_N[0]*wa(a, b, t) + f_N[n-1]*wb(a, b, t);

  return ans;
}

double* substitute_tN(double a, double b, double h, int N, int n)
{
  int i;
  double* t_N = AllocVec(n);

  for (i = -N; i <= N; i++) {
    t_N[i+N] = SE_trans(a, b, i*h);
  }

  return t_N;
}

double* substitute_AN(double a, double b, double h, int N, int n, double* t_N)
{
  int i, j;
  double* A_N = AllocVec(n*n); /* Column-major */
  double* V_N = AllocVec(n*n); /* Column-major */
  double* p_N = AllocVec(n);
  double* w_N = AllocVec(n);

  for (j = -N; j <= N; j++) {
    for (i = -N; i <= N; i++) {
      V_N[i+N + (j+N)*n]
        = kSE(t_N[i+N], a, b, j*h)*SE_trans_div(a, b, j*h)*h*(0.5+sigma[i-j]);
    }
  }

  for (i = -N; i <= N; i++) {
    w_N[i+N] = waSE(i*h);
  }

  /* p_N = V_N w_N */
  cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, 1, V_N, n, w_N, 1, 0, p_N, 1);

  /* j = -N */
  j = -N;
    for (i = -N; i <= N; i++) {
      A_N[i+N + (j+N)*n] = w_N[i+N] - p_N[i+N];
    }

  /* -N < j < N */
    for (i = -N+1; i <= N-1; i++) {
      A_N[i+N + (i+N)*n] = 1.0;
    }

  for (j = -N+1; j <= N-1; j++) {
    for (i = -N; i <= N; i++) {
      A_N[i+N + (j+N)*n] -= V_N[i+N + (j+N)*n];
    }
  }

  for (i = -N; i <= N; i++) {
    w_N[i+N] = wbSE(i*h);
  }

  /* p_N = V_N w_N */
  cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, 1, V_N, n, w_N, 1, 0, p_N, 1);

  /* j = N */
  j = N;
    for (i = -N; i <= N; i++) {
      A_N[i+N + (j+N)*n] = w_N[i+N] - p_N[i+N];
    }

  FreeVec(V_N);
  FreeVec(p_N);
  FreeVec(w_N);
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

int main()
{
  double a = 0.0;
  double b = 1.0;
  double d = 3.14;
  double alpha = 0.5;
  double *A_N, *f_N, *t_N;
  int i, n, N, info;
  double h, err, maxerr, t;
  int SAMPLE = 2048;
  double hh = (b - a)/SAMPLE;
  clock_t start, end;
  double time;

  for (N = 5; N <= 300; N += 5) {
    start = clock();

    n = 2*N+1;
    h = sqrt(M_PI*d / (alpha * N));

    t_N = substitute_tN(a, b, h, N, n);
    A_N = substitute_AN(a, b, h, N, n, t_N);
    f_N = substitute_fN(N, n, t_N);

    info = lapack_linsolve(A_N, f_N, n);

    if ( info == 0 ) {

      maxerr = 0;

      for (i = 1; i < SAMPLE; i++) {
        t = a + i*hh;

        err = fabs(u(t) - vRZn(a, b, t, h, N, f_N, n));

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
