#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <qmckl.h>
#include <trexio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <unistd.h>
#include <sys/time.h>

#include <time.h>

#define ITERMAX 4000000000

int main(int argc, char** argv)
{
  long start, end;
  struct timeval timecheck;

  qmckl_context context;
  context = qmckl_context_create();
  qmckl_exit_code rc;

  double* a = malloc(100*100*sizeof(double));
  double* b = malloc(100*100*sizeof(double));
  double det;
  for (int n=1 ; n<=10 ; ++n) {
    gettimeofday(&timecheck, NULL);
    start = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;
    
    const size_t niter = ITERMAX/(n*n*n) ;
    for (size_t k=0 ; k<niter ; ++k) {
      rc = qmckl_adjugate(context, n, a, n, b, n, &det);
    }

    gettimeofday(&timecheck, NULL);
    end = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

    printf("%3d %18.10e\n", n, (double) (end-start) / (double) niter);
  }


  rc = qmckl_context_destroy(context);

}
