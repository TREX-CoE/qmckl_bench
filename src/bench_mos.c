#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "qmckl.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <unistd.h>
#include <sys/time.h>
#include "Alz_small.h"
#include "h2o-sto3g.h"

#include <time.h>

#define ITERMAX 10

int main(int argc, char** argv)
{
  long start, end;
  struct timeval timecheck;

  qmckl_context context;
  context = qmckl_context_create();
  qmckl_exit_code rc;

  const char* file_name  = alz_small_file_name;
  const int64_t elec_num = alz_small_elec_num;
  const int64_t walk_num = alz_small_walk_num;

  printf("Reading %s.\n", file_name);
  rc = qmckl_trexio_read(context, file_name, 255);
  if (rc != QMCKL_SUCCESS) {
    printf("%s\n", qmckl_string_of_error(rc));
  }
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_electron_walk_num(context, walk_num);
  assert (rc == QMCKL_SUCCESS);

  int64_t mo_num;
  rc = qmckl_get_mo_basis_mo_num(context, &mo_num);
  assert (rc == QMCKL_SUCCESS);

  const int64_t size_max = 5*walk_num*elec_num*mo_num;
  double * mo_vgl = malloc (size_max * sizeof(double));
  assert (mo_vgl != NULL);

  rc = qmckl_set_electron_coord(context, 'T', alz_small_elec_coord, walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_mo_basis_mo_vgl(context, mo_vgl, size_max);
  assert (rc == QMCKL_SUCCESS);

  gettimeofday(&timecheck, NULL);
  start = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

  for (int i=0 ; i<ITERMAX ; ++i) {
    rc = qmckl_get_mo_basis_mo_vgl_inplace(context, mo_vgl, size_max);
  }
  gettimeofday(&timecheck, NULL);
  end = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

  printf("Time for the calculation of 1 step (ms): %10.1f\n", (double) (end-start) / (double) ITERMAX);


}
