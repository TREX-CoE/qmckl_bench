#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "qmckl.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include "Alz_small.h"
#include "h2o-sto3g.h"

#include <time.h>

#define ITERMAX 10000

int main(int argc, char** argv)
{
  clock_t start, end;
  double cpu_time_used;

  qmckl_context context;
  context = qmckl_context_create();
  qmckl_exit_code rc;

  const char* file_name  = alz_small_file_name;
  const int64_t elec_num = alz_small_elec_num;
/*
  const char* file_name  = h2o_sto_file_name;
  const int64_t elec_num = h2o_sto_elec_num;
*/
  const int64_t walk_num = 1; // alz_small_walk_num;

  printf("Reading %s.\n", file_name);
  rc = qmckl_trexio_read(context, file_name);
  if (rc != QMCKL_SUCCESS) {
    printf("%s\n", qmckl_string_of_error(rc));
  }
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_electron_walk_num(context, walk_num);
  assert (rc == QMCKL_SUCCESS);

  int64_t ao_num;
  rc = qmckl_get_ao_basis_ao_num(context, &ao_num);
  assert (rc == QMCKL_SUCCESS);

  const int64_t size_max = 5*elec_num*ao_num;
  double * ao_vgl = malloc (size_max * sizeof(double));
  assert (ao_vgl != NULL);

  rc = qmckl_set_electron_coord(context, 'N', h2o_sto_elec_coord, walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_ao_basis_ao_vgl(context, ao_vgl, size_max);
  assert (rc == QMCKL_SUCCESS);


  double * x;
/*
  x = alz_small_elec_coord;
  for (int l=0 ; l<3; ++l) {
    for (int i=0 ; i<elec_num; ++i) {
      printf("%f ", *x);
      ++x;
    }
  }
  printf("\n ");

  x = ao_vgl;
  for (int l=0 ; l<5; ++l) {
    for (int i=0 ; i<elec_num; ++i) {
      for (int j=0 ; j<ao_num; ++j) {
        printf("%3d %3d %3d %20.10e %20.10e\n", j+1,i+1,l,*x, h2o_sto_elec_coord[i]);
        ++x;
      }
    }
  }
*/

  start = clock();
  for (int i=0 ; i<ITERMAX ; ++i) {
    rc = qmckl_set_electron_coord(context, 'T', h2o_sto_elec_coord, walk_num*elec_num*3);
    rc = qmckl_get_ao_basis_ao_vgl(context, ao_vgl, size_max);
  }
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("Time for the calculation of 1 step (ms): %f\n", 1000. * cpu_time_used / ITERMAX);


}
