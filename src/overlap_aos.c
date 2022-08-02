#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "qmckl.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <time.h>

#define ITERMAX 10000

int main(int argc, char** argv)
{
  qmckl_context context;
  context = qmckl_context_create();
  qmckl_exit_code rc;

  if (argc < 2) {
    fprintf(stderr,"Syntax: %s FILE.    (FILEs are data/*.h5)\n", argv[0]);
    exit(-1);
  }
  char* file_name = argv[1];

  printf("Reading %s.\n", file_name);
  rc = qmckl_trexio_read(context, file_name, 255);
  if (rc != QMCKL_SUCCESS) {
    printf("%s\n", qmckl_string_of_error(rc));
  }
  assert (rc == QMCKL_SUCCESS);

  const int64_t nx = 80;
  const int64_t ny = 80;
  const int64_t nz = 80;
  const int64_t np = nx*ny*nz;

  const int64_t point_num = np;

  int64_t nucl_num;
  rc = qmckl_get_nucleus_num(context, &nucl_num);
  assert (rc == QMCKL_SUCCESS);

  double* nucl_coord = malloc(nucl_num * 3 * sizeof(double));
  rc = qmckl_get_nucleus_coord(context, 'N', nucl_coord, nucl_num*3);
  assert (rc == QMCKL_SUCCESS);

  printf("Nuclei:\n");
  for (int i=0 ; i<3*nucl_num ; i+=3) {
    printf("%f %f %f\n", nucl_coord[i], nucl_coord[i+1], nucl_coord[i+2]);
  }
  printf("\n");
  double rmin[3] = { nucl_coord[0], nucl_coord[1], nucl_coord[2] };
  double rmax[3] = { nucl_coord[0], nucl_coord[1], nucl_coord[2] };

  for (int64_t i=0 ; i<3*nucl_num ; i+=3) {
    rmin[0] = rmin[0] < nucl_coord[i] ? rmin[0] : nucl_coord[i];
    rmin[1] = rmin[1] < nucl_coord[i+1] ? rmin[1] : nucl_coord[i+1];
    rmin[2] = rmin[2] < nucl_coord[i+2] ? rmin[2] : nucl_coord[i+2];
    rmax[0] = rmax[0] > nucl_coord[i] ? rmax[0] : nucl_coord[i];
    rmax[1] = rmax[1] > nucl_coord[i+1] ? rmax[1] : nucl_coord[i+1];
    rmax[2] = rmax[2] > nucl_coord[i+2] ? rmax[2] : nucl_coord[i+2];
  }
  rmin[0] -= 5.;
  rmin[1] -= 5.;
  rmin[2] -= 5.;
  rmax[0] += 5.;
  rmax[1] += 5.;
  rmax[2] += 5.;

  printf("Box:\n");
  printf("%f %f %f \n", rmin[0], rmin[1], rmin[2]);
  printf("%f %f %f \n", rmax[0], rmax[1], rmax[2]);
  printf("\n");

  const double dx = (rmax[0] - rmin[0]) / (nx-1);
  const double dy = (rmax[1] - rmin[1]) / (ny-1);
  const double dz = (rmax[2] - rmin[2]) / (nz-1);

  printf("Integration steps:\n");
  printf("%f %f %f \n", dx, dy, dz);
  printf("\n");

  double * coord = malloc ( np * 3 * sizeof(double));
  double x, y, z;
  double *p = coord;
  for (int64_t iz=0 ; iz<nz ; iz++) {
    z = rmin[2] + iz * dz;
    for (int64_t iy=0 ; iy<ny ; iy++) {
      y = rmin[1] + iy * dy;
      for (int64_t ix=0 ; ix<nx ; ix++) {
        x = rmin[0] + ix * dx;
        *p = x; ++p;
        *p = y; ++p;
        *p = z; ++p;
      }
    }
  }

  int64_t ao_num;
  rc = qmckl_get_ao_basis_ao_num(context, &ao_num);
  assert (rc == QMCKL_SUCCESS);

  const int64_t size_max = 5*point_num*ao_num;
  double * ao_vgl = malloc (size_max * sizeof(double));
  assert (ao_vgl != NULL);


  rc = qmckl_set_point(context, 'N', point_num, coord, 3*point_num);
  assert (rc == QMCKL_SUCCESS);

  printf("DGEMM AOs\n"); fflush(stdout);
  rc = qmckl_get_ao_basis_ao_vgl(context, ao_vgl, 5*point_num*ao_num);
  assert (rc == QMCKL_SUCCESS);

  double * overlap = malloc(ao_num * ao_num * sizeof(double));

  printf("DGEMM grid\n"); fflush(stdout);
  rc = qmckl_dgemm(context, 'N', 'T', ao_num, ao_num, point_num, dx*dy*dz,
                   ao_vgl, ao_num*5, ao_vgl, ao_num*5, 0.0, overlap, ao_num);

  for (int j=0 ; j<ao_num ; ++j) {
    for (int i=0 ; i<ao_num ; ++i) {
      printf("%d %d %f\n", i, j, overlap[i + j*ao_num]);
    }
  }

}
