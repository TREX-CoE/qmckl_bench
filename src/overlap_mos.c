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
  int nx, ny, nz;
  nx = 80;
  if (argc > 2)
    sscanf(argv[2], "%d", &nx);
  ny = nx;
  nz = nx;
  if (argc > 3)
    sscanf(argv[3], "%d", &ny);
  if (argc > 4)
    sscanf(argv[4], "%d", &nz);

  printf("nx = %d  ny = %d  nz = %d\n", nx, ny, nz);

  printf("Reading %s.\n", file_name);
  rc = qmckl_trexio_read(context, file_name, 255);
  if (rc != QMCKL_SUCCESS) {
    printf("%s\n", qmckl_string_of_error(rc));
  }
  assert (rc == QMCKL_SUCCESS);

  const int64_t point_num = nx*ny;

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

  double * coord = malloc ( point_num * 3 * sizeof(double));
  double x, y, z;

  int64_t mo_num;
  rc = qmckl_get_mo_basis_mo_num(context, &mo_num);
  assert (rc == QMCKL_SUCCESS);

  const int64_t size_max = 5*point_num*mo_num;
  double * mo_vgl = malloc (size_max * sizeof(double));
  assert (mo_vgl != NULL);

  double * overlap = calloc(mo_num * mo_num, sizeof(double));
  assert (overlap != NULL);

  for (int64_t iz=0 ; iz<nz ; iz++) {
    double *p = coord;
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

    rc = qmckl_set_point(context, 'N', point_num, coord, 3*point_num);
    assert (rc == QMCKL_SUCCESS);

    printf("."); fflush(stdout);
    rc = qmckl_get_mo_basis_mo_vgl(context, mo_vgl, 5*point_num*mo_num);
    assert (rc == QMCKL_SUCCESS);

    rc = qmckl_dgemm(context, 'N', 'T', mo_num, mo_num, point_num, dx*dy*dz,
                    mo_vgl, mo_num*5, mo_vgl, mo_num*5, 1.0, overlap, mo_num);
  }
  printf("\n\n"); 

  double error = 0.;
  for (int j=0 ; j<mo_num ; ++j) {
    for (int i=0 ; i<j ; ++i) {
      const double x = overlap[i + j*mo_num];
      error += x*x;
    }
    const double x = overlap[j + j*mo_num] - 1.;
    error += x*x;
    for (int i=j+1 ; i<mo_num ; ++i) {
      const double x = overlap[i + j*mo_num];
      error += x*x;
    }
    for (int i=0 ; i<mo_num ; ++i) {
      printf("%d %d %f\n", i, j, overlap[i + j*mo_num]);
    }
  }
  error = sqrt(error);
  printf("\nError: %e\n", error);



}
