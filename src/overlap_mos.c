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

  printf("Reading %s.\n", file_name);
  rc = qmckl_trexio_read(context, file_name, 255);
  if (rc != QMCKL_SUCCESS) {
    printf("%s\n", qmckl_string_of_error(rc));
  }
  assert (rc == QMCKL_SUCCESS);

  const int64_t nx = 40;
  const int64_t ny = 40;
  const int64_t nz = 40;
  const int64_t np = nx*ny*nz;

  const int64_t elec_num = np;

  rc = qmckl_set_electron_walk_num(context, 1L);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_electron_num(context, np/2, np/2);
  assert (rc == QMCKL_SUCCESS);

  int64_t nucl_num;
  rc = qmckl_get_nucleus_num(context, &nucl_num);
  assert (rc == QMCKL_SUCCESS);

  double* nucl_coord = malloc(nucl_num * 3 * sizeof(double));
  rc = qmckl_get_nucleus_coord(context, 'T', nucl_coord, nucl_num*3);
  assert (rc == QMCKL_SUCCESS);

  double rmin[3] = { nucl_coord[0], nucl_coord[nucl_num], nucl_coord[2*nucl_num] };
  double rmax[3] = { nucl_coord[0], nucl_coord[nucl_num], nucl_coord[2*nucl_num] };

  for (int64_t i=0 ; i<nucl_num ; ++i) {
    rmin[0] = rmin[0] < nucl_coord[i] ? rmin[0] : nucl_coord[i];
    rmin[1] = rmin[1] < nucl_coord[i+nucl_num] ? rmin[1] : nucl_coord[i+nucl_num];
    rmin[2] = rmin[2] < nucl_coord[i+2*nucl_num] ? rmin[2] : nucl_coord[i+2*nucl_num];
    rmax[0] = rmax[0] > nucl_coord[i] ? rmax[0] : nucl_coord[i];
    rmax[1] = rmax[1] > nucl_coord[i+nucl_num] ? rmax[1] : nucl_coord[i+nucl_num];
    rmax[2] = rmax[2] > nucl_coord[i+2*nucl_num] ? rmax[2] : nucl_coord[i+2*nucl_num];
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

  const double dx = (rmax[0] - rmin[0]) / (nx+1);
  const double dy = (rmax[1] - rmin[1]) / (ny+1);
  const double dz = (rmax[2] - rmin[2]) / (nz+1);

  printf("Computing grid x,y,z\n"); fflush(stdout);
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

  int64_t mo_num;
  rc = qmckl_get_mo_basis_mo_num(context, &mo_num);
  assert (rc == QMCKL_SUCCESS);

  const int64_t size_max = 5*elec_num*mo_num;
  double * mo_vgl = malloc (size_max * sizeof(double));
  assert (mo_vgl != NULL);


  printf("%e %e %e  %e\n", coord[0], coord[1], coord[2], mo_vgl[0]);
  printf("Setting electron coordinates\n"); fflush(stdout);
  rc = qmckl_set_electron_coord(context, 'N', coord, elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  printf("DGEMM MOs\n"); fflush(stdout);
  rc = qmckl_get_mo_basis_mo_vgl(context, mo_vgl, 5*elec_num*mo_num);
  assert (rc == QMCKL_SUCCESS);

  double * overlap = malloc(mo_num * mo_num * sizeof(double));

  printf("DGEMM grid\n"); fflush(stdout);
  rc = qmckl_dgemm(context, 'N', 'T', mo_num, mo_num, elec_num, dx*dy*dz,
                   mo_vgl, mo_num, mo_vgl, mo_num, 0.0, overlap, mo_num);

  for (int j=0 ; j<mo_num ; ++j) {
    for (int i=0 ; i<mo_num ; ++i) {
      printf("%d %d %f\n", i, j, overlap[i + j*mo_num]);
    }
  }

}
