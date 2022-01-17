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


int main(int argc, char** argv)
{
  qmckl_context context;
  context = qmckl_context_create();
  qmckl_exit_code rc;

  const char* file_name  = alz_small_file_name;
  const int64_t walk_num = alz_small_walk_num;
  const int64_t elec_num = alz_small_elec_num;

  printf("Reading %s.\n", alz_small_file_name);
  rc = qmckl_trexio_read(context, file_name);
  if (rc != QMCKL_SUCCESS) {
    printf("%s\n", qmckl_string_of_error(rc));
  }
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_electron_walk_num(context, walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_electron_coord(context, 'N', alz_small_elec_coord, walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);


}
