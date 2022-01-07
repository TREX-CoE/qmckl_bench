#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "qmckl.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>


#define WALK_NUM  1000
int main(int argc, char** argv)
{
  const int64_t walk_num = WALK_NUM;

  if (argc < 2){
    fprintf(stderr,"syntax: %s [FILE]\n", argv[0]);
    return -1;
  }
  char* file_name = argv[1];

  qmckl_context context;
  context = qmckl_context_create();
  qmckl_exit_code rc;

  rc = qmckl_set_electron_walk_num(context, walk_num);
  assert (rc == QMCKL_SUCCESS);


  rc = qmckl_trexio_read(context, file_name);
  assert (rc == QMCKL_SUCCESS);

}
