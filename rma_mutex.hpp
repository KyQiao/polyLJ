#ifndef _RMA_MUTEX_H
#define _RMA_MUTEX_H

#include <mpi.h>

#define MAGIC_UNLOCK_TAG 0xBEEF

struct rma_mutex
{
  int owner, rank, size;
  char *req, *req_slice_buffer;
  MPI_Win win;
  MPI_Comm comm;
  MPI_Datatype req_slice_type;
};
typedef struct rma_mutex rma_mutex_t;

void rma_mutex_init(rma_mutex_t *, MPI_Comm, int);

void rma_mutex_free(rma_mutex_t *);

void rma_mutex_lock(rma_mutex_t *);

void rma_mutex_unlock(rma_mutex_t *);

#endif
