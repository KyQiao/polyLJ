/**
 * Mutual exclusion locks using MPI RMA, following the algorithm of Latham et
 * al., 2007.
 *
 * All errors in MPI routines are considered fatal (this may change by define 
 * marco NDEBUG).
 * 
 * copied from https://gist.github.com/aprell/1486197
 * modified into c++ version by me
 */

#include "rma_mutex.hpp"

#include <mpi.h>

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <cstring>

#define NDEBUG
#define MPI_ASSERT(fx)         \
  do {                         \
    assert(fx == MPI_SUCCESS); \
  } while (0)

/**
 * Non-public: Create the MPI derived type corresponding to the req array on
 * the owner rank, minus the element mapped to the origin rank (this is written
 * to the `req_slice_type` member of partially constructed mutex `new`).
 */
static void create_req_buffer_type(rma_mutex_t *new_COMM) {
  int ret, nblock, *len, *dsp;

  new_COMM->req_slice_type = MPI_DATATYPE_NULL;

  len = nullptr;
  dsp = nullptr;

  if (new_COMM->rank == 0 || new_COMM->rank == new_COMM->size - 1)
    nblock = 1;
  else
    nblock = 2;
  len = new int[nblock];
  assert(len != nullptr);
  dsp = new int[nblock];
  assert(dsp != nullptr);

  if (new_COMM->rank == 0) {
    len[0] = new_COMM->size - 1;
    dsp[0] = 1;
  } else if (new_COMM->rank == new_COMM->size - 1) {
    len[0] = new_COMM->size - 1;
    dsp[0] = 0;
  } else {
    len[0] = new_COMM->rank;
    dsp[0] = 0;
    len[1] = new_COMM->size - 1 - new_COMM->rank;
    dsp[1] = new_COMM->rank + 1;
  }

  MPI_ASSERT(MPI_Type_indexed(nblock, len, dsp, MPI_BYTE,
                              &new_COMM->req_slice_type));
  MPI_ASSERT(MPI_Type_commit(&new_COMM->req_slice_type));

  delete[] len;
  delete[] dsp;
}

/**
 * Initialize the RMA mutex `m` owned by rank `owner` of communcator `comm`.
 * This call is collective.
 */
void rma_mutex_init(rma_mutex_t *m, MPI_Comm comm, int owner) {
  int ret;
  MPI_Aint req_size;
  rma_mutex_t new_COMM;

  new_COMM.owner = owner;
  new_COMM.comm = comm;

  new_COMM.req = NULL;
  new_COMM.req_slice_buffer = NULL;

  MPI_ASSERT(MPI_Comm_rank(new_COMM.comm, &new_COMM.rank));
  MPI_ASSERT(MPI_Comm_size(new_COMM.comm, &new_COMM.size));

  new_COMM.req_slice_buffer = new char[new_COMM.size];
  // malloc(new_COMM.size);
  assert(new_COMM.req_slice_buffer != NULL);

  create_req_buffer_type(&new_COMM);

  if (new_COMM.rank == new_COMM.owner) {
    req_size = new_COMM.size;
    MPI_ASSERT(MPI_Alloc_mem(req_size, MPI_INFO_NULL, &new_COMM.req));
    bzero(new_COMM.req, req_size);
  } else {
    req_size = 0;
  }

  MPI_ASSERT(MPI_Win_create(new_COMM.req, req_size, 1, MPI_INFO_NULL,
                            new_COMM.comm, &new_COMM.win));

  /* lookin good - safe to modify argument */
  memcpy(m, &new_COMM, sizeof(rma_mutex_t));
}

/**
 * Free resources associated with the RMA mutex `m`. This call is collective.
 */
void rma_mutex_free(rma_mutex_t *m) {
  int ret;

  MPI_ASSERT(MPI_Win_free(&m->win)); /* collective */

  if (m->rank == m->owner)
    MPI_ASSERT(MPI_Free_mem(m->req));

  MPI_ASSERT(MPI_Type_free(&m->req_slice_type));

  // free(m->req_slice_buffer);
  delete[] m->req_slice_buffer;
}

/**
 * Acquire the lock associated with RMA mutex `m`.
 */
void rma_mutex_lock(rma_mutex_t *m) {
  int i, contested, ret;
  u_int8_t one = 1;

  /* exclusive lock
   * .. get full req_buffer
   * .. set my entry to 1 */
  MPI_ASSERT(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, m->owner, 0, m->win));

  MPI_ASSERT(MPI_Get(m->req_slice_buffer, m->size - 1, MPI_BYTE, m->owner, 0, 1,
                     m->req_slice_type, m->win));
  MPI_ASSERT(MPI_Put(&one, 1, MPI_BYTE, m->owner, m->rank, 1, MPI_BYTE,
                     m->win));

  MPI_ASSERT(MPI_Win_unlock(m->owner, m->win));

  /* lock obtained? */
  contested = 0;
  for (i = 0; i < m->size - 1; i++)
    if (m->req_slice_buffer[i]) {
      contested = 1;
      break;
    }

  /* no, it's contested - wait for it ... */
  if (contested)
    MPI_ASSERT(MPI_Recv(NULL, 0, MPI_BYTE, MPI_ANY_SOURCE, MAGIC_UNLOCK_TAG,
                        m->comm, MPI_STATUS_IGNORE));
}

/**
 * Non-public: Once the `req_slice_buffer` field of mutex `m` has been
 * populated, determine the next rank to receive the lock.
 * Return value: rank on the next process to get the lock (m->rank if no rank is
 * currently waiting).
 */
static int next_rank(rma_mutex_t *m) {
  int i, j, next;
  next = m->rank;
  for (i = 0; i < m->size - 1; i++) {
    j = (m->rank + i) % (m->size - 1);
    if (m->req_slice_buffer[j]) {
      if (j < m->rank)
        next = j;
      else
        next = j + 1;
      break;
    }
  }
  return next;
}

/**
 * Release the lock associated with RMA mutex `m`.
 */
void rma_mutex_unlock(rma_mutex_t *m) {
  int next, ret;
  u_int8_t zero = 0;

  /* exclusive lock
   * .. get new full req_buffer
   * .. set my entry to 0 */
  MPI_ASSERT(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, m->owner, 0, m->win));

  MPI_ASSERT(MPI_Get(m->req_slice_buffer, m->size - 1, MPI_BYTE, m->owner, 0, 1,
                     m->req_slice_type, m->win));
  MPI_ASSERT(MPI_Put(&zero, 1, MPI_BYTE, m->owner, m->rank, 1, MPI_BYTE,
                     m->win));

  MPI_ASSERT(MPI_Win_unlock(m->owner, m->win));

  /* who gets it next? */
  next = next_rank(m);

  /* someone ... give it to them ... */
  if (next != m->rank)
    MPI_ASSERT(MPI_Send(NULL, 0, MPI_BYTE, next, MAGIC_UNLOCK_TAG, m->comm));
}
