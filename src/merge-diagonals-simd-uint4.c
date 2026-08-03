#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "merge-diagonals-simd-uint4.h"
#include "assert.h"
#include "mem.h"
#include "popcount.h"		/* For clz_table */
#include "sedgesort.h"
#include "merge-uint4.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */

#include "simd.h"

#if 0
#if defined(HAVE_SSE4_1)
#include <smmintrin.h>
#endif
#if defined(HAVE_AVX2)
#include <immintrin.h>
#endif
#if defined(HAVE_AVX512)
#include <immintrin.h>
#endif
#endif


/* #define PYRAMID_SIZE 4 */

#define CUTOFF 1000
#define PYRAMID_SIZE 32


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* pyramid_merge */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


#ifdef CHECK_INPUTS
static void
check_ascending (unsigned int *diagonals, int ndiagonals) {
  int i;

  for (i = 1; i < ndiagonals; i++) {
    if (diagonals[i] < diagonals[i-1]) {
      abort();
    }
  }
  return;
}
#endif


#define PARENT(i) (i >> 1)
#define LEFT(i) (i << 1)
#define RIGHT(i) ((i << 1) | 1)

/* Each heap elt is aligned */
static int
pyramid_merge (unsigned int **_heap, int nstreams, int heapsize, int *nelts,
	       int pyramid_start, int pyramid_end) {
  int nodei;
#ifdef DEBUG
  int i;
#endif

  while (pyramid_end > pyramid_start) {
    debug(printf("Merging level: %d..%d for heapsize %d\n",pyramid_start,pyramid_end,heapsize));

    if (pyramid_end > heapsize) {
      nodei = heapsize;
    } else {
      nodei = pyramid_end;
    }

    while (nodei >= pyramid_start) {
      debug2(printf("Merging nodes %d (%d elts) and %d (%d elts) => %d\n",
		    nodei-1,nelts[nodei-1],nodei,nelts[nodei],PARENT(nodei)));
      _heap[PARENT(nodei)] = Merge_uint4(/*dest*/NULL,_heap[nodei-1],_heap[nodei],nelts[nodei-1],nelts[nodei]);
      CHECK_ALIGN(_heap[PARENT(nodei)]);
      nelts[PARENT(nodei)] = nelts[nodei-1] + nelts[nodei];
      debug2(printf("Created list %p of length %d at node %d\n",
		    _heap[PARENT(nodei)],nelts[PARENT(nodei)],PARENT(nodei)));

#ifdef DEBUG
      for (i = 0; i < nelts[PARENT(nodei)]; i++) {
	printf("%u\n",_heap[PARENT(nodei)][i]);
      }
#endif

      /* Don't free original lists (when nodei >= nstreams) */
      debug(printf("Freeing nodes %d and %d\n",nodei-1,nodei));
      if (nodei < nstreams) {
	FREE_ALIGN(_heap[nodei]);
      }
      if (nodei-1 < nstreams) {
	FREE_ALIGN(_heap[nodei-1]);
      }
      nodei -= 2;
    }

    pyramid_end = PARENT(pyramid_end);
    pyramid_start = PARENT(pyramid_start);
  }

  debug(printf("Returning ancestor %d\n\n",pyramid_start));
  return pyramid_start;
}


/* Fills mergeinfo->heap and mergeinfo->nelts */
static int
make_diagonals_heap (Mergeinfo_uint4_T mergeinfo, UINT4 **stream_array,
		     int *streamsize_array, int *diagterm_array, int nstreams) {
  int ncopied, *nelts;
  UINT4 **_heap, **_combined, *out, diagterm;
  int heapsize, heapi, ncombined;
  int *order, *totals, total, n, i, j, k, l;

  /* Put smallest streams at the end of the heap, to increase speed */

  /* Note: Sedgesort_order_int overwrites at array[n], so we need to copy */
  memcpy(mergeinfo->streamsize_copy,streamsize_array,nstreams*sizeof(int));
  order = Sedgesort_order_int(mergeinfo->streamsize_copy,nstreams);

  assert(nstreams < mergeinfo->querylength);
  _combined = mergeinfo->combined;
  totals = mergeinfo->totals;
  ncombined = 0;

  /* Combine the smallest streams, up to CUTOFF, and just use Sedgesort */
  i = 0;
  total = 1;			/* To start the loop */
  while (i < nstreams && total > 0) {
    total = 0;
    j = i;
    while (j < nstreams && total + streamsize_array[order[j]] < CUTOFF) {
      total += streamsize_array[order[j++]];
    }

    if (total > 0) {
      debug(printf("Merging from %d to %d with %d total elements\n",i,j-1,total));
      assert(ncombined <= mergeinfo->max_nstreams);
      totals[ncombined] = total;
      /* Need an extra value for Sedgesort_uint4 */
      _combined[ncombined] = (UINT4 *) MALLOC_ALIGN((total+1)*sizeof(UINT4));
      out = _combined[ncombined];
      debug_align(printf("Allocating %p to %p -- Malloc of 1 bytes requested from %s:%d\n",
			 out,out,__FILE__,__LINE__));
      for (k = i; k < j; k++) {
	n = streamsize_array[order[k]];
	if (n > 0) {
	  memcpy(out,stream_array[order[k]],n*sizeof(UINT4));
	}

	/* Add diagterm (Could use SIMD here) */
	diagterm = (UINT4) diagterm_array[order[k]];
	for (l = 0; l < n; l++) {
	  out[l] += diagterm;
	}

	out += n;
      }
      Sedgesort_uint4(_combined[ncombined++],total);
      i = j;
    }
  }
  
  ncopied = (nstreams - i) + ncombined;
  heapsize = 2*ncopied - 1;
  
  assert(heapsize < 2*mergeinfo->querylength - 1);
  _heap = mergeinfo->heap;
  nelts = mergeinfo->nelts;
  
  heapi = heapsize;
  /* Handle individual contents: Start with i value from before */
  while (i < nstreams) {
    if ((n = nelts[heapi] = streamsize_array[order[i]]) == 0) {
      _heap[heapi] = (UINT4 *) NULL;

    } else {
      /* Copy to make the merging process non-destructive */
      _heap[heapi] = (UINT4 *) MALLOC_ALIGN(n*sizeof(UINT4));
      out = _heap[heapi];
      debug_align(printf("Allocating %p to %p -- Malloc of 1 bytes requested from %s:%d\n",
			 out,out,__FILE__,__LINE__));
      CHECK_ALIGN(_heap[heapi]);
      memcpy(_heap[heapi],stream_array[order[i]],n*sizeof(UINT4));

      /* Add diagterm (Could use SIMD here) */
      diagterm = (UINT4) diagterm_array[order[i]];
      for (l = 0; l < n; l++) {
	out[l] += diagterm;
      }

#ifdef DEBUG
      printf("Assigning node %d with %d elts:",heapi,nelts[heapi]);
      for (k = 0; k < nelts[heapi]; k++) {
	printf(" %u",_heap[heapi][k]);
      }
      printf("\n");
#endif
    }

    heapi--;
    i++;
  }

  /* Handle combined contents */
  for (i = 0; i < ncombined; i++) {
    _heap[heapi] = _combined[i];
    nelts[heapi] = totals[i];
#ifdef DEBUG
    printf("Assigning node %d with %d elts:",heapi,nelts[heapi]);
    for (k = 0; k < nelts[heapi]; k++) {
      printf(" %u",_heap[heapi][k]);
    }
    printf("\n");
#endif
    heapi--;
  }

  FREE(order);

  return ncopied;
}


#ifdef USE_SIMD_MERGE
/* Input: streams might be aligned or not.  Caller will have to free appropriately */
/* Output: aligned */
UINT4 *
Merge_diagonals (int *nelts1, UINT4 **stream_array, int *streamsize_array,
		 int *diagterm_array, int nstreams, Mergeinfo_uint4_T mergeinfo) {
  /* UINT4 **heap; */
  /* int *nelts; */

  UINT4 *_result, *stream, diagterm;
  int l;
  int ncopied, heapi, heapsize, base, ancestori, pyramid_start, pyramid_end;
  int bits;
#ifdef DEBUG
  int i;
#endif


  if (nstreams == 0) {
    *nelts1 = 0;
    return (UINT4 *) NULL;

  } else if (nstreams == 1) {
    if ((*nelts1 = streamsize_array[0]) == 0) {
      return (UINT4 *) NULL;

    } else {
      stream = stream_array[0];
      diagterm = diagterm_array[0];

      _result = (UINT4 *) MALLOC_ALIGN((*nelts1)*sizeof(UINT4)); /* Output must be aligned */
      debug_align(printf("Allocating %p to %p -- Malloc of 1 bytes requested from %s:%d\n",
			 _result,_result,__FILE__,__LINE__));
      memcpy(_result,stream,(*nelts1)*sizeof(UINT4));

    /* Add diagterm (Could use SIMD here) */
      for (l = 0; l < *nelts1; l++) {
	_result[l] += diagterm;
      }

#ifdef DEBUG
      printf("Merge_diagonals returning result of length %d\n",*nelts1);
      for (i = 0; i < *nelts1; i++) {
	printf("%u\n",_result[i]);
      }
#endif

      return _result;
    }

  } else {
    ncopied = make_diagonals_heap(mergeinfo,stream_array,streamsize_array,diagterm_array,nstreams);
    if (ncopied == 0) {
      *nelts1 = 0;
      return (UINT4 *) NULL;

    } else if (ncopied == 1) {
      *nelts1 = mergeinfo->nelts[1];
      _result = mergeinfo->heap[1];

    } else {
      heapsize = 2*ncopied - 1;	/* also index of last node */
#ifdef HAVE_BUILTIN_CLZ
      bits = 31 - __builtin_clz((unsigned int) heapsize);
#elif defined(HAVE_ASM_BSR)
      asm("bsr %1,%0" : "=r"(bits) : "r"(heapsize));
#else
      bits = 31 - ((heapsize >> 16) ? clz_table[heapsize >> 16] : 16 + clz_table[heapsize]); 
#endif

      base = (1 << bits);
      debug(printf("nstreams %d, ncopied %d, heapsize %d, clz %d, bits %d, base %d\n",
		   nstreams,ncopied,heapsize,__builtin_clz(heapsize),bits,base));
      
      /* Middle pyramids */
      while (base > PYRAMID_SIZE) {
	for (pyramid_start = 2*base - PYRAMID_SIZE, pyramid_end = 2*base - 1; pyramid_start >= base;
	     pyramid_start -= PYRAMID_SIZE, pyramid_end -= PYRAMID_SIZE) {
	  debug(printf("diagonals: pyramid_start %d, pyramid_end %d, ncopied %d\n",pyramid_start,pyramid_end,ncopied));
	  ancestori = pyramid_merge(mergeinfo->heap,ncopied,heapsize,mergeinfo->nelts,pyramid_start,pyramid_end);
	}
	base = ancestori;
      }

      /* Last pyramid */
      pyramid_start = base;
      pyramid_end = 2*base - 1;
      debug(printf("diagonals: pyramid_start %d, pyramid_end %d, ncopied %d\n",pyramid_start,pyramid_end,ncopied));
      /* base = */ pyramid_merge(mergeinfo->heap,ncopied,heapsize,mergeinfo->nelts,pyramid_start,pyramid_end);

      *nelts1 = mergeinfo->nelts[1];
      _result = mergeinfo->heap[1];

      for (heapi = heapsize; heapi > heapsize - ncopied; heapi--) {
	FREE_ALIGN(mergeinfo->heap[heapi]);
      }
    }

#ifdef DEBUG
    printf("Merge_diagonals returning result of length %d\n",*nelts1);
    for (i = 0; i < *nelts1; i++) {
      printf("%u\n",_result[i]);
    }
#endif

    return _result;
  }
}
#endif


/* Fills mergeinfo->heap and mergeinfo->nelts */
static int
make_univdiagonals_heap (Mergeinfo_uint4_T mergeinfo, UINT4 **stream_array,
			 int *streamsize_array, int nstreams) {
  int ncopied, *nelts;
  UINT4 **_heap, **_combined, *out;
  int heapsize, heapi, ncombined;
  int *order, *totals, total, n, i, j, k;

  /* Put smallest streams at the end of the heap, to increase speed */

  /* Note: Sedgesort_order_int overwrites at array[n], so we need to copy */
  memcpy(mergeinfo->streamsize_copy,streamsize_array,nstreams*sizeof(int));
  order = Sedgesort_order_int(mergeinfo->streamsize_copy,nstreams);

  _combined = mergeinfo->combined;
  totals = mergeinfo->totals;
  ncombined = 0;

  /* Combine the smallest streams, up to CUTOFF, and just use Sedgesort */
  i = 0;
  total = 1;			/* To start the loop */
  while (i < nstreams && total > 0) {
    total = 0;
    j = i;
    while (j < nstreams && total + streamsize_array[order[j]] < CUTOFF) {
      total += streamsize_array[order[j++]];
    }

    if (total > 0) {
      debug(printf("Merging from %d to %d with %d total elements\n",i,j-1,total));
      assert(ncombined <= mergeinfo->max_nstreams);
      totals[ncombined] = total;
      /* Need an extra value for Sedgesort_uint4 */
      _combined[ncombined] = (UINT4 *) MALLOC_ALIGN((total+1)*sizeof(UINT4));
      out = _combined[ncombined];
      debug_align(printf("Allocating %p to %p -- Malloc of 1 bytes requested from %s:%d\n",
			 out,out,__FILE__,__LINE__));
      for (k = i; k < j; k++) {
	n = streamsize_array[order[k]];
	if (n > 0) {
	  memcpy(out,stream_array[order[k]],n*sizeof(UINT4));
	}
	out += n;
      }
      Sedgesort_uint4(_combined[ncombined++],total);
#ifdef CHECK_INPUTS
      check_ascending(_combined[ncombined-1],total);
#endif
      i = j;
    }
  }
  
  ncopied = (nstreams - i) + ncombined;
  heapsize = 2*ncopied - 1;
  
  _heap = mergeinfo->heap;
  nelts = mergeinfo->nelts;
  
  heapi = heapsize;
  /* Handle individual contents: Start with i value from before */
  while (i < nstreams) {
    if ((n = nelts[heapi] = streamsize_array[order[i]]) == 0) {
      _heap[heapi] = (UINT4 *) NULL;

    } else {
#ifdef CHECK_INPUTS
      check_ascending(stream_array[order[i]],n);
#endif
    
      /* Copy to make the merging process non-destructive */
      _heap[heapi] = (UINT4 *) MALLOC_ALIGN(n*sizeof(UINT4));
      out = _heap[heapi];
      debug_align(printf("Allocating %p to %p -- Malloc of 1 bytes requested from %s:%d\n",
			 out,out,__FILE__,__LINE__));
      CHECK_ALIGN(_heap[heapi]);
      memcpy(_heap[heapi],stream_array[order[i]],n*sizeof(UINT4));

#ifdef DEBUG
      printf("Assigning node %d with %d elts:",heapi,nelts[heapi]);
      for (k = 0; k < nelts[heapi]; k++) {
	printf(" %u",_heap[heapi][k]);
      }
      printf("\n");
#endif
    }
    heapi--;
    i++;
  }

  /* Handle combined contents */
  for (i = 0; i < ncombined; i++) {
    _heap[heapi] = _combined[i];
    nelts[heapi] = totals[i];
#ifdef DEBUG
    printf("Assigning node %d with %d elts:",heapi,nelts[heapi]);
    for (k = 0; k < nelts[heapi]; k++) {
      printf(" %u",_heap[heapi][k]);
    }
    printf("\n");
#endif
    heapi--;
  }

  FREE(order);

  return ncopied;
}


#ifdef USE_SIMD_MERGE
/* Used for Kmer_anypair_univdiagonals and localdb */
/* Input: elements of stream_array do not have to be aligned */
/* Output: aligned */
UINT4 *
Merge_diagonals_uint4 (int *nelts1, UINT4 **stream_array, int *streamsize_array,
		       int nstreams, Mergeinfo_uint4_T mergeinfo) {
  /* UINT4 **heap; */
  /* int *nelts; */
  UINT4 *_result, *stream;
  int ncopied, heapi, heapsize, base, ancestori, pyramid_start, pyramid_end;
  int bits;
#ifdef DEBUG
  int i;
#endif


  if (nstreams == 0) {
    *nelts1 = 0;
    return (UINT4 *) NULL;

  } else if (nstreams == 1) {
    if ((*nelts1 = streamsize_array[0]) == 0) {
      return (UINT4 *) NULL;
      
    } else {
      stream = stream_array[0];
#ifdef CHECK_INPUTS
      check_ascending(stream,*nelts1);
#endif

      _result = (UINT4 *) MALLOC_ALIGN((*nelts1)*sizeof(UINT4)); /* Output must be aligned */
      debug_align(printf("Allocating %p to %p -- Malloc of 1 bytes requested from %s:%d\n",
			 _result,_result,__FILE__,__LINE__));
      memcpy(_result,stream,(*nelts1)*sizeof(UINT4));
      return _result;
    }

  } else {
    ncopied = make_univdiagonals_heap(mergeinfo,stream_array,streamsize_array,nstreams);
    if (ncopied == 0) {
      *nelts1 = 0;
      return (UINT4 *) NULL;
      
    } else if (ncopied == 1) {
      *nelts1 = mergeinfo->nelts[1];
      _result = mergeinfo->heap[1];

    } else {
      heapsize = 2*ncopied - 1;	/* also index of last node */
#ifdef HAVE_BUILTIN_CLZ
      bits = 31 - __builtin_clz((unsigned int) heapsize);
#elif defined(HAVE_ASM_BSR)
      asm("bsr %1,%0" : "=r"(bits) : "r"(heapsize));
#else
      bits = 31 - ((heapsize >> 16) ? clz_table[heapsize >> 16] : 16 + clz_table[heapsize]); 
#endif

      base = (1 << bits);
      debug(printf("nstreams %d, ncopied %d, heapsize %d, clz %d, bits %d, base %d\n",
		   nstreams,ncopied,heapsize,__builtin_clz(heapsize),bits,base));
      
      /* Middle pyramids */
      while (base > PYRAMID_SIZE) {
	for (pyramid_start = 2*base - PYRAMID_SIZE, pyramid_end = 2*base - 1; pyramid_start >= base;
	     pyramid_start -= PYRAMID_SIZE, pyramid_end -= PYRAMID_SIZE) {
	  debug(printf("diagonals: pyramid_start %d, pyramid_end %d, ncopied %d\n",pyramid_start,pyramid_end,ncopied));
	  ancestori = pyramid_merge(mergeinfo->heap,ncopied,heapsize,mergeinfo->nelts,pyramid_start,pyramid_end);
	}
	base = ancestori;
      }

      /* Last pyramid */
      pyramid_start = base;
      pyramid_end = 2*base - 1;
      debug(printf("diagonals: pyramid_start %d, pyramid_end %d, ncopied %d\n",pyramid_start,pyramid_end,ncopied));
      /* base = */ pyramid_merge(mergeinfo->heap,ncopied,heapsize,mergeinfo->nelts,pyramid_start,pyramid_end);

      *nelts1 = mergeinfo->nelts[1];
      _result = mergeinfo->heap[1];

      for (heapi = heapsize; heapi > heapsize - ncopied; heapi--) {
	FREE_ALIGN(mergeinfo->heap[heapi]);
      }
    }
  }


#ifdef DEBUG
  printf("Merge_diagonals returning result of length %d\n",*nelts1);
  for (i = 0; i < *nelts1; i++) {
    printf("%u\n",_result[i]);
  }
#endif

  return _result;
}
#endif

