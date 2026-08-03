#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "merge-diagonals-simd-uint8.h"
#include "assert.h"
#include "mem.h"
#include "popcount.h"		/* For clz_table */
#include "sedgesort.h"
#include "merge-uint8.h"

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

#define GETPOS(high,low) (((Univcoord_T) high << 32) + low)


#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


#define PARENT(i) (i >> 1)
#define LEFT(i) (i << 1)
#define RIGHT(i) ((i << 1) | 1)

#if defined(HAVE_AVX512) || defined(HAVE_AVX2)
/* Each help elt is aligned */
static int
pyramid_merge (UINT8 **_heap, int nstreams, int heapsize, int *nelts,
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
      _heap[PARENT(nodei)] = Merge_uint8(/*dest*/NULL,_heap[nodei-1],_heap[nodei],nelts[nodei-1],nelts[nodei]);
      CHECK_ALIGN(_heap[PARENT(nodei)]);
      nelts[PARENT(nodei)] = nelts[nodei-1] + nelts[nodei];
      debug2(printf("Created list %p of length %d at node %d\n",
		    _heap[PARENT(nodei)],nelts[PARENT(nodei)],PARENT(nodei)));

#ifdef DEBUG
      for (i = 0; i < nelts[PARENT(nodei)]; i++) {
	printf("%llu\n",_heap[PARENT(nodei)][i]);
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
#endif


#if defined(HAVE_AVX512) || defined(HAVE_AVX2)
/* Fills mergeinfo->heap and mergeinfo->nelts */
static int
make_diagonals_heap (Mergeinfo_uint8_T mergeinfo,
		     unsigned char **stream_high_array, UINT4 **stream_low_array,
		     int *streamsize_array, int *diagterm_array, int nstreams) {
  int ncopied, *nelts;
  UINT8 **_heap, **_combined, *out, diagterm;
  int heapsize, heapi, ncombined;
  int *order, *totals, total, n, i, j, k, l;

  debug(printf("Entered make_diagonals_heap\n"));

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
      totals[ncombined] = total;
      /* Need an extra value for Sedgesort_uint8 */
      _combined[ncombined] = (UINT8 *) MALLOC_ALIGN((total+1)*sizeof(UINT8));
      out = _combined[ncombined];
      for (k = i; k < j; k++) {
	n = streamsize_array[order[k]];

	/* Combine stream_high and stream_low and add diagterm (Could use SIMD here) */
	diagterm = (UINT8) diagterm_array[order[k]];
	for (l = 0; l < n; l++) {
	  out[l] = GETPOS(stream_high_array[order[k]][l],stream_low_array[order[k]][l]) + diagterm;
	}

	out += n;
      }
      Sedgesort_uint8(_combined[ncombined++],total);
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
      _heap[heapi] = (UINT8 *) NULL;
      
    } else {
      /* Copy to make the merging process non-destructive */
      _heap[heapi] = (UINT8 *) MALLOC_ALIGN(n*sizeof(UINT8));
      out = _heap[heapi];
      CHECK_ALIGN(_heap[heapi]);
    }

    /* Combine stream_high and stream_low and add diagterm (Could use SIMD here) */
    diagterm = (UINT8) diagterm_array[order[i]];
    for (l = 0; l < n; l++) {
      out[l] = GETPOS(stream_high_array[order[i]][l],stream_low_array[order[i]][l]) + diagterm;
    }

#ifdef DEBUG
    printf("Assigning node %d with %d elts:",heapi,nelts[heapi]);
    for (k = 0; k < nelts[heapi]; k++) {
      printf(" %llu",_heap[heapi][k]);
    }
    printf("\n");
#endif
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
      printf(" %llu",_heap[heapi][k]);
    }
    printf("\n");
#endif
    heapi--;
  }

  FREE(order);

  return ncopied;
}
#endif


/* For heap code, see merge-diagonals-heap-simd8.c */
#ifdef USE_SIMD_MERGE
/* Input: streams might be aligned or not.  Caller will have to free appropriately */
/* Output: aligned */
Univcoord_T *
Merge_diagonals_large (int *nelts1, unsigned char **stream_high_array, UINT4 **stream_low_array,
		       int *streamsize_array, int *diagterm_array, int nstreams,
		       Mergeinfo_uint8_T mergeinfo) {
  /* Univcoord_T **heap; */
  /* int *nelts; */

  Univcoord_T *result, diagterm;
  unsigned char *stream_high;
  UINT4 *stream_low;
  int l;
  int ncopied, heapi, heapsize, base, ancestori, pyramid_start, pyramid_end;
  int bits;
#ifdef DEBUG
  int i;
#endif


  if (nstreams == 0) {
    *nelts1 = 0;
    return (Univcoord_T *) NULL;

  } else if (nstreams == 1) {
    if ((*nelts1 = streamsize_array[0]) == 0) {
      return (Univcoord_T *) NULL;
      
    } else {
      stream_high = stream_high_array[0];
      stream_low = stream_low_array[0];
      diagterm = diagterm_array[0];

      result = (UINT8 *) MALLOC_ALIGN((*nelts1)*sizeof(UINT8)); /* Output must be aligned */
      /* Combine stream_high and stream_low and add diagterm (Could use SIMD here) */
      for (l = 0; l < *nelts1; l++) {
	result[l] = GETPOS(stream_high[l],stream_low[l]) + diagterm;
      }

      return result;
    }

  } else {
    ncopied = make_diagonals_heap(mergeinfo,stream_high_array,stream_low_array,
				  streamsize_array,diagterm_array,nstreams);
    if (ncopied == 0) {
      *nelts1 = 0;
      return (Univcoord_T *) NULL;

    } else if (ncopied == 1) {
      *nelts1 = mergeinfo->nelts[1];
      result = mergeinfo->heap[1];

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
      result = mergeinfo->heap[1];

      for (heapi = heapsize; heapi > heapsize - ncopied; heapi--) {
	FREE_ALIGN(mergeinfo->heap[heapi]);
      }
    }
  }


#ifdef DEBUG
  printf("Merge_diagonals_large returning result of length %d\n",*nelts1);
  for (i = 0; i < *nelts1; i++) {
    printf("%llu\n",result[i]);
  }
#endif

  return result;
}
#endif


/* Fills mergeinfo->heap and mergeinfo->nelts */
static int
make_univdiagonals_heap (Mergeinfo_uint8_T mergeinfo, Univcoord_T **stream_array,
			 int *streamsize_array, int nstreams) {
  int ncopied, *nelts;
  Univcoord_T **_heap, **_combined, *out;
  int heapsize, heapi, ncombined;
  int *order, *totals, total, n, i, j, k;

  debug(printf("Entered make_univdiagonals_heap\n"));
#ifdef DEBUG
  for (i = 0; i < nstreams; i++) {
    printf("Stream %d:",i);
    for (j = 0; j < streamsize_array[i]; j++) {
      printf(" %llu",stream_array[i][j]);
    }
    printf("\n");
  }
#endif


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
      totals[ncombined] = total;
      /* Need an extra value for Sedgesort_uint8 */
      _combined[ncombined] = (Univcoord_T *) MALLOC_ALIGN((total+1)*sizeof(Univcoord_T));
      out = _combined[ncombined];
      for (k = i; k < j; k++) {
	n = streamsize_array[order[k]];
	if (n > 0) {
	  memcpy(out,stream_array[order[k]],n*sizeof(Univcoord_T));
	}
	out += n;
      }
      Sedgesort_uint8(_combined[ncombined++],total);
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
      _heap[heapi] = (Univcoord_T *) NULL;

    } else {
      /* Copy to make the merging process non-destructive */
      _heap[heapi] = (Univcoord_T *) MALLOC_ALIGN(n*sizeof(Univcoord_T));
      out = _heap[heapi];
      CHECK_ALIGN(_heap[heapi]);
      memcpy(_heap[heapi],stream_array[order[i]],n*sizeof(Univcoord_T));
    }
 
#ifdef DEBUG
    printf("Assigning node %d with %d elts:",heapi,nelts[heapi]);
    for (k = 0; k < nelts[heapi]; k++) {
      printf(" %llu",_heap[heapi][k]);
    }
    printf("\n");
#endif
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
      printf(" %llu",_heap[heapi][k]);
    }
    printf("\n");
#endif
    heapi--;
  }

  FREE(order);

  return ncopied;
}


/* For heap version, see merge-diagonals-heap-uint8.c */

#ifdef USE_SIMD_MERGE
/* Used for Kmer_anypair_univdiagonals */
/* kmer-search.c calls Merge_diagonals_large with stream_high_list and
   stream_low_list.  localdb.c calls Merge_diagonals_uint4. */
Univcoord_T *
Merge_diagonals_uint8 (int *nelts1, Univcoord_T **stream_array, int *streamsize_array,
		       int nstreams, Mergeinfo_uint8_T mergeinfo) {
  /* Univcoord_T **heap; */
  /* int *nelts; */

  Univcoord_T *_result, *stream;
  int ncopied, heapi, heapsize, base, ancestori, pyramid_start, pyramid_end;
  int bits;
#ifdef DEBUG
  int i;
#endif


  if (nstreams == 0) {
    *nelts1 = 0;
    return (Univcoord_T *) NULL;

  } else if (nstreams == 1) {
    if ((*nelts1 = streamsize_array[0]) == 0) {
      return (Univcoord_T *) NULL;

    } else {
      stream = stream_array[0];

      _result = (Univcoord_T *) MALLOC_ALIGN((*nelts1)*sizeof(Univcoord_T)); /* Output must be aligned */
      memcpy(_result,stream,(*nelts1)*sizeof(Univcoord_T));

#ifdef DEBUG
      printf("Merge_diagonals_uint8 returning result of length %d\n",*nelts1);
      for (i = 0; i < *nelts1; i++) {
	printf("%llu\n",_result[i]);
      }
#endif

      return _result;
    }

  } else {
    ncopied = make_univdiagonals_heap(mergeinfo,stream_array,streamsize_array,nstreams);
    if (ncopied == 0) {
      *nelts1 = 0;
      return (Univcoord_T *) NULL;

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
    printf("Merge_diagonals_uint8 returning result of length %d\n",*nelts1);
    for (i = 0; i < *nelts1; i++) {
      printf("%llu\n",_result[i]);
    }
#endif

    return _result;
  }
}
#endif


