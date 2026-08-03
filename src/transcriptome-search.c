#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "transcriptome-search.h"

#include "assert.h"
#include "trpath.h"
#include "trpath-solve.h"

#include "path.h"
#include "path-eval.h"

#include "mem.h"
#include "types.h"
#include "chrnum.h"
#include "reader.h"
#include "oligo.h"

#include "sedgesort.h"
#include "intersect-small.h"
#include "intersect-approx-uint4.h" /* For Intersect_approx_uint4 */
#include "intersect-indices-small.h" /* For Transcriptome_anypair */

#include "merge-method.h"
#ifdef USE_HEAP_MERGE
#include "merge-diagonals-heap-uint4.h"
#elif defined(USE_SIMD_MERGE)
#include "merge-diagonals-simd-uint4.h"
#endif

#include "transcript.h"
#include "junction.h"
#include "genomebits_count.h"
#include "genomebits_kmer.h"
#include "genomebits_trim.h"
#include "spliceends.h"


/* Merging faster than count table */
#define USE_MERGE 1
#define SUBOPT 3


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Processing merged diagonals */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* Exon and path overlap */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Transcriptome_anypair */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Transcriptome alignment details*/
#ifdef DEBUG2A
#define debug2a(x) x
#else
#define debug2a(x)
#endif

/* Transcriptome_exact */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* Transcriptome_search_approx */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif

/* Transcriptome_search_end */
#ifdef DEBUG5
#define debug5(x) x
#else
#define debug5(x)
#endif

/* Transcriptome_prevalent */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif


static int index1interval_tr;
static int index1part_tr;

static Indexdb_T tr_indexdb;

static Transcriptome_T transcriptome;
static EF64_T transcript_ef64;

static Genomebits_T genomebits;
static Genomebits_T genomebits_alt;
static Genomebits_T transcriptomebits;

static int leftreadshift_tr;
static Oligospace_T oligobase_mask_tr;


#if 0
/* Overwrite values array */
/* Repetitive oligos can calse false indel diagonals */
static int
most_prevalent_trcoord (int *nloci, int **counts, Trcoord_T *values, int nvalues) {
  int max_count, count, *count_ptr;
  Trcoord_T *out, *ptr, *end, *first;

  assert(nvalues > 0);
  count_ptr = *counts = (int *) MALLOC(nvalues*sizeof(int));

  ptr = out = &(values[0]);	/* Reset */
  end = &(values[nvalues]);

  max_count = 1;
  while (ptr < end) {
    first = ptr;
    debug0(printf("DIAGONAL %u",*first));
    if (ptr + max_count - 1 >= end) {
      /* End of list fails */
      debug0(printf(" => Goes beyond end of list\n"));
      ptr = end;

    } else if (ptr[max_count - 1] != *first) {
      /* Fails */
      debug0(printf(" => Fails\n"));
      if (ptr[max_count - 1] != ptr[max_count - 2]) {
	/* Can jump immediately */
	ptr += max_count - 1;
      } else {
	/* Advance forward until we see a new value */
	while (ptr < end && *ptr == *first) {
	  ptr++;
	}
      }

    } else {
      /* Contender */
      ptr += max_count;
      while (ptr < end && *ptr == *first) {
	ptr++;
      }
      debug0(printf(" => Count %d, index %d => printing\n",ptr - first,first - values));
      if ((count = ptr - first) > max_count + SUBOPT) {
	out = &(values[0]);	/* Reset */
	count_ptr = *counts;	/* Reset */
	max_count = count;
      } else if (count > max_count) {
	max_count = count;
      }
      *out++ = *first;
      *count_ptr++ = count;
    }
  }

  *nloci = out - &(values[0]);
  return max_count;
}
#endif


#define T Trpath_T


#if 0
/* Expecting transcriptome to not be a large genome */
/* Use oligos computed by Stage1_fill_all_oligos_tr */
/* Assigns own positions, so does not need a procedure for Stage1_fill_all_positions_tr */
void
Transcriptome_search_complete (int *found_score,

			       List_T *sense_trpaths, List_T *antisense_trpaths,

			       Stage1_T stage1, int querylength, Indelinfo_T indelinfo,

			       Mergeinfo_uint4_T mergeinfo,
			       Trcoord_T **tplus_stream_array, int *tplus_streamsize_array, int *tplus_diagterm_array,
			       Trcoord_T **tminus_stream_array, int *tminus_streamsize_array, int *tminus_diagterm_array,
			       int *mismatch_positions_alloc, Compress_T query_compress_fwd, Compress_T query_compress_rev,
			       int max_insertionlen, int max_deletionlen,
			       Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			       Listpool_T listpool, Trpathpool_T trpathpool, Pathpool_T pathpool,
			       Hitlistpool_T hitlistpool, Method_T method) {
  Trpath_T trpath;

  int max_count, max_tplus_count, max_tminus_count, *tplus_counts, *tminus_counts;

  Trcoord_T trdiagonal;
  Trcoord_T *_tplus_diagonals, *_tminus_diagonals; /* aligned */
  int n_tplus_loci, n_tminus_loci, n_tplus_diagonals, n_tminus_diagonals,
    total_ndiagonals_tplus = 0, total_ndiagonals_tminus = 0, i;
  
  Trcoord_T *positions;
  int npositions;
  int tplus_streami = 0, tminus_streami = 0;

  int querypos, query_lastpos;

  Trnum_T trnum;
  Trcoord_T troffset, trhigh;


  query_lastpos = querylength - index1part_tr;
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    if (stage1->tr_validp[querypos] == false) {
      /* Skip */
    } else {
      if (stage1->tr_plus_retrievedp[querypos] == true) {
	positions = stage1->tr_plus_positions[querypos];
	npositions = stage1->tr_plus_npositions[querypos];
      } else {
	npositions = stage1->tr_plus_npositions[querypos] =
	  Indexdb_ptr(&stage1->tr_plus_positions[querypos],tr_indexdb,stage1->tr_forward_oligos[querypos]);
	positions = stage1->tr_plus_positions[querypos];
	stage1->tr_plus_retrievedp[querypos] = true;
      }
      if (npositions > 0) {
	tplus_stream_array[tplus_streami] = positions;
	tplus_streamsize_array[tplus_streami] = npositions;
	tplus_diagterm_array[tplus_streami] = querylength - querypos;
	total_ndiagonals_tplus += npositions;
	tplus_streami++;
      }

      if (stage1->tr_minus_retrievedp[querypos] == true) {
	positions = stage1->tr_minus_positions[querypos];
	npositions = stage1->tr_minus_npositions[querypos];
      } else {
	npositions = stage1->tr_minus_npositions[querypos] =
	  Indexdb_ptr(&stage1->tr_minus_positions[querypos],tr_indexdb,stage1->tr_revcomp_oligos[querypos]);
	positions = stage1->tr_minus_positions[querypos];
	stage1->tr_minus_retrievedp[querypos] = true;
      }
      if (npositions > 0) {
	tminus_stream_array[tminus_streami] = positions;
	tminus_streamsize_array[tminus_streami] = npositions;
	tminus_diagterm_array[tminus_streami] = querypos + index1part_tr;
	total_ndiagonals_tminus += npositions;
	tminus_streami++;
      }
    }
  }


  _tplus_diagonals = Merge_diagonals(&n_tplus_diagonals,tplus_stream_array,
				     tplus_streamsize_array,tplus_diagterm_array,/*nstreams*/tplus_streami,
				     mergeinfo);

  if (n_tplus_diagonals == 0) {
    max_tplus_count = 0;
    n_tplus_loci = 0;
    tplus_counts = (int *) NULL;
  } else {
    max_tplus_count = most_prevalent_trcoord(&n_tplus_loci,&tplus_counts,_tplus_diagonals,n_tplus_diagonals);
  }

#ifdef DEBUG2
  printf("max_tplus_count: %d\n",max_tplus_count);
  for (i = 0; i < n_tplus_loci; i++) {
    printf("PLUS_LOCUS %u, count %d\n",_tplus_diagonals[i],tplus_counts[i]);
  }
#endif

  _tminus_diagonals = Merge_diagonals(&n_tminus_diagonals,tminus_stream_array,
				     tminus_streamsize_array,tminus_diagterm_array,/*nstreams*/tminus_streami,
				     mergeinfo);

  if (n_tminus_diagonals == 0) {
    max_tminus_count = 0;
    n_tminus_loci = 0;
    tminus_counts = (int *) NULL;
  } else {
    max_tminus_count = most_prevalent_trcoord(&n_tminus_loci,&tminus_counts,_tminus_diagonals,n_tminus_diagonals);
  }

#ifdef DEBUG2
  printf("max_tminus_count: %d\n",max_tminus_count);
  for (i = 0; i < n_tminus_loci; i++) {
    printf("MINUS_LOCUS %u, count %d\n",_tminus_diagonals[i],tminus_counts[i]);
  }
#endif


  if (max_tplus_count > max_tminus_count) {
    max_count = max_tplus_count;
  } else {
    max_count = max_tminus_count;
  }

  /* tplusp == true */
  debug2(printf("tplus loci %d\n",n_tplus_loci));
  if (max_tplus_count < max_count - SUBOPT) {
    debug2(printf("tplus not good enough\n"));
  } else {
    for (i = 0; i < n_tplus_loci; i++) {
      trdiagonal = _tplus_diagonals[i];
      if (trdiagonal /*+qstart:0*/ < (Trcoord_T) querylength) {
	/* Skip */
      } else if (tplus_counts[i] < max_count - SUBOPT) {
	debug2(printf("Skipping trdiagonal %u with %d counts\n",trdiagonal,tplus_counts[i]));
      } else {
	trnum = EF64_trnum(&troffset,&trhigh,transcript_ef64,trdiagonal - querylength,trdiagonal);
	chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
	if ((trpath = Trpath_solve_from_diagonals(&(*found_score),trdiagonal,
						  /*middle_trdiagonal_qstart*/0,/*middle_trdiagonal_qend*/querylength,
						  /*middle_nmismatches*/0,/*qstart_trdiag*/NULL,/*qend_trdiag*/NULL,
						  /*tplusp*/true,querylength,/*query_compress_tr*/query_compress_fwd,
						  mismatch_positions_alloc,trnum,troffset,trhigh,
						  chrnum,transcript_genestrand,indelinfo,
						  intlistpool,uintlistpool,listpool,
						  trpathpool,pathpool,method)) != NULL) {
	  *sense_trpaths = Hitlist_push(*sense_trpaths,hitlistpool,(void *) trpath
					hitlistpool_trace(__FILE__,__LINE__));
	}
      }
    }
  }
    
  FREE(tplus_counts);
  FREE_ALIGN(_tplus_diagonals);


  /* tplusp == false */
  debug2(printf("tminus loci %d\n",n_tminus_loci));
  if (max_tminus_count < max_count - SUBOPT) {
    debug2(printf("tminus not good enough\n"));
  } else {
    for (i = 0; i < n_tminus_loci; i++) {
      trdiagonal = _tminus_diagonals[i];
      if (trdiagonal /*+qstart:0*/ < (Trcoord_T) querylength) {
	/* Skip */
      } else if (tminus_counts[i] < max_count - SUBOPT) {
	debug2(printf("Skipping trdiagonal %u with %d counts\n",trdiagonal,tminus_counts[i]));
      } else {
	trnum = EF64_trnum(&troffset,&trhigh,transcript_ef64,trdiagonal - querylength,trdiagonal);
	chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
	if ((trpath = Trpath_solve_from_diagonals(&(*found_score),trdiagonal,
						  /*middle_trdiagonal_qstart*/0,/*middle_trdiagonal_qend*/querylength,
						  /*middle_nmismatches*/0,/*qstart_trdiag*/NULL,/*qend_trdiag*/NULL,
						  /*tplusp*/false,querylength,/*query_compress_tr*/query_compress_rev,
						  mismatch_positions_alloc,trnum,troffset,trhigh,
						  chrnum,transcript_genestrand,indelinfo,
						  intlistpool,uintlistpool,listpool,
						  trpathpool,pathpool,method)) != NULL) {
	  *antisense_trpaths = Hitlist_push(*antisense_trpaths,hitlistpool,(void *) trpath
					    hitlistpool_trace(__FILE__,__LINE__));
	}
      }
    }
  }

  FREE(tminus_counts);
  FREE_ALIGN(_tminus_diagonals);

  return;
}
#endif


int
Transcriptome_exact1 (Trnum_T **sense_trnums, Trcoord_T **sense_troffsets, Trcoord_T **sense_trhighs,
		      Trcoord_T **_sense_trdiagonals, int *n_sense_trdiagonals,
		      Trnum_T **antisense_trnums, Trcoord_T **antisense_troffsets, Trcoord_T **antisense_trhighs,
		      Trcoord_T **_antisense_trdiagonals, int *n_antisense_trdiagonals,
		      Stage1_T stage1, int querylength) {

  Trcoord_T *tplus_positions_end5 = NULL, *tminus_positions_end5 = NULL,
    *tplus_positions_end3 = NULL, *tminus_positions_end3 = NULL;
  int n_tplus_positions_end5 = 0, n_tminus_positions_end5 = 0,
    n_tplus_positions_end3 = 0, n_tminus_positions_end3 = 0;
  int tplus_diagterm_end5, tminus_diagterm_end5, tplus_diagterm_end3, tminus_diagterm_end3;

  Trcoord_T trdiagonal;
  int query_lastpos;
  int k, i;


  debug3(printf("\n\n***Transcriptome_exact\n"));
  Stage1_fill_all_oligos_tr(stage1,querylength);
  Stage1_fill_all_positions_tr(stage1,querylength);

  query_lastpos = querylength - index1part_tr;
  if (stage1->tr_validp[0] == false || stage1->tr_validp[query_lastpos] == false) {

    *sense_trnums = (Trnum_T *) NULL;
    *_sense_trdiagonals = (Trcoord_T *) NULL;
    *n_sense_trdiagonals = 0;

    *antisense_trnums = (Trnum_T *) NULL;
    *_antisense_trdiagonals = (Trcoord_T *) NULL;
    *n_antisense_trdiagonals = 0;

  } else {
    tplus_positions_end5 = stage1->tr_plus_positions[0];
    n_tplus_positions_end5 = stage1->tr_plus_npositions[0];
    tplus_diagterm_end5 = querylength /* - 0, querypos */;

    tminus_positions_end5 = stage1->tr_minus_positions[0];
    n_tminus_positions_end5 = stage1->tr_minus_npositions[0];
    tminus_diagterm_end5 = /* 0, querypos + */ index1part_tr;

    tplus_positions_end3 = stage1->tr_plus_positions[query_lastpos];
    n_tplus_positions_end3 = stage1->tr_plus_npositions[query_lastpos];
    tplus_diagterm_end3 = index1part_tr /* querylength - querypos */;

    tminus_positions_end3 = stage1->tr_minus_positions[query_lastpos];
    n_tminus_positions_end3 = stage1->tr_minus_npositions[query_lastpos];
    tminus_diagterm_end3 = querylength /* querypos + index1part_tr */;

    *_sense_trdiagonals = Intersect_small(&(*n_sense_trdiagonals),
					  tplus_positions_end5,n_tplus_positions_end5,tplus_diagterm_end5,
					  tplus_positions_end3,n_tplus_positions_end3,tplus_diagterm_end3,
					  /*alignp*/true);
    *_antisense_trdiagonals = Intersect_small(&(*n_antisense_trdiagonals),
					      tminus_positions_end5,n_tminus_positions_end5,tminus_diagterm_end5,
					      tminus_positions_end3,n_tminus_positions_end3,tminus_diagterm_end3,
					      /*alignp*/true);

    if ((*n_sense_trdiagonals) == 0) {
      *sense_trnums = (Trnum_T *) NULL;
      *sense_troffsets = (Trcoord_T *) NULL;
      *sense_trhighs = (Trcoord_T *) NULL;

    } else {
      *sense_trnums = (Trnum_T *) MALLOC((*n_sense_trdiagonals) * sizeof(Trnum_T));
      *sense_troffsets = (Trcoord_T *) MALLOC((*n_sense_trdiagonals) * sizeof(Trcoord_T));
      *sense_trhighs = (Trcoord_T *) MALLOC((*n_sense_trdiagonals) * sizeof(Trcoord_T));

      k = 0;
      for (i = 0; i < (*n_sense_trdiagonals); i++) {
	trdiagonal = (*_sense_trdiagonals)[i];
	if (trdiagonal /*+qstart:0*/ < (Trcoord_T) querylength) {
	  /* Skip */
	} else {
	  (*_sense_trdiagonals)[k] = trdiagonal;
	  (*sense_trnums)[k] = EF64_trnum(&(*sense_troffsets)[k],&(*sense_trhighs)[k],
					  transcript_ef64,trdiagonal - querylength,trdiagonal);
	  debug3(printf("sense trnum: %u  trdiagonal: %u\n",(*sense_trnums)[k],trdiagonal));
	  k++;
	}
      }

      if ((*n_sense_trdiagonals = k) == 0) {
	FREE(*sense_trhighs);
	FREE(*sense_troffsets);
	FREE(*sense_trnums);
	*sense_trnums = (Trnum_T *) NULL;
	*sense_troffsets = (Trcoord_T *) NULL;
	*sense_trhighs = (Trcoord_T *) NULL;
      }
    }

    if ((*n_antisense_trdiagonals) == 0) {
      *antisense_trnums = (Trnum_T *) NULL;
      *antisense_troffsets = (Trcoord_T *) NULL;
      *antisense_trhighs = (Trcoord_T *) NULL;

    } else {
      *antisense_trnums = (Trnum_T *) MALLOC((*n_antisense_trdiagonals) * sizeof(Trnum_T));
      *antisense_troffsets = (Trcoord_T *) MALLOC((*n_antisense_trdiagonals) * sizeof(Trcoord_T));
      *antisense_trhighs = (Trcoord_T *) MALLOC((*n_antisense_trdiagonals) * sizeof(Trcoord_T));

      k = 0;
      for (i = 0; i < (*n_antisense_trdiagonals); i++) {
	trdiagonal = (*_antisense_trdiagonals)[i];
	debug3(printf("tminus trdiagonal: %u\n",trdiagonal));
	if (trdiagonal /*+qstart:0*/ < (Trcoord_T) querylength) {
	  /* Skip */
	} else {
	  (*_antisense_trdiagonals)[k] = trdiagonal;
	  (*antisense_trnums)[k] = EF64_trnum(&(*antisense_troffsets)[k],&(*antisense_trhighs)[k],
					      transcript_ef64,trdiagonal - querylength,trdiagonal);
	  debug3(printf("antisense trnum: %u  trdiagonal: %u\n",(*antisense_trnums)[k],trdiagonal));
	  k++;
	}
      }

      if ((*n_antisense_trdiagonals = k) == 0) {
	FREE(*antisense_trhighs);
	FREE(*antisense_troffsets);
	FREE(*antisense_trnums);
	*antisense_trnums = (Trnum_T *) NULL;
	*antisense_troffsets = (Trcoord_T *) NULL;
	*antisense_trhighs = (Trcoord_T *) NULL;
      }
    }
  }

  debug3(printf("***Transcriptome_exact1: %d sense and %d antisense trdiagonals\n",
		*n_sense_trdiagonals,*n_antisense_trdiagonals));
	 
  return (*n_sense_trdiagonals) + (*n_antisense_trdiagonals);
}


#if 0
static int
make_unique (Trcoord_T *diagonals, int ndiagonals) {
  int k = 0, i, j;

  i = 0;
  while (i < ndiagonals) {
    j = i + 1;
    while (j < ndiagonals && diagonals[j] == diagonals[i]) {
      j++;
    }
    diagonals[k++] = diagonals[i];
    
    i = j;
  }

  return k;
}
#endif


#if 0
static int
keep_duplicates (Trcoord_T *diagonals, int ndiagonals) {
  int k = 0, i, j;

  i = 0;
  while (i < ndiagonals) {
    j = i + 1;
    while (j < ndiagonals && diagonals[j] == diagonals[i]) {
      j++;
    }
    if (j - i > 1) {
      diagonals[k++] = diagonals[i];
    }
    
    i = j;
  }

  return k;
}
#endif


#if 0
/* Can give rise to ends that don't match, so we need to trim the ends */
/* _sense_trdiagonals and _antisense_trdiagonals are aligned */
int
Transcriptome_exact2 (Trnum_T **sense_trnums, Trcoord_T **sense_troffsets, Trcoord_T **sense_trhighs,
		      Trcoord_T **_sense_trdiagonals, int *n_sense_trdiagonals,
		      Trnum_T **antisense_trnums, Trcoord_T **antisense_troffsets, Trcoord_T **antisense_trhighs,
		      Trcoord_T **_antisense_trdiagonals, int *n_antisense_trdiagonals,
		      Stage1_T stage1, int querylength) {

  int ndiagonals;
  int nstreams;
  Trcoord_T trdiagonal;
  int query_lastpos;
  int k, i;


  debug3(printf("\n\n***Transcriptome_exact2\n"));
  /* assertions fail, but we have called Stage1_trdiagonals_gc */
  /* assert(*_sense_trdiagonals == NULL); */
  /* assert(*antisense_trdiagonals == NULL); */

  Stage1_dump_tr(stage1,querylength);


  query_lastpos = querylength - index1part_tr;

  /* plus */
  nstreams = 0;
  stage1->tplus_stream_array[nstreams] = stage1->tr_plus_positions[0];
  stage1->tplus_streamsize_array[nstreams] = stage1->tr_plus_npositions[0];
  stage1->tplus_diagterm_array[nstreams++] = /*diagterm*/querylength - 0;
  if (index1part_tr <= query_lastpos) {
    stage1->tplus_stream_array[nstreams] = stage1->tr_plus_positions[index1part_tr];
    stage1->tplus_streamsize_array[nstreams] = stage1->tr_plus_npositions[index1part_tr];
    stage1->tplus_diagterm_array[nstreams++] = /*diagterm*/querylength - index1part_tr;
  }
  
  stage1->tplus_stream_array[nstreams] = stage1->tr_plus_positions[query_lastpos];
  stage1->tplus_streamsize_array[nstreams] = stage1->tr_plus_npositions[query_lastpos];
  stage1->tplus_diagterm_array[nstreams++] = /*diagterm*/querylength - query_lastpos;
  if (query_lastpos - index1part_tr >= 0) {
    stage1->tplus_stream_array[nstreams] = stage1->tr_plus_positions[query_lastpos - index1part_tr];
    stage1->tplus_streamsize_array[nstreams] = stage1->tr_plus_npositions[query_lastpos - index1part_tr];
    stage1->tplus_diagterm_array[nstreams++] = /*diagterm*/querylength - (query_lastpos - index1part_tr);
  }

  *_sense_trdiagonals = Merge_diagonals(&ndiagonals,
				       stage1->tplus_stream_array,
				       stage1->tplus_streamsize_array,
				       stage1->tplus_diagterm_array,
				       nstreams,stage1->mergeinfo_tr); 

#if 0
  /* ? This looks for 30-mers at the ends */
  *n_sense_trdiagonals = keep_duplicates(*_sense_trdiagonals,ndiagonals);
#else
  *n_sense_trdiagonals = make_unique(*_sense_trdiagonals,ndiagonals);
#endif


  if ((*n_sense_trdiagonals) == 0) {
    *sense_trnums = (Trnum_T *) NULL;
    *sense_troffsets = (Trcoord_T *) NULL;
    *sense_trhighs = (Trcoord_T *) NULL;

  } else {
    *sense_trnums = (Trnum_T *) MALLOC((*n_sense_trdiagonals) * sizeof(Trnum_T));
    *sense_troffsets = (Trcoord_T *) MALLOC((*n_sense_trdiagonals) * sizeof(Trcoord_T));
    *sense_trhighs = (Trcoord_T *) MALLOC((*n_sense_trdiagonals) * sizeof(Trcoord_T));

    k = 0;
    for (i = 0; i < (*n_sense_trdiagonals); i++) {
      trdiagonal = (*_sense_trdiagonals)[i];
      if (trdiagonal /*+qstart:0*/ < (Trcoord_T) querylength) {
	/* Skip */
	debug3(printf("Skipping trdiagonal %u < querylength %d\n",trdiagonal,querylength));
      } else {
	(*sense_trdiagonals)[k] = trdiagonal;
	(*sense_trnums)[k] = EF64_trnum(&(*sense_troffsets)[k],&(*sense_trhighs)[k],
					transcript_ef64,trdiagonal - querylength,trdiagonal);
	debug3(printf("sense trnum: %u  trdiagonal: %u\n",(*sense_trnums)[k],trdiagonal));
	k++;
      }
    }
    if ((*n_sense_trdiagonals = k) == 0) {
      FREE(*sense_trhighs);
      FREE(*sense_troffsets);
      FREE(*sense_trnums);
      *sense_trnums = (Trnum_T *) NULL;
      *sense_troffsets = (Trcoord_T *) NULL;
      *sense_trhighs = (Trcoord_T *) NULL;
    }
  }


  /* minus */
  nstreams = 0;
  stage1->tminus_stream_array[nstreams] = stage1->tr_minus_positions[0];
  stage1->tminus_streamsize_array[nstreams] = stage1->tr_minus_npositions[0];
  stage1->tminus_diagterm_array[nstreams++] = /*diagterm*/0 + index1part_tr;
  if (index1part_tr <= query_lastpos) {
    stage1->tminus_stream_array[nstreams] = stage1->tr_minus_positions[index1part_tr];
    stage1->tminus_streamsize_array[nstreams] = stage1->tr_minus_npositions[index1part_tr];
    stage1->tminus_diagterm_array[nstreams++] = /*diagterm*/(index1part_tr + index1part_tr);
  }
  
  stage1->tminus_stream_array[nstreams] = stage1->tr_minus_positions[query_lastpos];
  stage1->tminus_streamsize_array[nstreams] = stage1->tr_minus_npositions[query_lastpos];
  stage1->tminus_diagterm_array[nstreams++] = /*diagterm*/(query_lastpos + index1part_tr);
  if (query_lastpos - index1part_tr >= 0) {
    stage1->tminus_stream_array[nstreams] = stage1->tr_minus_positions[query_lastpos - index1part_tr];
    stage1->tminus_streamsize_array[nstreams] = stage1->tr_minus_npositions[query_lastpos - index1part_tr];
    stage1->tminus_diagterm_array[nstreams++] = /*diagterm*/(query_lastpos - index1part_tr) + index1part_tr;
  }

  *_antisense_trdiagonals = Merge_diagonals(&ndiagonals,
					   stage1->tminus_stream_array,
					   stage1->tminus_streamsize_array,
					   stage1->tminus_diagterm_array,
					   nstreams,stage1->mergeinfo_tr);

#if 0
  /* ? This looks for 30-mers at the ends */
  *n_antisense_trdiagonals = keep_duplicates(*_antisense_trdiagonals,ndiagonals);
#else
  *n_antisense_trdiagonals = make_unique(*_antisense_trdiagonals,ndiagonals);
#endif

  if ((*n_antisense_trdiagonals) == 0) {
    *antisense_trnums = (Trnum_T *) NULL;
    *antisense_troffsets = (Trcoord_T *) NULL;
    *antisense_trhighs = (Trcoord_T *) NULL;

  } else {
    *antisense_trnums = (Trnum_T *) MALLOC((*n_antisense_trdiagonals) * sizeof(Trnum_T));
    *antisense_troffsets = (Trcoord_T *) MALLOC((*n_antisense_trdiagonals) * sizeof(Trcoord_T));
    *antisense_trhighs = (Trcoord_T *) MALLOC((*n_antisense_trdiagonals) * sizeof(Trcoord_T));

    k = 0;
    for (i = 0; i < (*n_antisense_trdiagonals); i++) {
      trdiagonal = (*_antisense_trdiagonals)[i];
      debug3(printf("tminus trdiagonal: %u\n",trdiagonal));
      if (trdiagonal /*+qstart:0*/ < (Trcoord_T) querylength) {
	/* Skip */
	debug3(printf("Skipping trdiagonal %u < querylength %d\n",trdiagonal,querylength));
      } else {
	(*_antisense_trdiagonals)[k] = trdiagonal;
	(*antisense_trnums)[k] = EF64_trnum(&(*antisense_troffsets)[k],&(*antisense_trhighs)[k],
					    transcript_ef64,trdiagonal - querylength,trdiagonal);
	debug3(printf("antisense trnum: %u  trdiagonal: %u\n",(*antisense_trnums)[k],trdiagonal));
	k++;
      }
    }
    if ((*n_antisense_trdiagonals = k) == 0) {
      FREE(*antisense_trhighs);
      FREE(*antisense_troffsets);
      FREE(*antisense_trnums);
      *antisense_trnums = (Trnum_T *) NULL;
      *antisense_troffsets = (Trcoord_T *) NULL;
      *antisense_trhighs = (Trcoord_T *) NULL;
    }
  }


  debug3(printf("***Transcriptome_exact2: %d sense and %d antisense trdiagonals\n",
		*n_sense_trdiagonals,*n_antisense_trdiagonals));
	 
  return (*n_sense_trdiagonals) + (*n_antisense_trdiagonals);
}
#endif


#if 0
/* Want a pair with sufficient separation */
static int
check_indexpairs (bool **separatedp, int *indexpairs, int nindexpairs,
		  int querypos_j, int *accumulated_tmin, int *accumulated_tmax,
		  int min_separation) {
  int nseparated = 0;
  int i, k;
  int index2;

  *separatedp = (bool *) CALLOC(nindexpairs,sizeof(bool));
  for (i = 0, k = 0; i < nindexpairs; i++, k += 2) {
    index2 = indexpairs[k+1];
    debug2(printf("Comparing querypos_j %d with querypos_i %d..%d",
		  querypos_j,accumulated_tmin[index2],accumulated_tmax[index2]));
    if (querypos_j > /*querypos_i*/accumulated_tmin[index2] + min_separation) {
      debug2(printf(" => Sufficient\n"));
      (*separatedp)[i] = true;
      nseparated++;
    } else if (/*querypos_i*/accumulated_tmax[index2] > querypos_j + min_separation) {
      debug2(printf(" => Sufficient\n"));
      (*separatedp)[i] = true;
      nseparated++;
    } else {
      debug2(printf("\n"));
    }
  }

  if (nseparated == 0) {
    FREE(*separatedp);
  }

  return nseparated;
}
#endif


#if 0
static void
revise_accumulated (Trcoord_T **accumulated_diagonals, int **accumulated_tmin,
		    int **accumulated_tmax, int *accumulated_ndiagonals,
		    Trcoord_T *positions_j,
		    int npositions_j, int querypos_j, int diagterm_j) {

  Trcoord_T *new_diagonals, trdiagonal_i, trdiagonal_j;
  int *new_tmin, *new_tmax;
  int new_ndiagonals, i, j, k = 0;
  

  new_ndiagonals = (*accumulated_ndiagonals) + npositions_j; /* Includes duplicates */
  new_diagonals = (Trcoord_T *) MALLOC(new_ndiagonals * sizeof(Trcoord_T));
  new_tmin = (int *) MALLOC(new_ndiagonals * sizeof(int));
  new_tmax = (int *) MALLOC(new_ndiagonals * sizeof(int));

  i = j = 0;
  while (i < (*accumulated_ndiagonals) && j < npositions_j) {
    trdiagonal_i = (*accumulated_diagonals)[i];
    trdiagonal_j = positions_j[j] + diagterm_j;

    if (trdiagonal_i < trdiagonal_j) {
      new_diagonals[k] = trdiagonal_i;
      new_tmin[k] = (*accumulated_tmin)[i];
      new_tmax[k++] = (*accumulated_tmax)[i];
      i++;
    } else if (trdiagonal_j < trdiagonal_i) {
      new_diagonals[k] = trdiagonal_j;
      new_tmin[k] = querypos_j;
      new_tmax[k++] = querypos_j;
      j++;
    } else {
      new_diagonals[k] = trdiagonal_i;
      if (querypos_j < (*accumulated_tmin)[i]) {
	new_tmin[k] = querypos_j;
      } else {
	new_tmin[k] =  (*accumulated_tmin)[i];
      }
      if (querypos_j > (*accumulated_tmax)[i]) {
	new_tmax[k++] = querypos_j;
      } else {
	new_tmax[k++] =  (*accumulated_tmax)[i];
      }
      i++; j++;
    }
  }
	
  while (i < (*accumulated_ndiagonals)) {
    new_diagonals[k] = /*trdiagonal_i =*/ (*accumulated_diagonals)[i];
    new_tmin[k] = (*accumulated_tmin)[i];
    new_tmax[k++] = (*accumulated_tmax)[i];
    i++;
  }

  while (j < npositions_j) {
    new_diagonals[k] = /*trdiagonal_j =*/ positions_j[j] + diagterm_j;
    new_tmin[k] = querypos_j;
    new_tmax[k++] = querypos_j;
    j++;
  }
	
  FREE(*accumulated_tmax);
  FREE(*accumulated_tmin);
  FREE(*accumulated_diagonals);
  *accumulated_diagonals = new_diagonals;
  *accumulated_tmin = new_tmin;
  *accumulated_tmax = new_tmax;
  *accumulated_ndiagonals = k;	/* Excludes duplicates */


#ifdef DEBUG2
  for (k = 0; k < *accumulated_ndiagonals; k++) {
    printf("%d: %u %d..%d\n",
	   k,(*accumulated_diagonals)[k],(*accumulated_tmin)[k],(*accumulated_tmax)[k]);
    }
#endif

  return;
}
#endif


#if 0
/* Iterate until we find our first match */
/* _sense_trdiagonals and _antisense_trdiagonals need to be aligned */
int
Transcriptome_anypair (Trnum_T **sense_trnums, Trcoord_T **sense_troffsets, Trcoord_T **sense_trhighs,
		       Trcoord_T **_sense_trdiagonals, int **sense_tstarts, int **sense_tends, int *n_sense_trdiagonals,
		       Trnum_T **antisense_trnums, Trcoord_T **antisense_troffsets, Trcoord_T **antisense_trhighs,
		       Trcoord_T **_antisense_trdiagonals, int **antisense_tstarts, int **antisense_tends,
		       int *n_antisense_trdiagonals, Stage1_T stage1, int querylength) {

  int *indexpairs_tplus = NULL, *indexpairs_tminus = NULL;
  int nindexpairs_tplus = 0, nindexpairs_tminus = 0;

  Trcoord_T *accumulated_diagonals_tplus, *accumulated_diagonals_tminus;
  Trcoord_T *positions_j_tplus, *positions_j_tminus;
  Trcoord_T trdiagonal;

  int accumulated_ndiagonals_tplus, accumulated_ndiagonals_tminus,
    npositions_j_tplus, npositions_j_tminus, i, j, k_tplus, k_tminus;
  int *accumulated_tmin_tplus, *accumulated_tmax_tplus, *accumulated_tmin_tminus, *accumulated_tmax_tminus,
    tpos_j_tplus, diagterm_j_tplus, tpos_j_tminus, diagterm_j_tminus;

  int *order_tplus, *order_tminus;

  /* npositions should be filled from 0 through query_lastpos, inclusive */
  int len = querylength - index1part_tr + 1, index2, k;
  bool *separatedp_tplus = NULL, *separatedp_tminus = NULL, donep;
  int nseparated_tplus = 0, nseparated_tminus = 0;


  debug2(printf("Entered Transcriptome_anypair\n"));
  debug2(Stage1_dump_tr(stage1,querylength));

  /* These assertions should hold since we have called Stage1_trdiagonals_gc */
  /* assert(*_sense_trdiagonals == NULL); */
  /* assert(*_antisense_trdiagonals == NULL); */

  /* Initialize accumulated_ndiagonals, plus */
  /* Note: This will overwrite stage1->tr_plus_npositions[query_lastpos+1] */
  order_tplus = Sedgesort_order_int(stage1->tr_plus_npositions,len);

  accumulated_diagonals_tplus = (Trcoord_T *) NULL;
  accumulated_tmin_tplus = (int *) NULL;
  accumulated_tmax_tplus = (int *) NULL;
  accumulated_ndiagonals_tplus = 0;

  k_tplus = 0;
  while (k_tplus < len && accumulated_ndiagonals_tplus == 0) {
    tpos_j_tplus = order_tplus[k_tplus];
    if ((accumulated_ndiagonals_tplus = stage1->tr_plus_npositions[tpos_j_tplus]) == 0) {
      /* Skip */
      debug2(printf("Skipping initial plus tpos %d with no positions\n",tpos_j_tplus));

    } else {
      accumulated_diagonals_tplus = (Trcoord_T *) MALLOC(accumulated_ndiagonals_tplus * sizeof(Trcoord_T));
      accumulated_tmin_tplus = (int *) MALLOC(accumulated_ndiagonals_tplus * sizeof(int));
      accumulated_tmax_tplus = (int *) MALLOC(accumulated_ndiagonals_tplus * sizeof(int));
      positions_j_tplus = stage1->tr_plus_positions[tpos_j_tplus];
      diagterm_j_tplus = querylength - tpos_j_tplus; /* plus */
      for (j = 0; j < accumulated_ndiagonals_tplus; j++) {
	accumulated_diagonals_tplus[j] = positions_j_tplus[j] + diagterm_j_tplus;
	accumulated_tmin_tplus[j] = tpos_j_tplus;
	accumulated_tmax_tplus[j] = tpos_j_tplus;
      }
    }
    k_tplus++;
  }


  /* Initialize accumulated_ndiagonals, minus */
  /* Note: This will overwrite stage1->tr_minus_npositions[query_lastpos+1] */
  order_tminus = Sedgesort_order_int(stage1->tr_minus_npositions,len);

  accumulated_diagonals_tminus = (Trcoord_T *) NULL;
  accumulated_tmin_tminus = (int *) NULL;
  accumulated_tmax_tminus = (int *) NULL;
  accumulated_ndiagonals_tminus = 0;

  k_tminus = 0;
  while (k_tminus < len && accumulated_ndiagonals_tminus == 0) {
    tpos_j_tminus = order_tminus[k_tminus];
    if ((accumulated_ndiagonals_tminus = stage1->tr_minus_npositions[tpos_j_tminus]) == 0) {
      /* Skip */
      debug2(printf("Skipping initial minus tpos %d with no positions\n",tpos_j_tminus));

    } else {
      accumulated_diagonals_tminus = (Trcoord_T *) MALLOC(accumulated_ndiagonals_tminus * sizeof(Trcoord_T));
      accumulated_tmin_tminus = (int *) MALLOC(accumulated_ndiagonals_tminus * sizeof(int));
      accumulated_tmax_tminus = (int *) MALLOC(accumulated_ndiagonals_tminus * sizeof(int));
      positions_j_tminus = stage1->tr_minus_positions[tpos_j_tminus];
      diagterm_j_tminus = tpos_j_tminus + index1part_tr; /* minus */
      for (j = 0; j < accumulated_ndiagonals_tminus; j++) {
	accumulated_diagonals_tminus[j] = positions_j_tminus[j] + diagterm_j_tminus;
	accumulated_tmin_tminus[j] = tpos_j_tminus;
	accumulated_tmax_tminus[j] = tpos_j_tminus;
      }
    }
    k_tminus++;
  }

  /* Perform intersections with accumulated diagonals.  Search plus and minus in parallel. */
  donep = false;
  indexpairs_tplus = indexpairs_tminus = (int *) NULL;
  nindexpairs_tplus = nindexpairs_tminus = 0;
  while (donep == false && k_tplus < len && k_tminus < len) {
    /* plus */
    tpos_j_tplus = order_tplus[k_tplus];
    npositions_j_tplus = stage1->tr_plus_npositions[tpos_j_tplus]; /* Must be > 0 after initial loop */
    positions_j_tplus = stage1->tr_plus_positions[tpos_j_tplus];
    diagterm_j_tplus = querylength - tpos_j_tplus; /* plus */

    FREE(indexpairs_tplus);
    indexpairs_tplus = (int *) MALLOC(2 * npositions_j_tplus * sizeof(int));
    debug2(printf("Performing plus intersection of tpos %d (%d positions) with %d accumulated diagonals\n",
		  tpos_j_tplus,npositions_j_tplus,accumulated_ndiagonals_tplus));

    nseparated_tplus = 0;
    if ((nindexpairs_tplus = Intersect_indices_small(indexpairs_tplus,positions_j_tplus,npositions_j_tplus,
						     /*diagterm1*/diagterm_j_tplus,
						     accumulated_diagonals_tplus,accumulated_ndiagonals_tplus,
						     /*diagterm2*/0)) > 0 &&
	(nseparated_tplus = check_indexpairs(&separatedp_tplus,indexpairs_tplus,nindexpairs_tplus,tpos_j_tplus,
					     accumulated_tmin_tplus,accumulated_tmax_tplus,
					     /*min_separation*/index1part_tr)) > 0) {
      donep = true;

    } else {
      FREE(indexpairs_tplus);
    }
	

    /* minus */
    tpos_j_tminus = order_tminus[k_tminus];
    npositions_j_tminus = stage1->tr_minus_npositions[tpos_j_tminus]; /* Must be > 0 after initial loop */
    positions_j_tminus = stage1->tr_minus_positions[tpos_j_tminus];
    diagterm_j_tminus = tpos_j_tminus + index1part_tr; /* minus */

    FREE(indexpairs_tminus);
    indexpairs_tminus = (int *) MALLOC(2 * npositions_j_tminus * sizeof(int));
    debug2(printf("Performing minus intersection of tpos %d (%d positions) with %d accumulated diagonals\n",
		  tpos_j_tminus,npositions_j_tminus,accumulated_ndiagonals_tminus));

    nseparated_tminus = 0;
    if ((nindexpairs_tminus = Intersect_indices_small(indexpairs_tminus,positions_j_tminus,npositions_j_tminus,
						      /*diagterm1*/diagterm_j_tminus,
						      accumulated_diagonals_tminus,accumulated_ndiagonals_tminus,
						      /*diagterm2*/0)) > 0 &&
	(nseparated_tminus = check_indexpairs(&separatedp_tminus,indexpairs_tminus,nindexpairs_tminus,tpos_j_tminus,
					      accumulated_tmin_tminus,accumulated_tmax_tminus,
					      /*min_separation*/index1part_tr)) > 0) {
      donep = true;

    } else {
      FREE(indexpairs_tminus);
    }


    if (donep == false) {
      revise_accumulated(&accumulated_diagonals_tplus,&accumulated_tmin_tplus,&accumulated_tmax_tplus,
			 &accumulated_ndiagonals_tplus,
			 positions_j_tplus,npositions_j_tplus,tpos_j_tplus,diagterm_j_tplus);
      revise_accumulated(&accumulated_diagonals_tminus,&accumulated_tmin_tminus,&accumulated_tmax_tminus,
			 &accumulated_ndiagonals_tminus,
			 positions_j_tminus,npositions_j_tminus,
			 tpos_j_tminus,diagterm_j_tminus);
    }

    k_tplus++;
    k_tminus++;
  }


  while (donep == false && k_tplus < len) {
    /* plus */
    tpos_j_tplus = order_tplus[k_tplus];
    npositions_j_tplus = stage1->tr_plus_npositions[tpos_j_tplus]; /* Must be > 0 after initial loop */
    positions_j_tplus = stage1->tr_plus_positions[tpos_j_tplus];
    diagterm_j_tplus = querylength - tpos_j_tplus; /* plus */

    FREE(indexpairs_tplus);
    indexpairs_tplus = (int *) MALLOC(2 * npositions_j_tplus * sizeof(int));
    debug2(printf("Performing plus intersection of tpos %d (%d positions) with %d accumulated diagonals\n",
		  tpos_j_tplus,npositions_j_tplus,accumulated_ndiagonals_tplus));

    nseparated_tplus = 0;
    if ((nindexpairs_tplus = Intersect_indices_small(indexpairs_tplus,positions_j_tplus,npositions_j_tplus,
						     /*diagterm1*/diagterm_j_tplus,
						     accumulated_diagonals_tplus,accumulated_ndiagonals_tplus,
						     /*diagterm2*/0)) > 0 &&
	(nseparated_tplus = check_indexpairs(&separatedp_tplus,indexpairs_tplus,nindexpairs_tplus,tpos_j_tplus,
					     accumulated_tmin_tplus,accumulated_tmax_tplus,
					     /*min_separation*/index1part_tr)) > 0) {
      donep = true;

    } else {
      FREE(indexpairs_tplus);
      revise_accumulated(&accumulated_diagonals_tplus,&accumulated_tmin_tplus,
			 &accumulated_tmax_tplus,&accumulated_ndiagonals_tplus,
			 positions_j_tplus,npositions_j_tplus,
			 tpos_j_tplus,diagterm_j_tplus);
    }

    k_tplus++;
  }


  while (donep == false && k_tminus < len) {
    /* minus */
    tpos_j_tminus = order_tminus[k_tminus];
    npositions_j_tminus = stage1->tr_minus_npositions[tpos_j_tminus]; /* Must be > 0 after initial loop */
    positions_j_tminus = stage1->tr_minus_positions[tpos_j_tminus];
    diagterm_j_tminus = tpos_j_tminus + index1part_tr; /* minus */

    FREE(indexpairs_tminus);
    indexpairs_tminus = (int *) MALLOC(2 * npositions_j_tminus * sizeof(int));
    debug2(printf("Performing minus intersection of tpos %d (%d positions) with %d accumulated diagonals\n",
		  tpos_j_tminus,npositions_j_tminus,accumulated_ndiagonals_tminus));

    nseparated_tminus = 0;
    if ((nindexpairs_tminus = Intersect_indices_small(indexpairs_tminus,positions_j_tminus,npositions_j_tminus,
						      /*diagterm1*/diagterm_j_tminus,
						      accumulated_diagonals_tminus,accumulated_ndiagonals_tminus,
						      /*diagterm2*/0)) > 0 &&
	(nseparated_tminus = check_indexpairs(&separatedp_tminus,indexpairs_tminus,nindexpairs_tminus,tpos_j_tminus,
					      accumulated_tmin_tminus,accumulated_tmax_tminus,
					      /*min_separation*/index1part_tr)) > 0) {
      donep = true;

    } else {
      FREE(indexpairs_tminus);
      revise_accumulated(&accumulated_diagonals_tminus,&accumulated_tmin_tminus,
			 &accumulated_tmax_tminus,&accumulated_ndiagonals_tminus,
			 positions_j_tminus,npositions_j_tminus,
			 tpos_j_tminus,diagterm_j_tminus);
    }

    k_tminus++;
  }

  FREE(order_tminus);
  FREE(order_tplus);


  /* Prepare output, tplus */
  if (nseparated_tplus == 0) {
    *sense_trnums = (Trnum_T *) NULL;
    *sense_troffsets = (Trcoord_T *) NULL;
    *sense_trhighs = (Trcoord_T *) NULL;
    *_sense_trdiagonals = (Trcoord_T *) NULL;
    *sense_tstarts = (int *) NULL;
    *sense_tends = (int *) NULL;
    *n_sense_trdiagonals = 0;

  } else {
    *n_sense_trdiagonals = nseparated_tplus;
    *sense_trnums = (Trnum_T *) MALLOC(nseparated_tplus * sizeof(int));
    *sense_troffsets = (Trcoord_T *) MALLOC(nseparated_tplus * sizeof(Trcoord_T));
    *sense_trhighs = (Trcoord_T *) MALLOC(nseparated_tplus * sizeof(Trcoord_T));
    *_sense_trdiagonals = (Trcoord_T *) MALLOC_ALIGN(nseparated_tplus * sizeof(Trcoord_T));
    *sense_tstarts = (int *) MALLOC(nseparated_tplus * sizeof(int));
    *sense_tends = (int *) MALLOC(nseparated_tplus * sizeof(int));
    
    j = 0;
    for (i = 0, k = 0; i < nindexpairs_tplus; i++, k += 2) {
      if (separatedp_tplus[i] == true) {
	index2 = indexpairs_tplus[k+1];
	trdiagonal = accumulated_diagonals_tplus[index2];
	if (trdiagonal /*+tstart:0*/ < (Trcoord_T) querylength) {
	  /* Skip */
	} else {
	  (*_sense_trdiagonals)[j] = trdiagonal;
	  (*sense_trnums)[j] = EF64_trnum(&(*sense_troffsets)[j],&(*sense_trhighs)[j],
					  transcript_ef64,trdiagonal - querylength,trdiagonal);
	  (*sense_tstarts)[j] = accumulated_tmin_tplus[index2];
	  (*sense_tends)[j] = accumulated_tmax_tplus[index2] + index1part_tr;
	  j++;
	}
      }
    }
    FREE(separatedp_tplus);

    if ((*n_sense_trdiagonals = j) == 0) {
      FREE(*sense_tends);
      FREE(*sense_tstarts);
      FREE_ALIGN(*_sense_trdiagonals);
      FREE(*sense_trhighs);
      FREE(*sense_troffsets);
      FREE(*sense_trnums);
      *sense_trnums = (Trnum_T *) NULL;
      *sense_troffsets = (Trcoord_T *) NULL;
      *sense_trhighs = (Trcoord_T *) NULL;
      *_sense_trdiagonals = (Trcoord_T *) NULL;
      *sense_tstarts = (int *) NULL;
      *sense_tends = (int *) NULL;
    }
  }
  FREE(indexpairs_tplus);

  FREE(accumulated_diagonals_tplus);
  FREE(accumulated_tmin_tplus);
  FREE(accumulated_tmax_tplus);


  /* Prepare output, tminus */
  if (nseparated_tminus == 0) {
    *antisense_trnums = (Trnum_T *) NULL;
    *antisense_troffsets = (Trcoord_T *) NULL;
    *antisense_trhighs = (Trcoord_T *) NULL;
    *_antisense_trdiagonals = (Trcoord_T *) NULL;
    *antisense_tstarts = (int *) NULL;
    *antisense_tends = (int *) NULL;
    *n_antisense_trdiagonals = 0;

  } else {
    *n_antisense_trdiagonals = nseparated_tminus;
    *antisense_trnums = (Trnum_T *) MALLOC(nseparated_tminus * sizeof(Trnum_T));
    *antisense_troffsets = (Trcoord_T *) MALLOC(nseparated_tminus * sizeof(Trcoord_T));
    *antisense_trhighs = (Trcoord_T *) MALLOC(nseparated_tminus * sizeof(Trcoord_T));
    *_antisense_trdiagonals = (Trcoord_T *) MALLOC_ALIGN(nseparated_tminus * sizeof(Trcoord_T));
    *antisense_tstarts = (int *) MALLOC(nseparated_tminus * sizeof(int));
    *antisense_tends = (int *) MALLOC(nseparated_tminus * sizeof(int));
    
    j = 0;
    for (i = 0, k = 0; i < nindexpairs_tminus; i++, k += 2) {
      if (separatedp_tminus[i] == true) {
	index2 = indexpairs_tminus[k+1];
	trdiagonal = accumulated_diagonals_tminus[index2];
	if (trdiagonal /*+tstart:0*/ < (Trcoord_T) querylength) {
	  /* Skip */
	} else {
	  (*_antisense_trdiagonals)[j] = trdiagonal;
	  (*antisense_trnums)[j] = EF64_trnum(&(*antisense_troffsets)[j],&(*antisense_trhighs)[j],
					      transcript_ef64,trdiagonal - querylength,trdiagonal);
	  (*antisense_tstarts)[j] = querylength - (accumulated_tmax_tminus[index2] + index1part_tr);
	  (*antisense_tends)[j] = querylength - accumulated_tmin_tminus[index2];
	  j++;
	}
      }
    }
    FREE(separatedp_tminus);

    if ((*n_antisense_trdiagonals = j) == 0) {
      FREE(*antisense_tends);
      FREE(*antisense_tstarts);
      FREE_ALIGN(*_antisense_trdiagonals);
      FREE(*antisense_trhighs);
      FREE(*antisense_troffsets);
      FREE(*antisense_trnums);
      *antisense_trnums = (Trnum_T *) NULL;
      *antisense_troffsets = (Trcoord_T *) NULL;
      *antisense_trhighs = (Trcoord_T *) NULL;
      *_antisense_trdiagonals = (Trcoord_T *) NULL;
      *antisense_tstarts = (int *) NULL;
      *antisense_tends = (int *) NULL;
    }
  }
  FREE(indexpairs_tminus);

  FREE(accumulated_diagonals_tminus);
  FREE(accumulated_tmin_tminus);
  FREE(accumulated_tmax_tminus);
  
  debug2(printf("Transcriptome_anypair returning %d sense and %d antisense hits\n",
		(*n_sense_trdiagonals),(*n_antisense_trdiagonals)));

#ifdef DEBUG2
  for (k = 0; k < *n_sense_trdiagonals; k++) {
    printf("sense %d: %u %d..%d\n",
	   k,(*_sense_trdiagonals)[k],(*sense_tstarts)[k],(*sense_tends)[k]);
  }
  for (k = 0; k < *n_antisense_trdiagonals; k++) {
    printf("antisense %d: %u %d..%d\n",
	   k,(*_antisense_trdiagonals)[k],(*antisense_tstarts)[k],(*antisense_tends)[k]);
  }
#endif

  return (*n_sense_trdiagonals) + (*n_antisense_trdiagonals);
}
#endif


static int
most_prevalent (Trcoord_T *diagonals, int ndiagonals) {
  int most_prevalent = 0;
  int i, j;

  i = 0;
  while (i < ndiagonals) {
    j = i + 1;
    while (j < ndiagonals && diagonals[j] == diagonals[i]) {
      j++;
    }
    if (j - i > most_prevalent) {
      most_prevalent = j - i;
    }
    
    i = j;
  }

  debug9(printf("most_prevalent found a trdiagonal with %d kmers\n",most_prevalent));
  return most_prevalent;
}


static Trcoord_T *
keep_prevalence (int *ndiagonals, Trcoord_T *all_diagonals, int all_ndiagonals,
		int prevalence) {
  Trcoord_T *diagonals;
  int k, i, j;

  /* Count number of diagonals with given prevalence */
  k = 0;
  i = 0;
  while (i < all_ndiagonals) {
    j = i + 1;
    while (j < all_ndiagonals && all_diagonals[j] == all_diagonals[i]) {
      j++;
    }
    if (j - i >= prevalence) {
      k++;
    }
    
    i = j;
  }

  debug9(printf("keep_prevalence found %d diagonals with %d kmers\n",k,prevalence));

  /* Allocate new array.  Needs to be aligned for later merging with Stage1_T all_trdiagonals */
  if (k == 0) {
    *ndiagonals = 0;
    return (Trcoord_T *) NULL;

  } else {
    diagonals = (Trcoord_T *) MALLOC_ALIGN(k*sizeof(Trcoord_T));

    /* Save diagonals if there are prevalence or more of them */
    k = 0;
    i = 0;
    while (i < all_ndiagonals) {
      j = i + 1;
      while (j < all_ndiagonals && all_diagonals[j] == all_diagonals[i]) {
	j++;
      }
      if (j - i >= prevalence) {
	diagonals[k++] = all_diagonals[i];
      }
      
      i = j;
    }

    *ndiagonals = k;
    return diagonals;
  }
}


static Trcoord_T *
get_trdiagonals_sense (int *ndiagonals, int *max_npositions,
		       Stage1_T this, int querylength) {
  int querypos, query_lastpos;

  query_lastpos = querylength - index1part_tr;

  *max_npositions = 0;
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    this->tr_plus_diagterms[querypos] = querylength - querypos; /* plus */
    if (this->tr_plus_npositions[querypos] > *max_npositions) {
      *max_npositions = this->tr_plus_npositions[querypos];
    }
  }
  
  return Merge_diagonals(&(*ndiagonals),this->tr_plus_positions,this->tr_plus_npositions,
			 this->tr_plus_diagterms,/*nstreams*/query_lastpos+1,this->mergeinfo_tr);
}


static Trcoord_T *
get_trdiagonals_antisense (int *ndiagonals, int *max_npositions,
			Stage1_T this, int querylength) {
  int querypos, query_lastpos;

  query_lastpos = querylength - index1part_tr;

  *max_npositions = 0;
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    this->tr_minus_diagterms[querypos] = querypos + index1part_tr; /* minus */
    if (this->tr_minus_npositions[querypos] > *max_npositions) {
      *max_npositions = this->tr_minus_npositions[querypos];
    }
  }
  
  return Merge_diagonals(&(*ndiagonals),this->tr_minus_positions,this->tr_minus_npositions,
			 this->tr_minus_diagterms,/*nstreams*/query_lastpos+1,this->mergeinfo_tr);
}


/* Results are aligned because of Merge_diagonals */
/* Essentially a streamlined version of Kmer_search_complete */
/* Saves all_trdiagonals for later calls to Kmer_widest or Kmer_anchored */
int
Transcriptome_prevalent (Trnum_T **sense_trnums, Trcoord_T **sense_troffsets, Trcoord_T **sense_trhighs,
			 Trcoord_T **_sense_trdiagonals, int **sense_tstarts, int **sense_tends, int *n_sense_trdiagonals,
			 Trnum_T **antisense_trnums, Trcoord_T **antisense_troffsets, Trcoord_T **antisense_trhighs,
			 Trcoord_T **_antisense_trdiagonals, int **antisense_tstarts, int **antisense_tends,
			 int *n_antisense_trdiagonals,
			 Stage1_T stage1, int querylength, Compress_T query_compress_fwd, Compress_T query_compress_rev) {

  Trcoord_T *_all_trdiagonals_sense, *_all_trdiagonals_antisense, trdiagonal;
  int max_npositions_sense, max_npositions_antisense;

  int plus_ntrdiagonals, minus_ntrdiagonals, plus_nprevalent, minus_nprevalent, i, k;
  int tstart, tend, nmismatches_to_trimpos;


  debug9(printf("Entered Transcriptome_prevalent\n"));
  debug9(Stage1_dump_tr(stage1,querylength));

  /* all_trdiagonals has duplicates so we can find most prevalent */
  /* tplus or sense */
  _all_trdiagonals_sense =
    get_trdiagonals_sense(&plus_ntrdiagonals,&max_npositions_sense,
			  stage1,querylength);
  plus_nprevalent = most_prevalent(_all_trdiagonals_sense,plus_ntrdiagonals);
  /* plus_nprevalent -= 1; -- ? Need to allow for suboptimal */

  /* tminus or antisense */
  _all_trdiagonals_antisense =
    get_trdiagonals_antisense(&minus_ntrdiagonals,&max_npositions_antisense,
			      stage1,querylength);
  minus_nprevalent = most_prevalent(_all_trdiagonals_antisense,minus_ntrdiagonals);
  /* minus_nprevalent -= 1; -- ? Need to allow for suboptimal*/


  /* Pick best direction */
  if (plus_nprevalent > minus_nprevalent + index1interval_tr) {
    *_sense_trdiagonals = keep_prevalence(&(*n_sense_trdiagonals),
					 _all_trdiagonals_sense,plus_ntrdiagonals,
					 plus_nprevalent - (index1interval_tr - 1));
    *_antisense_trdiagonals = (Trcoord_T *) NULL;
    *n_antisense_trdiagonals = 0;

  } else if (minus_nprevalent > plus_nprevalent + index1interval_tr) {
    *_antisense_trdiagonals = keep_prevalence(&(*n_antisense_trdiagonals),
					     _all_trdiagonals_antisense,minus_ntrdiagonals,
					     minus_nprevalent - (index1interval_tr - 1));
    *_sense_trdiagonals = (Trcoord_T *) NULL;
    *n_sense_trdiagonals = 0;

  } else {
    *_sense_trdiagonals = keep_prevalence(&(*n_sense_trdiagonals),
					 _all_trdiagonals_sense,plus_ntrdiagonals,
					 plus_nprevalent - (index1interval_tr - 1));
    *_antisense_trdiagonals = keep_prevalence(&(*n_antisense_trdiagonals),
					     _all_trdiagonals_antisense,minus_ntrdiagonals,
					     minus_nprevalent - (index1interval_tr - 1));
  }
    
  /* all_ntrdiagonals_sense = make_unique(_all_trdiagonals_sense,plus_ntrdiagonals); */
  /* all_ntrdiagonals_antisense = make_unique(_all_trdiagonals_antisense,minus_ntrdiagonals); */
  FREE_ALIGN(_all_trdiagonals_antisense);
  FREE_ALIGN(_all_trdiagonals_sense);


  if (*n_sense_trdiagonals == 0) {
    *sense_trnums = (Trnum_T *) NULL;
    *sense_troffsets = (Trcoord_T *) NULL;
    *sense_trhighs = (Trcoord_T *) NULL;
    *sense_tstarts = (int *) NULL;
    *sense_tends = (int *) NULL;

  } else {
    *sense_trnums = (Trnum_T *) MALLOC((*n_sense_trdiagonals) * sizeof(Trnum_T));
    *sense_troffsets = (Trcoord_T *) MALLOC((*n_sense_trdiagonals) * sizeof(Trcoord_T));
    *sense_trhighs = (Trcoord_T *) MALLOC((*n_sense_trdiagonals) * sizeof(Trcoord_T));
    *sense_tstarts = (int *) MALLOC((*n_sense_trdiagonals) * sizeof(int));
    *sense_tends = (int *) MALLOC((*n_sense_trdiagonals) * sizeof(int));

    k = 0;
    for (i = 0; i < *n_sense_trdiagonals; i++) {
      trdiagonal = (*_sense_trdiagonals)[i];
      if (trdiagonal /*+qstart:0*/ < (Trcoord_T) querylength) {
	/* Skip */
      } else {
	tstart = Genomebits_trim_qstart(&nmismatches_to_trimpos,/*query_compress*/query_compress_fwd,
					transcriptomebits,(Univcoord_T) trdiagonal,querylength,
					/*pos5*/0,/*pos3*/querylength,/*plusp*/true,/*genestrand*/0);
	tend = Genomebits_trim_qend(&nmismatches_to_trimpos,/*query_compress*/query_compress_fwd,
				    transcriptomebits,(Univcoord_T) trdiagonal,querylength,
				    /*pos5*/0,/*pos3*/querylength,/*plusp*/true,/*genestrand*/0);
	if (tstart < tend) {
	  (*_sense_trdiagonals)[k] = trdiagonal;
	  (*sense_trnums)[k] = EF64_trnum(&(*sense_troffsets)[k],&(*sense_trhighs)[k],
					transcript_ef64,trdiagonal - querylength,trdiagonal);
	  (*sense_tstarts)[k] = tstart;
	  (*sense_tends)[k] = tend;
	  debug3(printf("sense trnum: %u  trdiagonal: %u, trims: %d..%d\n",
			(*sense_trnums)[k],trdiagonal,(*sense_tstarts)[k],(*sense_tends)[k]));
	  k++;
	}
      }
    }

    if ((*n_sense_trdiagonals = k) == 0) {
      FREE(*sense_tends);
      FREE(*sense_tstarts);
      FREE(*sense_trhighs);
      FREE(*sense_troffsets);
      FREE(*sense_trnums);

      *sense_trnums = (Trnum_T *) NULL;
      *sense_troffsets = (Trcoord_T *) NULL;
      *sense_trhighs = (Trcoord_T *) NULL;
      *sense_tstarts = (int *) NULL;
      *sense_tends = (int *) NULL;
    }
  }


  if (*n_antisense_trdiagonals == 0) {
    *antisense_trnums = (Trnum_T *) NULL;
    *antisense_troffsets = (Trcoord_T *) NULL;
    *antisense_trhighs = (Trcoord_T *) NULL;
    *antisense_tstarts = (int *) NULL;
    *antisense_tends = (int *) NULL;

  } else {
    *antisense_trnums = (Trnum_T *) MALLOC((*n_antisense_trdiagonals) * sizeof(Trnum_T));
    *antisense_troffsets = (Trcoord_T *) MALLOC((*n_antisense_trdiagonals) * sizeof(Trcoord_T));
    *antisense_trhighs = (Trcoord_T *) MALLOC((*n_antisense_trdiagonals) * sizeof(Trcoord_T));
    *antisense_tstarts = (int *) MALLOC((*n_antisense_trdiagonals) * sizeof(int));
    *antisense_tends = (int *) MALLOC((*n_antisense_trdiagonals) * sizeof(int));

    k = 0;
    for (i = 0; i < *n_antisense_trdiagonals; i++) {
      trdiagonal = (*_antisense_trdiagonals)[i];
      if (trdiagonal /*+qstart:0*/ < (Trcoord_T) querylength) {
	/* Skip */
      } else {
	tstart = Genomebits_trim_qstart(&nmismatches_to_trimpos,/*query_compress*/query_compress_rev,
					transcriptomebits,(Univcoord_T) trdiagonal,querylength,
					/*pos5*/0,/*pos3*/querylength,/*plusp*/false,/*genestrand*/0);
	tend = Genomebits_trim_qend(&nmismatches_to_trimpos,/*query_compress*/query_compress_rev,
				    transcriptomebits,(Univcoord_T) trdiagonal,querylength,
				    /*pos5*/0,/*pos3*/querylength,/*plusp*/false,/*genestrand*/0);

	if (tstart < tend) {
	  (*_antisense_trdiagonals)[k] = trdiagonal;
	  (*antisense_trnums)[k] = EF64_trnum(&(*antisense_troffsets)[k],&(*antisense_trhighs)[k],
					      transcript_ef64,trdiagonal - querylength,trdiagonal);
	  (*antisense_tstarts)[k] = tstart;
	  (*antisense_tends)[k] = tend;
	  debug3(printf("antisense trnum: %u  trdiagonal: %u, trims: %d..%d\n",
			(*antisense_trnums)[k],trdiagonal,(*antisense_tstarts)[k],(*antisense_tends)[k]));
	  k++;
	}
      }
    }

    if ((*n_antisense_trdiagonals = k) == 0) {
      FREE(*antisense_tends);
      FREE(*antisense_tstarts);
      FREE(*antisense_trhighs);
      FREE(*antisense_troffsets);
      FREE(*antisense_trnums);
      *antisense_trnums = (Trnum_T *) NULL;
      *antisense_troffsets = (Trcoord_T *) NULL;
      *antisense_trhighs = (Trcoord_T *) NULL;
      *antisense_tstarts = (int *) NULL;
      *antisense_tends = (int *) NULL;
    }
  }

  return (*n_sense_trdiagonals) + (*n_antisense_trdiagonals);
}


#if 0
/* Creates a trpath from trdiagonals at each end.  Should handle the
   case where exact search fails due to a mismatch in one of the end
   k-mers */
void
Transcriptome_search_end (int *found_score,

			  List_T *sense_trpaths, List_T *antisense_trpaths,
				  
			  Stage1_T stage1, int querylength,
			  Compress_T query_compress_fwd, Compress_T query_compress_rev,
			  Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			  Trpathpool_T trpathpool, Hitlistpool_T hitlistpool, Method_T method) {
  Trpath_T trpath;
  
  Trcoord_T *tplus_positions_end5 = NULL, *tminus_positions_end5 = NULL,
    *tplus_positions_end3 = NULL, *tminus_positions_end3 = NULL;
  int n_tplus_positions_end5 = 0, n_tminus_positions_end5 = 0,
    n_tplus_positions_end3 = 0, n_tminus_positions_end3 = 0;
  int tplus_diagterm_end5, tminus_diagterm_end5, tplus_diagterm_end3, tminus_diagterm_end3;

  int i;
  Trcoord_T trdiagonal;

  int query_lastpos;



  debug5(printf("\n\n***Transcriptome_search_end\n"));

  query_lastpos = querylength - index1part_tr;
  if (stage1->tr_validp[0] == false || stage1->tr_validp[query_lastpos] == false) {
    return;
  } else {
    tplus_positions_end5 = stage1->tr_plus_positions[0];
    n_tplus_positions_end5 = stage1->tr_plus_npositions[0];
    tplus_diagterm_end5 = querylength /* - 0, querypos */;

    tminus_positions_end5 = stage1->tr_minus_positions[0];
    n_tminus_positions_end5 = stage1->tr_minus_npositions[0];
    tminus_diagterm_end5 = /* 0, querypos + */ index1part_tr;

    tplus_positions_end3 = stage1->tr_plus_positions[query_lastpos];
    n_tplus_positions_end3 = stage1->tr_plus_npositions[query_lastpos];
    tplus_diagterm_end3 = index1part_tr /* querylength - querypos */;

    tminus_positions_end3 = stage1->tr_minus_positions[query_lastpos];
    n_tminus_positions_end3 = stage1->tr_minus_npositions[query_lastpos];
    tminus_diagterm_end3 = querylength /* querypos + index1part_tr */;
  }

  debug5(printf("***Transcriptome end: plus 5' %d\n",n_tplus_positions_end5));
  for (i = 0; i < n_tplus_positions_end5; i++) {
    trdiagonal = tplus_positions_end5[i] + tplus_diagterm_end5;
    debug5(printf("tplus trdiagonal: %u\n",trdiagonal));
    if (trdiagonal /*+qstart:0*/ < (Trcoord_T) querylength) {
      /* Skip */
    } else {
      if ((trpath = Trpath_solve_from_trstart(trdiagonal,/*tplusp*/true,
					      querylength,/*query_compress_tr*/query_compress_fwd,
					      intlistpool,uintlistpool,
					      trpathpool,method)) != NULL) {
	if (trpath->found_score < *found_score) {
	  *found_score = trpath->found_score;
	}
	*sense_trpaths = Hitlist_push(*sense_trpaths,hitlistpool,(void *) trpath
				      hitlistpool_trace(__FILE__,__LINE__));
      }
    }
  }
    

  debug5(printf("***Transcriptome end: plus 3' %d trdiagonals\n",n_tplus_positions_end3));
  for (i = 0; i < n_tplus_positions_end3; i++) {
    trdiagonal = tplus_positions_end3[i] + tplus_diagterm_end3;
    debug5(printf("tplus trdiagonal: %u\n",trdiagonal));
    if (trdiagonal /*+qtart:0*/ < (Trcoord_T) querylength) {
      /* Skip */
    } else {
      if ((trpath = Trpath_solve_from_trend(trdiagonal,/*tplusp*/true,
					    querylength,/*query_compress_tr*/query_compress_fwd,
					    intlistpool,uintlistpool,
					    trpathpool,method)) != NULL) {
	if (trpath->found_score < *found_score) {
	  *found_score = trpath->found_score;
	}
	*sense_trpaths = Hitlist_push(*sense_trpaths,hitlistpool,(void *) trpath
				      hitlistpool_trace(__FILE__,__LINE__));
      }
    }
  }


  debug5(printf("***Transcriptome end: minus 5' %dtrdiagonals\n",n_tminus_positions_end5));
  for (i = 0; i < n_tminus_positions_end5; i++) {
    trdiagonal = tminus_positions_end5[i] + tminus_diagterm_end5;
    debug5(printf("tminus trdiagonal: %u\n",trdiagonal));
    if (trdiagonal /*+qstart:0*/ < (Trcoord_T) querylength) {
      /* Skip */
    } else {
      if ((trpath = Trpath_solve_from_trend(trdiagonal,/*tplusp*/false,
					    querylength,/*query_compress_tr*/query_compress_rev,
					    intlistpool,uintlistpool,
					    trpathpool,method)) != NULL) {
	if (trpath->found_score < *found_score) {
	  *found_score = trpath->found_score;
	}
	*antisense_trpaths = Hitlist_push(*antisense_trpaths,hitlistpool,(void *) trpath
					  hitlistpool_trace(__FILE__,__LINE__));
      }
    }
  }


  debug5(printf("***Transcriptome end: minus 3' %d trdiagonals\n",n_tminus_positions_end3));
  for (i = 0; i < n_tminus_positions_end3; i++) {
    trdiagonal = tminus_positions_end3[i] + tminus_diagterm_end3;
    debug5(printf("tminus trdiagonal: %u\n",trdiagonal));
    if (trdiagonal /*+qstart:0*/ < (Trcoord_T) querylength) {
      /* Skip */
    } else {
      if ((trpath = Trpath_solve_from_trstart(trdiagonal,/*tplusp*/false,
					      querylength,/*query_compress_tr*/query_compress_rev,
					      intlistpool,uintlistpool,
					      trpathpool,method)) != NULL) {
	if (trpath->found_score < *found_score) {
	  *found_score = trpath->found_score;
	}
	*antisense_trpaths = Hitlist_push(*antisense_trpaths,hitlistpool,(void *) trpath
					  hitlistpool_trace(__FILE__,__LINE__));
      }
    }
  }

  return;
}
#endif


#if 0
/* Intersect_approx does not seem to reduce speed and can find some
   alignments that Intersect_small cannot.  However, this was before
   we implemented a SIMD version of Intersect_small. */
void
Transcriptome_search_approx (int *found_score,

			     List_T *sense_trpaths, List_T *antisense_trpaths,
				  
			     Stage1_T stage1, int querylength,

			     Compress_T query_compress_fwd, Compress_T query_compress_rev,
			     int max_insertionlen, int max_deletionlen,
			     Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			     Listpool_T listpool, Trpathpool_T trpathpool, Pathpool_T pathpool,
			     Hitlistpool_T hitlistpool, Method_T method) {
  Trpath_T trpath;
  
  Trcoord_T *tplus_positions_end5 = NULL, *tminus_positions_end5 = NULL,
    *tplus_positions_end3 = NULL, *tminus_positions_end3 = NULL;
  int n_tplus_positions_end5 = 0, n_tminus_positions_end5 = 0,
    n_tplus_positions_end3 = 0, n_tminus_positions_end3 = 0;
  int tplus_diagterm_end5, tminus_diagterm_end5, tplus_diagterm_end3, tminus_diagterm_end3;

  Trcoord_T trdiagonala, trdiagonalb;
  Trcoord_T min_distance, distance;

  Trcoord_T *tplus_diagpairs, *tminus_diagpairs;
  int n_tplus_diagpairs, n_tminus_diagpairs, k, j;
#ifdef INCLUDE_EXACT
  bool tplus_exactp, tminus_exactp;
#endif
  int pos5a, pos3a, pos5b, pos3b;
  int nmismatches_a5, nmismatches_a3, nmismatches_b5, nmismatches_b3;

  int query_lastpos;

  Trnum_T trnum;
  Trcoord_T troffset, trhigh;


#if 0
  *tplus_positions_5 = *tminus_positions_5 = *tplus_positions_3 = *tminus_positions_3 = (Trcoord_T *) NULL;
  *n_tplus_positions_5 = *n_tminus_positions_5 = *n_tplus_positions_3 = *n_tminus_positions_3 = 0;
#endif
			
  debug4(printf("\n\n***Transcriptome_search_approx\n"));

  /* Stage1_dump_tr(stage1,querylength); */

  query_lastpos = querylength - index1part_tr;
  if (stage1->tr_validp[0] == false || stage1->tr_validp[query_lastpos] == false) {
    return;
  } else {
    tplus_positions_end5 = stage1->tr_plus_positions[0];
    n_tplus_positions_end5 = stage1->tr_plus_npositions[0];
    tplus_diagterm_end5 = querylength /* - 0, querypos */;

    tminus_positions_end5 = stage1->tr_minus_positions[0];
    n_tminus_positions_end5 = stage1->tr_minus_npositions[0];
    tminus_diagterm_end5 = /* 0, querypos + */ index1part_tr;

    tplus_positions_end3 = stage1->tr_plus_positions[query_lastpos];
    n_tplus_positions_end3 = stage1->tr_plus_npositions[query_lastpos];
    tplus_diagterm_end3 = index1part_tr /* querylength - querypos */;

    tminus_positions_end3 = stage1->tr_minus_positions[query_lastpos];
    n_tminus_positions_end3 = stage1->tr_minus_npositions[query_lastpos];
    tminus_diagterm_end3 = querylength /* querypos + index1part_tr */;
  }

  /* Output of Intersect_approx is to have rare first, and frequent second */
  tplus_diagpairs = Intersect_approx_uint4(&n_tplus_diagpairs,
					   /*set1*/tplus_positions_end5,/*length1*/n_tplus_positions_end5,/*diagterm1*/tplus_diagterm_end5,
					   /*set2*/tplus_positions_end3,/*length2*/n_tplus_positions_end3,/*diagterm2*/tplus_diagterm_end3,
					   /*below_slop*/max_insertionlen,/*above_slop*/max_deletionlen);
  tminus_diagpairs = Intersect_approx_uint4(&n_tminus_diagpairs,
					    /*set1*/tminus_positions_end5,/*length1*/n_tminus_positions_end5,/*diagterm1*/tminus_diagterm_end5,
					    /*set2*/tminus_positions_end3,/*length2*/n_tminus_positions_end3,/*diagterm2*/tminus_diagterm_end3,
					    /*below_slop*/max_deletionlen,/*above_slop*/max_insertionlen);

  debug4(printf("***Transcriptome approx: %d plus and %d minus diagpairs\n",
		n_tplus_diagpairs,n_tminus_diagpairs));

  k = 0;
  while (k < 2*n_tplus_diagpairs) {
    trdiagonala = tplus_diagpairs[k];
    trdiagonalb = tplus_diagpairs[k+1];
    if (trdiagonalb <= trdiagonala) {
      min_distance = trdiagonala - trdiagonalb;
    } else {
      min_distance = trdiagonalb - trdiagonala;
    }

    /* Pick the closest pair */
    j = k + 2;
    while (j < n_tplus_diagpairs && tplus_diagpairs[j] == trdiagonala) {
      if (tplus_diagpairs[j+1] <= trdiagonala) {
	distance = trdiagonala - tplus_diagpairs[j+1];
      } else {
	distance = tplus_diagpairs[j+1] - trdiagonala;
      }
      if (distance < min_distance) {
	min_distance = distance;
	trdiagonalb = tplus_diagpairs[j+1];
      }
      j += 2;
    }

    debug4(printf("tplus diagpairs: %u and %u\n",trdiagonala,trdiagonalb));
    if (trdiagonalb == trdiagonala) {
      /* Skip exact results, which were covered by TR_EXACT */
    } else {
      pos5a = Genomebits_first_kmer_left(&nmismatches_a5,transcriptomebits,query_compress_fwd,
					 /*univdiagonal*/(Univcoord_T) trdiagonala,querylength,
					 /*pos5*/0,/*pos3*/querylength,/*plusp*/true,/*genestrand*/0,
					 /*query_unk_mismatch_local_p*/false,/*kmer*/index1part_tr);
      pos3a = Genomebits_first_kmer_right(&nmismatches_a3,transcriptomebits,query_compress_fwd,
					  /*univdiagonal*/(Univcoord_T) trdiagonala,querylength,
					  /*pos5*/0,/*pos3*/querylength,/*plusp*/true,/*genestrand*/0,
					  /*query_unk_mismatch_local_p*/false,/*kmer*/index1part_tr);

      pos5b = Genomebits_first_kmer_left(&nmismatches_b5,transcriptomebits,query_compress_fwd,
					 /*univdiagonal*/(Univcoord_T) trdiagonalb,querylength,
					 /*pos5*/0,/*pos3*/querylength,/*plusp*/true,/*genestrand*/0,
					 /*query_unk_mismatch_local_p*/false,/*kmer*/index1part_tr);
      pos3b = Genomebits_first_kmer_right(&nmismatches_b3,transcriptomebits,query_compress_fwd,
					  /*univdiagonal*/(Univcoord_T) trdiagonalb,querylength,
					  /*pos5*/0,/*pos3*/querylength,/*plusp*/true,/*genestrand*/0,
					  /*query_unk_mismatch_local_p*/false,/*kmer*/index1part_tr);
	  
      debug4(printf("pos5a, pos3a %d..%d.  pos5b, pos3b %d..%d\n",pos5a,pos3a,pos5b,pos3b));
      if (trdiagonala /*+qstart:0*/ < (Trcoord_T) querylength) {
	/* Skip */

      } else if (pos5a == 0 || pos3b == querylength) {
	/* Setting a first, b second */
	trnum = EF64_trnum(&troffset,&trhigh,transcript_ef64,trdiagonala - querylength,trdiagonala);
	chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
	if ((trpath = Trpath_solve_from_ends(&(*found_score),
					     trdiagonala,pos5a,pos3a,trdiagonalb,pos5b,pos3b,
					     /*tplusp*/true,querylength,/*query_compress_tr*/query_compress_fwd,
					     trnum,troffset,trhigh,chrnum,transcript_genestrand,
					     stage1->indelinfo,max_insertionlen,max_deletionlen,
					     intlistpool,uintlistpool,listpool,
					     trpathpool,pathpool,method)) != NULL) {
	  *sense_trpaths = Hitlist_push(*sense_trpaths,hitlistpool,(void *) trpath
					hitlistpool_trace(__FILE__,__LINE__));
	}

      } else if (pos5b == 0 || pos3a == querylength) {
	/* Setting b first, a second */
	trnum = EF64_trnum(&troffset,&trhigh,transcript_ef64,trdiagonala - querylength,trdiagonala);
	chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
	if ((trpath = Trpath_solve_from_ends(&(*found_score),
					     trdiagonalb,pos5b,pos3b,trdiagonala,pos5a,pos3a,
					     /*tplusp*/true,querylength,/*query_compress_tr*/query_compress_fwd,
					     trnum,troffset,trhigh,chrnum,transcript_genestrand,
					     stage1->indelinfo,max_insertionlen,max_deletionlen,
					     intlistpool,uintlistpool,listpool,
					     trpathpool,pathpool,method)) != NULL) {
	  debug4(printf("Created trpath\n"));
	  *sense_trpaths = Hitlist_push(*sense_trpaths,hitlistpool,(void *) trpath
					hitlistpool_trace(__FILE__,__LINE__));
	}

      } else {
	/* Skip.  Alignment does not go to the ends of the read */
      }
    }

    k = j;
  }

  /* Case for n_tminus_diagpairs == 0 already handled above */
  k = 0;
  while (k < 2*n_tminus_diagpairs) {
    trdiagonala = tminus_diagpairs[k];
    trdiagonalb = tminus_diagpairs[k+1];
    if (trdiagonalb <= trdiagonala) {
      min_distance = trdiagonala - trdiagonalb;
    } else {
      min_distance = trdiagonalb - trdiagonala;
    }

    /* Pick the closest pair */
    j = k + 2;
    while (j < n_tminus_diagpairs && tminus_diagpairs[j] == trdiagonala) {
      if (tminus_diagpairs[j+1] <= trdiagonala) {
	distance = trdiagonala - tminus_diagpairs[j+1];
      } else {
	distance = tminus_diagpairs[j+1] - trdiagonala;
      }
      if (distance < min_distance) {
	min_distance = distance;
	trdiagonalb = tminus_diagpairs[j+1];
      }
      j += 2;
    }

    debug4(printf("tminus diagpairs: %u and %u\n",trdiagonala,trdiagonalb));
    if (trdiagonalb == trdiagonala) {
      /* Skip exact results, which were covered by TR_EXACT */
    } else {
      pos5a = Genomebits_first_kmer_left(&nmismatches_a5,transcriptomebits,query_compress_rev,
					 /*univdiagonal*/(Univcoord_T) trdiagonala,querylength,
					 /*pos5*/0,/*pos3*/querylength,/*plusp*/false,/*genestrand*/0,
					 /*query_unk_mismatch_local_p*/false,/*kmer*/index1part_tr);
      pos3a = Genomebits_first_kmer_right(&nmismatches_a3,transcriptomebits,query_compress_rev,
					  /*univdiagonal*/(Univcoord_T) trdiagonala,querylength,
					  /*pos5*/0,/*pos3*/querylength,/*plusp*/false,/*genestrand*/0,
					  /*query_unk_mismatch_local_p*/false,/*kmer*/index1part_tr);
	  
      pos5b = Genomebits_first_kmer_left(&nmismatches_b5,transcriptomebits,query_compress_rev,
					 /*univdiagonal*/(Univcoord_T) trdiagonalb,querylength,
					 /*pos5*/0,/*pos3*/querylength,/*plusp*/false,/*genestrand*/0,
					 /*query_unk_mismatch_local_p*/false,/*kmer*/index1part_tr);
      pos3b = Genomebits_first_kmer_right(&nmismatches_b3,transcriptomebits,query_compress_rev,
					  /*univdiagonal*/(Univcoord_T) trdiagonalb,querylength,
					  /*pos5*/0,/*pos3*/querylength,/*plusp*/false,/*genestrand*/0,
					  /*query_unk_mismatch_local_p*/false,/*kmer*/index1part_tr);

      debug4(printf("pos5a, pos3a %d..%d.  pos5b, pos3b %d..%d\n",pos5a,pos3a,pos5b,pos3b));
      if (trdiagonala /*+qstart:0*/ < (Trcoord_T) querylength) {
	/* Skip */

      } else if (pos5a == 0 || pos3b == querylength) {
	/* Setting a first, b second */
	trnum = EF64_trnum(&troffset,&trhigh,transcript_ef64,trdiagonala - querylength,trdiagonala);
	chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
	if ((trpath = Trpath_solve_from_ends(&(*found_score),
					     trdiagonala,pos5a,pos3a,trdiagonalb,pos5b,pos3b,
					     /*tplusp*/false,querylength,/*query_compress_tr*/query_compress_rev,
					     trnum,troffset,trhigh,chrnum,transcript_genestrand,
					     stage1->indelinfo,max_insertionlen,max_deletionlen,
					     intlistpool,uintlistpool,listpool,
					     trpathpool,pathpool,method)) != NULL) {
	  debug4(printf("Created trpath\n"));
	  *antisense_trpaths = Hitlist_push(*antisense_trpaths,hitlistpool,(void *) trpath
					    hitlistpool_trace(__FILE__,__LINE__));
	}

      } else if (pos5b == 0 || pos3a == querylength) {
	/* Setting b first, a second */
	trnum = EF64_trnum(&troffset,&trhigh,transcript_ef64,trdiagonala - querylength,trdiagonala);
	chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
	if ((trpath = Trpath_solve_from_ends(&(*found_score),
					     trdiagonalb,pos5b,pos3b,trdiagonala,pos5a,pos3a,
					     /*tplusp*/false,querylength,/*query_compress_tr*/query_compress_rev,
					     trnum,troffset,trhigh,chrnum,transcript_genestrand,
					     stage1->indelinfo,max_insertionlen,max_deletionlen,
					     intlistpool,uintlistpool,listpool,
					     trpathpool,pathpool,method)) != NULL) {
	  debug4(printf("Created trpath\n"));
	  *antisense_trpaths = Hitlist_push(*antisense_trpaths,hitlistpool,(void *) trpath
					    hitlistpool_trace(__FILE__,__LINE__));
	}

      } else {
	/* Skip.  Read does not extend to the ends, so this method is not applicable */
      }
    }

    k = j;
  }

  FREE(tminus_diagpairs);
  FREE(tplus_diagpairs);

  return;
}
#endif


void
Transcriptome_search_setup (int index1part_tr_in, int index1interval_tr_in,
			    Indexdb_T tr_indexdb_in, Transcriptome_T transcriptome_in,
			    EF64_T transcript_ef64_in, Genomebits_T genomebits_in,
			    Genomebits_T genomebits_alt_in, Genomebits_T transcriptomebits_in) {

  index1interval_tr = index1interval_tr_in;
  index1part_tr = index1part_tr_in;
  
  tr_indexdb = tr_indexdb_in;

  transcriptome = transcriptome_in;
  transcript_ef64 = transcript_ef64_in;

  genomebits = genomebits_in;
  genomebits_alt = genomebits_alt_in;
  transcriptomebits = transcriptomebits_in;

#ifdef HAVE_64_BIT
  leftreadshift_tr = 64 - index1part_tr - index1part_tr;
  oligobase_mask_tr = ~(~ (Oligospace_T) 0 << 2*index1part_tr);
#else
  leftreadshift_tr = 32 - index1part_tr - index1part_tr;
  oligobase_mask_tr = ~(~ (Oligospace_T) 0 << 2*index1part_tr);
#endif

  return;
}

