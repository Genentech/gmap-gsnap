static char rcsid[] = "$Id: 46aaf008199b78d8c408f9c1103f1ed1476fe420 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "altsplice.h"
#include "genomicpos.h"
#include "univdiagdef.h"

#include <stdio.h>
#include <string.h>		/* For memcpy */
#include "assert.h"
#include "splice.h"		/* For splice probs */
#include "sedgesort.h"
#include "genomebits_count.h"


static int max_insertlength;

static Genomebits_T genomebits;
static Genomebits_T genomebits_alt;


/* Altsplice creation */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* Altsplice best */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


/* Altsplice_trim */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif


/* Altsplice_select procedures */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

/* Altsplice_resolve procedures */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif



#define T Altsplice_T


void
Altsplice_free (T *old, Pathpool_T pathpool) {
#ifdef DEBUG0
  static int call_i = 0;
#endif

  if (*old) {
    debug0(printf("%d: Altsplice_free of %p\n",++call_i,*old));

    /* Allocated by Vectorpool_T, and not reclaiming due to variable lengths */
    /* FREE_OUT((*old)->coords); */
    /* FREE_OUT((*old)->distal_qpos); */
    /* FREE_OUT((*old)->distal_trimpos); */
    /* FREE_OUT((*old)->nmismatches); */
    /* FREE_OUT((*old)->ref_nmismatches); */
    /* FREE_OUT((*old)->distal_probs); */
    Pathpool_free_altsplice(&(*old),pathpool
			    pathpool_trace(__FILE__,__LINE__)); /* Allocated by Pathpool_new_altsplice */
  }
  return;
}


void
Altsplice_print (T this) {
  int i;
  Chrpos_T splice_distance;
  int distal_length;

  if (this != NULL) {
    printf("%p splice_qpos:%d",this,this->splice_qpos);
    if (this->boundedp == true) {
      printf(" BOUNDED");
    }

    for (i = 0; i < this->nunivdiagonals; i++) {
      printf(" %f..%f",this->medial_prob,this->distal_probs[i]);

      if (this->univdiagonals[i] > this->medial_univdiagonal) {
	splice_distance = this->univdiagonals[i] - this->medial_univdiagonal;
      } else {
	splice_distance = this->medial_univdiagonal - this->univdiagonals[i];
      }

      if (this->distal_trimpos[i] < this->splice_qpos) {
	distal_length = this->splice_qpos - this->distal_trimpos[i];
      } else {
	distal_length = this->distal_trimpos[i] - this->splice_qpos;
      }

#ifdef LARGE_GENOMES
      /* printf(" %u,%llu(%f-%f)->%d",splice_distance,this->univdiagonals[i],this->medial_prob,this->distal_probs[i],
	 this->distal_trimpos[i]); */
      printf(":%u,len:%d,%llu(%d,%d)",
	     splice_distance,this->univdiagonals[i],distal_length,this->medial_nmismatches,this->distal_nmismatches[i]);
#else
      /* printf(" %u,len:%u(%f-%f)->%d",splice_distance,this->univdiagonals[i],this->medial_prob,this->distal_probs[i],
	 this->distal_trimpos[i]); */
      printf(":%u,%d,%u(%d,%d)",splice_distance,distal_length,this->univdiagonals[i],this->medial_nmismatches,this->distal_nmismatches[i]);
#endif
    }
  }
  return;
}


T
Altsplice_copy (T old, Pathpool_T pathpool, Vectorpool_T vectorpool) {
  T new;
  int npartners;

#ifdef DEBUG0
  static int call_i = 0;
#endif

  if (old == NULL) {
    return (T) NULL;
  } else {
    new = (T) Pathpool_new_altsplice(pathpool
				     pathpool_trace(__FILE__,__LINE__));

    debug0(printf("%d: Altsplice_copy of %p => %p\n",++call_i,old,new));

    new->boundedp = old->boundedp;
    new->innerp = old->innerp;
    new->nmismatches_known_p = old->nmismatches_known_p;

    new->splice_qpos = old->splice_qpos;
    new->anchor_qpos = old->anchor_qpos;
    new->medial_univdiagonal = old->medial_univdiagonal;
    new->medial_nmismatches = old->medial_nmismatches;
    new->medial_prob = old->medial_prob;
    new->medial_support = old->medial_support;

    new->best_distali = old->best_distali;

    npartners = new->nunivdiagonals = old->nunivdiagonals;

    new->univdiagonals = Vectorpool_new_univcoordvector(vectorpool,npartners);
    new->distal_lengths = Vectorpool_new_intvector(vectorpool,npartners);
    new->distal_qpos = Vectorpool_new_intvector(vectorpool,npartners);
    new->distal_trimpos = Vectorpool_new_intvector(vectorpool,npartners);
    new->distal_nmismatches = Vectorpool_new_intvector(vectorpool,npartners);
    new->distal_probs = Vectorpool_new_doublevector(vectorpool,npartners);
    new->distal_support = Vectorpool_new_intvector(vectorpool,npartners);

    memcpy(new->univdiagonals,old->univdiagonals,npartners*sizeof(Univcoord_T));
    memcpy(new->distal_lengths,old->distal_lengths,npartners*sizeof(int));
    memcpy(new->distal_qpos,old->distal_qpos,npartners*sizeof(int));
    memcpy(new->distal_trimpos,old->distal_trimpos,npartners*sizeof(int));
    memcpy(new->distal_nmismatches,old->distal_nmismatches,npartners*sizeof(int));
    memcpy(new->distal_probs,old->distal_probs,npartners*sizeof(double));
    memcpy(new->distal_support,old->distal_support,npartners*sizeof(int));

    return new;
  }
}


static bool
check_ascending (Univcoord_T *univdiagonals, int *order_medial_qpos, int n) {
  Univcoord_T prev_univdiagonal;
  int i, ii;

  i = (order_medial_qpos == NULL) ? 0 : order_medial_qpos[0];
  prev_univdiagonal = univdiagonals[i];

  for (ii = 1; ii < n; ii++) {
    i = (order_medial_qpos == NULL) ? ii : order_medial_qpos[ii];
    if (univdiagonals[i] < prev_univdiagonal) {
      return false;
    }
    prev_univdiagonal = univdiagonals[i];
  }
 
  return true;
}

static bool
check_descending (Univcoord_T *univdiagonals, int *order_medial_qpos, int n) {
  Univcoord_T prev_univdiagonal;
  int i, ii;

  i = (order_medial_qpos == NULL) ? 0 : order_medial_qpos[0];
  prev_univdiagonal = univdiagonals[i];

  for (ii = 1; ii < n; ii++) {
    i = (order_medial_qpos == NULL) ? ii : order_medial_qpos[ii];
    if (univdiagonals[i] > prev_univdiagonal) {
      return false;
    }
    prev_univdiagonal = univdiagonals[i];
  }
 
  return true;
}


Univcoord_T
Altsplice_qstart_best (int *distal_trimpos, int splice_qpos, int querylength,
		       Univcoord_T medial_univdiagonal,
		       double medial_prob, int medial_support, int medial_nmismatches,
		       int *distal_support, Univdiag_T *distal_univdiags,
		       double *distal_probs, int npartners,
		       bool plusp, bool sense_forward_p, bool innerp) {
  Univcoord_T distal_univdiagonal = 0, first_dist = 0, second_dist = 0, dist;
  int i;

#ifdef DEBUG1
  printf("Entering Altsplice_qstart_best with medial univdiagonal %u\n",medial_univdiagonal);
  for (i = 0; i < npartners; i++) {
    printf("%u (dist:%u) %d..%d %f\n",
	   distal_univdiags[i]->univdiagonal,medial_univdiagonal - distal_univdiags[i]->univdiagonal,
	   distal_univdiags[i]->qstart,distal_univdiags[i]->qend,distal_probs[i]);
  }
#endif

  for (i = 0; i < npartners; i++) {
    assert(distal_univdiags[i]->univdiagonal < medial_univdiagonal);
    /* Need to use distal_support when a distal indel is found by Spliceends */
    if (Splice_accept_p(splice_qpos,querylength,medial_prob,distal_probs[i],
			/*univdiagonal_i*/distal_univdiags[i]->univdiagonal,
			/*univdiagonal_j*/medial_univdiagonal,
			medial_support,medial_nmismatches,distal_support[i],
			/* distal_univdiags[i]->qend - distal_univdiags[i]->qstart, */
			distal_univdiags[i]->nmismatches,
			plusp,sense_forward_p,innerp) == false) {
      /* Skip */
    } else {
      dist = medial_univdiagonal - distal_univdiags[i]->univdiagonal;
      if (first_dist == 0 || dist < first_dist) {
	second_dist = first_dist;
	first_dist = dist;
	distal_univdiagonal = distal_univdiags[i]->univdiagonal;
	*distal_trimpos = distal_univdiags[i]->qstart;

      } else if (second_dist == 0 || dist < second_dist) {
	second_dist = dist;
      }
    }
  }
	
  if (second_dist > 100 * first_dist) {
    /* first_dist is clearly shorter */
    debug1(printf("Returning %u with distal_trimpos %d\n",distal_univdiagonal,*distal_trimpos));
    return distal_univdiagonal;
  } else {
    debug1(printf("Returning 0\n"));
    return 0;
  }
}


Univcoord_T
Altsplice_qend_best (int *distal_trimpos, int splice_qpos, int querylength,
		     Univcoord_T medial_univdiagonal,
		     double medial_prob, int medial_support, int medial_nmismatches,
		     int *distal_support, Univdiag_T *distal_univdiags,
		     double *distal_probs, int npartners,
		     bool plusp, bool sense_forward_p, bool innerp) {
  Univcoord_T distal_univdiagonal = 0, first_dist = 0, second_dist = 0, dist;
  int i;

#ifdef DEBUG1
  printf("Entering Altsplice_qend_best with medial univdiagonal %u\n",medial_univdiagonal);
  for (i = 0; i < npartners; i++) {
    printf("%u (dist:%u) %d..%d %f\n",
	   distal_univdiags[i]->univdiagonal,distal_univdiags[i]->univdiagonal - medial_univdiagonal,
	   distal_univdiags[i]->qstart,distal_univdiags[i]->qend,distal_probs[i]);
  }
#endif

  for (i = 0; i < npartners; i++) {
    assert(medial_univdiagonal < distal_univdiags[i]->univdiagonal);
    /* Need to use distal_support when a distal indel is found by Spliceends */
    if (Splice_accept_p(splice_qpos,querylength,medial_prob,distal_probs[i],
			/*univdiagonal_i*/medial_univdiagonal,
			/*univdiagonal_j*/distal_univdiags[i]->univdiagonal,
			medial_support,medial_nmismatches,distal_support[i],
			/* distal_univdiags[i]->qend - distal_univdiags[i]->qstart, */
			distal_univdiags[i]->nmismatches,plusp,sense_forward_p,innerp) == false) {
      /* Skip */
    } else {
      dist = distal_univdiags[i]->univdiagonal - medial_univdiagonal;
      if (first_dist == 0 || dist < first_dist) {
	second_dist = first_dist;
	first_dist = dist;
	distal_univdiagonal = distal_univdiags[i]->univdiagonal;
	*distal_trimpos = distal_univdiags[i]->qend;

      } else if (second_dist == 0 || dist < second_dist) {
	second_dist = dist;
      }
    }
  }

  if (second_dist > 100 * first_dist) {
    /* first_dist is clearly shorter */
    debug1(printf("Returning %u with distal_trimpos %d\n",distal_univdiagonal,*distal_trimpos));
    return distal_univdiagonal;
  } else {
    debug1(printf("Returning 0\n"));
    return 0;
  }
}


T
Altsplice_qstart_new (int *endpoint, bool boundedp, int splice_qpos, int anchor_qpos, int querylength,
		      Univcoord_T medial_univdiagonal, int medial_nmismatches,
		      double medial_prob, int medial_support, Univcoord_T *distal_positions,
		      int *distal_support, Univdiag_T *distal_univdiags,
		      double *distal_probs, int *order_medial_qpos, int npartners,
		      bool plusp, bool sense_forward_p, bool innerp,
		      Pathpool_T pathpool, Vectorpool_T vectorpool, bool sort_bydistal_p) {
  T new;
  int *order_distal_positions, i, ii, oldi, best_distali, k;
  double best_distal_prob;


#ifdef DEBUG0
  static int call_i = 0;
#endif
  
  debug0(printf("%d: Creating Altsplice_qstart_new %p with splice_qpos %d, prob %f\n",
		++call_i,new,splice_qpos,medial_prob));
#ifdef DEBUG0
  for (ii = 0; ii < npartners; ii++) {
    i = (order_medial_qpos == NULL) ? ii : order_medial_qpos[ii];
    printf("%u (dist:%u) %d..%d nmismatches:%d %f\n",
	   distal_univdiags[i]->univdiagonal,medial_univdiagonal - distal_univdiags[i]->univdiagonal,
	   distal_univdiags[i]->qstart,distal_univdiags[i]->qend,distal_univdiags[i]->nmismatches,distal_probs[i]);
  }
#endif

  /* Find reasonable values for Path_eval_nmatches based on (1) best
     medial prob and (2) best distal prob for that */

  best_distal_prob = 0.0;
  k = 0;
  for (ii = 0; ii < npartners; ii++) {
    i = (order_medial_qpos == NULL) ? ii : order_medial_qpos[ii];
    assert(distal_univdiags[i]->univdiagonal < medial_univdiagonal);
    /* Need to use distal_support when a distal indel is found by Spliceends */
    if (Splice_accept_p(splice_qpos,querylength,medial_prob,distal_probs[i],
			/*univdiagonal_i*/distal_univdiags[i]->univdiagonal,
			/*univdiagonal_j*/medial_univdiagonal,
			medial_support,medial_nmismatches,distal_support[i],
			/* distal_univdiags[i]->qend - distal_univdiags[i]->qstart, */
			distal_univdiags[i]->nmismatches,plusp,sense_forward_p,innerp) == false) {
      /* Skip */
    } else {
      if (k == 0 || distal_probs[i] > best_distal_prob) {
	best_distali = k;
      }
      k++;
    }
  }

  if (k == 0) {
    return (T) NULL;
  } else {
    new = Pathpool_new_altsplice(pathpool
				 pathpool_trace(__FILE__,__LINE__));
  }

  new->boundedp = boundedp;
  new->innerp = innerp;
  new->nmismatches_known_p = true;

  new->splice_qpos = splice_qpos;
  new->anchor_qpos = anchor_qpos;
  new->medial_univdiagonal = medial_univdiagonal;
  new->medial_nmismatches = medial_nmismatches;
  new->medial_prob = medial_prob;
  new->medial_support = medial_support;

  new->best_distali = best_distali;

  new->nunivdiagonals = k;

  new->univdiagonals = Vectorpool_new_univcoordvector(vectorpool,npartners);
  new->distal_lengths = Vectorpool_new_intvector(vectorpool,npartners);
  new->distal_qpos = Vectorpool_new_intvector(vectorpool,npartners);
  new->distal_trimpos = Vectorpool_new_intvector(vectorpool,npartners);
  new->distal_nmismatches = Vectorpool_new_intvector(vectorpool,npartners);
  new->distal_probs = Vectorpool_new_doublevector(vectorpool,npartners);
  new->distal_support = Vectorpool_new_intvector(vectorpool,npartners);

  if (sort_bydistal_p == false || new->nunivdiagonals == 1 ||
      check_descending(distal_positions,order_medial_qpos,npartners) == true) {

    k = 0;
    for (ii = 0; ii < npartners; ii++) {
      i = (order_medial_qpos == NULL) ? ii : order_medial_qpos[ii];
      /* Need to use distal_support when a distal indel is found by Spliceends */
      if (Splice_accept_p(splice_qpos,querylength,medial_prob,distal_probs[i],
			  /*univdiagonal_i*/distal_univdiags[i]->univdiagonal,
			  /*univdiagonal_j*/medial_univdiagonal,
			  medial_support,medial_nmismatches,distal_support[i],
			  /* distal_univdiags[i]->qend - distal_univdiags[i]->qstart, */
			  distal_univdiags[i]->nmismatches,plusp,sense_forward_p,innerp) == false) {
	/* Skip */
      } else {
	new->univdiagonals[k] = distal_univdiags[i]->univdiagonal;
	new->distal_lengths[k] = distal_support[i];
	new->distal_qpos[k] = distal_univdiags[i]->qend;
	new->distal_trimpos[k] = distal_univdiags[i]->qstart;
	if (distal_univdiags[i]->qend == splice_qpos) {
	  new->distal_nmismatches[k] = distal_univdiags[i]->nmismatches;
#if 0
	} else if (distal_univdiags[i]->qend < splice_qpos &&
		   distal_univdiags[i]->nmismatches == 0) {
	  new->distal_nmismatches[k] = 0;
#endif
	} else {
	  new->distal_nmismatches[k] = -1;
	  new->nmismatches_known_p = false;
	}
	new->distal_probs[k] = distal_probs[i];
	new->distal_support[k] = distal_support[i];
	k++;
      }
    }

  } else {
    assert(order_medial_qpos == NULL);
#ifdef LARGE_GENOMES
    order_distal_positions = Sedgesort_order_uint8(distal_positions,npartners);
#else
    order_distal_positions = Sedgesort_order_uint4(distal_positions,npartners);
#endif
    
    k = 0;
    for (i = npartners - 1; i >= 0; i--) {
      oldi = order_distal_positions[i];

      /* Need to use distal_support when a distal indel is found by Spliceends */
      if (Splice_accept_p(splice_qpos,querylength,medial_prob,distal_probs[oldi],
			  /*univdiagonal_i*/distal_univdiags[oldi]->univdiagonal,
			  /*univdiagonal_j*/medial_univdiagonal,
			  medial_support,medial_nmismatches,distal_support[oldi],
			  /* distal_univdiags[oldi]->qend - distal_univdiags[oldi]->qstart, */
			  distal_univdiags[oldi]->nmismatches,plusp,sense_forward_p,innerp) == false) {
	/* Skip */
      } else {
	new->univdiagonals[k] = distal_univdiags[oldi]->univdiagonal;
	new->distal_lengths[k] = distal_univdiags[oldi]->qend - distal_univdiags[oldi]->qstart;
	new->distal_qpos[k] = distal_univdiags[oldi]->qend;
	new->distal_trimpos[k] = distal_univdiags[oldi]->qstart;
	if (distal_univdiags[oldi]->qend == splice_qpos) {
	  new->distal_nmismatches[k] = distal_univdiags[oldi]->nmismatches;
#if 0
	} else if (distal_univdiags[oldi]->qend < splice_qpos &&
		   distal_univdiags[oldi]->nmismatches == 0) {
	  new->distal_nmismatches[k] = 0;
#endif
	} else {
	  new->distal_nmismatches[k] = -1;
	  new->nmismatches_known_p = false;
	}
	new->distal_probs[k] = distal_probs[oldi];
	new->distal_support[k] = distal_support[oldi];
	k++;
      }
    }

    FREE(order_distal_positions);
  }

  *endpoint = splice_qpos;

#ifdef DEBUG0
  printf("best distali is %d\n",new->best_distali);
  printf("best univdiagonal is %u\n",Altsplice_best_univdiagonal(new));
  printf("best distal prob is %f\n",Altsplice_best_distal_prob(new));
  printf("best nmatches is %d\n",Altsplice_best_nmatches(new));
  printf("found score is %d\n",Altsplice_found_score(new));
#endif

  return new;
}


T
Altsplice_qend_new (int *endpoint, bool boundedp, int splice_qpos, int anchor_qpos, int querylength,
		    Univcoord_T medial_univdiagonal, int medial_nmismatches,
		    double medial_prob, int medial_support, Univcoord_T *distal_positions,
		    int *distal_support, Univdiag_T *distal_univdiags,
		    double *distal_probs, int *order_medial_qpos, int npartners,
		    bool plusp, bool sense_forward_p, bool innerp,
		    Pathpool_T pathpool, Vectorpool_T vectorpool, bool sort_bydistal_p) {
  T new;
  int *order_distal_positions, i, ii, oldi, best_distali, k;
  double best_distal_prob;


#ifdef DEBUG0
  static int call_i = 0;
#endif
  
  debug0(printf("%d: Creating Altsplice_qend_new %p with splice_qpos %d, prob %f\n",
		++call_i,new,splice_qpos,medial_prob));
#ifdef DEBUG0
  for (ii = 0; ii < npartners; ii++) {
    i = (order_medial_qpos == NULL) ? ii : order_medial_qpos[ii];
    printf("%u (dist:%u) %d..%d nmismatches:%d %f\n",
	   distal_univdiags[i]->univdiagonal,distal_univdiags[i]->univdiagonal - medial_univdiagonal,
	   distal_univdiags[i]->qstart,distal_univdiags[i]->qend,distal_univdiags[i]->nmismatches,distal_probs[i]);
  }
#endif


  /* Find reasonable values for Path_eval_nmatches based on (1) best
     medial prob and (2) best distal prob for that */

  best_distal_prob = 0.0;
  k = 0;
  for (ii = 0; ii < npartners; ii++) {
    i = (order_medial_qpos == NULL) ? ii : order_medial_qpos[ii];
    assert(medial_univdiagonal < distal_univdiags[i]->univdiagonal);
    /* Need to use distal_support when a distal indel is found by Spliceends */
    if (Splice_accept_p(splice_qpos,querylength,medial_prob,distal_probs[i],
			/*univdiagonal_i*/medial_univdiagonal,
			/*univdiagonal_j*/distal_univdiags[i]->univdiagonal,
			medial_support,medial_nmismatches,distal_support[i],
			/* distal_univdiags[i]->qend - distal_univdiags[i]->qstart, */
			distal_univdiags[i]->nmismatches,plusp,sense_forward_p,innerp) == false) {
      /* Skip */
    } else {
      if (k == 0 || distal_probs[i] > best_distal_prob) {
	best_distali = k;
      }
      k++;
    }
  }

  if (k == 0) {
    return (T) NULL;

  } else {
    new = Pathpool_new_altsplice(pathpool
				 pathpool_trace(__FILE__,__LINE__));
  }

  new->boundedp = boundedp;
  new->innerp = innerp;
  new->nmismatches_known_p = true;

  new->splice_qpos = splice_qpos;
  new->anchor_qpos = anchor_qpos;
  new->medial_univdiagonal = medial_univdiagonal;
  new->medial_nmismatches = medial_nmismatches;
  new->medial_prob = medial_prob;
  new->medial_support = medial_support;
  new->best_distali = best_distali;

  new->nunivdiagonals = k;

  new->univdiagonals = Vectorpool_new_univcoordvector(vectorpool,npartners);
  new->distal_lengths = Vectorpool_new_intvector(vectorpool,npartners);
  new->distal_qpos = Vectorpool_new_intvector(vectorpool,npartners);
  new->distal_trimpos = Vectorpool_new_intvector(vectorpool,npartners);
  new->distal_nmismatches = Vectorpool_new_intvector(vectorpool,npartners);
  new->distal_probs = Vectorpool_new_doublevector(vectorpool,npartners);
  new->distal_support = Vectorpool_new_intvector(vectorpool,npartners);

  if (sort_bydistal_p == false || new->nunivdiagonals == 1 ||
      check_ascending(distal_positions,order_medial_qpos,npartners) == true) {

    k = 0;
    for (ii = 0; ii < npartners; ii++) {
      i = (order_medial_qpos == NULL) ? ii : order_medial_qpos[ii];
      /* Need to use distal_support when a distal indel is found by Spliceends */
      if (Splice_accept_p(splice_qpos,querylength,medial_prob,distal_probs[i],
			  /*univdiagonal_i*/medial_univdiagonal,
			  /*univdiagonal_j*/distal_univdiags[i]->univdiagonal,
			  medial_support,medial_nmismatches,distal_support[i],
			  /* distal_univdiags[i]->qend - distal_univdiags[i]->qstart,*/
			  distal_univdiags[i]->nmismatches,plusp,sense_forward_p,innerp) == false) {
	/* Skip */
      } else {
	new->univdiagonals[k] = distal_univdiags[i]->univdiagonal;
	new->distal_lengths[k] = distal_univdiags[i]->qend - distal_univdiags[i]->qstart;
	new->distal_qpos[k] = distal_univdiags[i]->qstart;
	new->distal_trimpos[k] = distal_univdiags[i]->qend;
	if (distal_univdiags[i]->qstart == splice_qpos) {
	  new->distal_nmismatches[k] = distal_univdiags[i]->nmismatches;
#if 0
	} else if (distal_univdiags[i]->qstart > splice_qpos &&
		   distal_univdiags[i]->nmismatches == 0) {
	  new->distal_nmismatches[k] = 0;
#endif
	} else {
	  new->distal_nmismatches[k] = -1;
	  new->nmismatches_known_p = false;
	}
	new->distal_probs[k] = distal_probs[i];
	new->distal_support[k] = distal_support[i];
	k++;
      }
    }

  } else {
    assert(order_medial_qpos == NULL);
#ifdef LARGE_GENOMES
    order_distal_positions = Sedgesort_order_uint8(distal_positions,npartners);
#else
    order_distal_positions = Sedgesort_order_uint4(distal_positions,npartners);
#endif
    
    k = 0;
    for (i = 0; i < npartners; i++) {
      oldi = order_distal_positions[i];

      /* Need to use distal_support when a distal indel is found by Spliceends */
      if (Splice_accept_p(splice_qpos,querylength,medial_prob,distal_probs[oldi],
			  /*univdiagonal_i*/medial_univdiagonal,
			  /*univdiagonal_j*/distal_univdiags[oldi]->univdiagonal,
			  medial_support,medial_nmismatches,distal_support[oldi],
			  /* distal_univdiags[oldi]->qend - distal_univdiags[oldi]->qstart, */
			  distal_univdiags[oldi]->nmismatches,plusp,sense_forward_p,innerp) == false) {
	/* Skip */
      } else {
	new->univdiagonals[k] = distal_univdiags[oldi]->univdiagonal;
	new->distal_lengths[k] = distal_univdiags[oldi]->qend - distal_univdiags[oldi]->qstart;
	new->distal_qpos[k] = distal_univdiags[oldi]->qstart;
	new->distal_trimpos[k] = distal_univdiags[oldi]->qend;
	if (distal_univdiags[oldi]->qstart == splice_qpos) {
	  new->distal_nmismatches[k] = distal_univdiags[oldi]->nmismatches;
#if 0
	} else if (distal_univdiags[oldi]->qstart > splice_qpos &&
		   distal_univdiags[oldi]->nmismatches == 0) {
	  new->distal_nmismatches[k] = 0;
#endif
	} else {
	  new->distal_nmismatches[k] = -1;
	  new->nmismatches_known_p = false;
	}
	new->distal_probs[k] = distal_probs[oldi];
	new->distal_support[k] = distal_support[oldi];
	k++;
      }
    }

    FREE(order_distal_positions);
  }

  *endpoint = splice_qpos;

#ifdef DEBUG0
  printf("best distali is %d\n",new->best_distali);
  printf("best univdiagonal is %u\n",Altsplice_best_univdiagonal(new));
  printf("best distal prob is %f\n",Altsplice_best_distal_prob(new));
  printf("best nmatches is %d\n",Altsplice_best_nmatches(new));
  printf("found score is %d\n",Altsplice_found_score(new));
#endif

  return new;
}


Univcoord_T
Altsplice_first_univdiagonal (T this) {
  return this->univdiagonals[0];
}

Univcoord_T
Altsplice_last_univdiagonal (T this) {
  return this->univdiagonals[this->nunivdiagonals - 1];
}


/* Returns whether anything is left */
bool
Altsplice_trim_qstart_chrbounds (T this, Univcoord_T chroffset, int querylength) {
  int endi;
  
  /* Traverse from genome start and go upstream */
  endi = this->nunivdiagonals - 1;
  debug3(printf("Entered Altsplice_trim_qstart_chrbounds from %d to 0\n",endi));
  while (endi >= 0 && this->univdiagonals[endi] < chroffset + querylength) {
    debug3(printf("Univdiagonal %u < chroffset %u + querylength %d, so trimming\n",
		  this->univdiagonals[endi],chroffset,querylength));
    endi--;
  }

  if (endi < 0) {
    return false;
  } else {
    this->nunivdiagonals = endi + 1;
    debug3(printf("Resetting nunivdiagonals to be %d\n",this->nunivdiagonals));
    return true;
  }
}


/* Returns whether anything is left */
bool
Altsplice_trim_qend_chrbounds (T this, Univcoord_T chrhigh) {
  int endj;
  
  /* Traverse from genome end and go downstream */
  endj = this->nunivdiagonals - 1;
  debug3(printf("Entered Altsplice_trim_qend_chrbounds from %d to 0\n",endj));
  while (endj >= 0 && this->univdiagonals[endj] >= chrhigh) {
    debug3(printf("Univdiagonal %u >= chrhigh %u, so trimming\n",
		  this->univdiagonals[endj],chrhigh));
    endj--;
  }

  if (endj < 0) {
    return false;
  } else {
    this->nunivdiagonals = endj + 1;
    debug3(printf("Resetting nunivdiagonals to be %d\n",this->nunivdiagonals));
    return true;
  }
}


/* Should look for insertlengths within max_insertlength (continuous
   exon), and choose among those with the best nmatches and first
   one with a sufficiently high prob.  Could look for sufficiently
   high prob.  Otherwise (intervening exon), insertlength doesn't
   matter, and we should look for the best nmatches and prob */

/* For inner altsplices */
/* Diagonals are in order from closest to farthest from pathL, so need to check array from n-1 to 0 */
bool
Altsplice_resolve_qend (Univcoord_T *univdiagonal_L, int *splice_qpos_L, int *distal_trimpos_L,
			int *medial_nmismatches_L, int *distal_nmismatches_L, 
			double *medial_prob_L, double *distal_prob_L,
			T this, int anchor_qpos_L, int querylengthL, int querylengthH,
			Univcoord_T genomicstartH, Compress_T query_compress,
			bool plusp, int genestrand, bool sense_forward_p) {

  /* Candidates within max_insertlength are in starti..endi */
  int besti = -1, starti, endi, i;
  int ncandidates;
  int insertlength;
  Univcoord_T genomicendL;
  int best_nmatches_continuous, best_nmatches_intervening, nmatches, ref_nmismatches;
  double best_prob, prob;


#ifdef DEBUG9
  printf("Altsplice_resolve_qend for univdiagonal %u\n",this->medial_univdiagonal);
  Altsplice_print(this);
  printf("\n");
#endif

#ifdef CHECK_ASSERTIONS
  check_ascending(this->univdiagonals,/*order_medial_qpos*/NULL,this->nunivdiagonals);
#endif

  best_nmatches_continuous = 0;
  best_nmatches_intervening = 0;

  endi = this->nunivdiagonals - 1;
  /* Test for insertlength <= 0, or genomicstartH + querylengthL + querylengthH <= genomicendL */
  while (endi >= 0 && genomicstartH + querylengthL + querylengthH <= this->univdiagonals[endi]) {
    debug9(genomicendL = this->univdiagonals[endi]);
    debug9(insertlength = genomicstartH - genomicendL + querylengthL + querylengthH);
    debug9(printf(" Overreach: (insertlength %d, nmismatches %d+%d, prob %f+%f)\n",
		  insertlength,this->medial_nmismatches,this->distal_nmismatches[endi],
		  this->medial_prob,this->distal_probs[endi]));
    endi--;
  }


  i = starti = endi;
  while (i >= 0) {
    genomicendL = this->univdiagonals[i];
    insertlength = genomicstartH - genomicendL + querylengthL + querylengthH;
    if (insertlength <= max_insertlength) {
      debug9(printf(" Continuous: (insertlength %d, nmatches (%d-%d)-%d, prob %f+%f)\n",
		    insertlength,this->distal_trimpos[i],this->splice_qpos,this->distal_nmismatches[i],
		    this->medial_prob,this->distal_probs[i]));
      if (this->distal_nmismatches[i] < 0) {
	this->distal_nmismatches[i] =
	  Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress,
						/*univdiagonal*/this->univdiagonals[i],querylengthL,
						/*pos5*/this->splice_qpos,/*pos3*/this->distal_trimpos[i],
						plusp,genestrand);
      }
      if ((nmatches = (this->distal_trimpos[i] - this->splice_qpos) - this->distal_nmismatches[i]) > best_nmatches_continuous) {
	best_nmatches_continuous = nmatches;
      }
      starti = i;
    } else {
      debug9(printf(" Intervening: (insertlength %d, nmatches (%d-%d)-%d, prob %f+%f)\n",
		    insertlength,this->distal_trimpos[i],this->splice_qpos,this->distal_nmismatches[i],
		    this->medial_prob,this->distal_probs[i]));
      if (this->distal_nmismatches[i] < 0) {
	this->distal_nmismatches[i] =
	  Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress,
						/*univdiagonal*/this->univdiagonals[i],querylengthL,
						/*pos5*/this->splice_qpos,/*pos3*/this->distal_trimpos[i],
						plusp,genestrand);
      }
      if ((nmatches = (this->distal_trimpos[i] - this->splice_qpos) - this->distal_nmismatches[i]) > best_nmatches_intervening) {
	best_nmatches_intervening = nmatches;
      }
    }
    i--;
  }
  
  debug9(printf("best_nmatches_continuous %d\n",best_nmatches_continuous));
  debug9(printf("best_nmatches_intervening %d\n",best_nmatches_intervening));

  if (best_nmatches_continuous > 0) {
    /* Find shortest (first) good continuous solution */
    ncandidates = 0;
    for (i = endi; i >= starti; i--) {
      assert(this->medial_univdiagonal < this->univdiagonals[i]);
      if ((this->distal_trimpos[i] - this->splice_qpos) - this->distal_nmismatches[i] == best_nmatches_continuous &&
	  Splice_accept_p(this->splice_qpos,querylengthL,this->medial_prob,this->distal_probs[i],
			  /*univdiagonal_i*/this->medial_univdiagonal,
			  /*univdiagonal_j*/this->univdiagonals[i],
			  this->medial_support,this->medial_nmismatches,
			  this->distal_support[i],this->distal_nmismatches[i],
			  plusp,sense_forward_p,this->innerp) == true) {
	if (ncandidates++ == 0) {
	  besti = i;
	}
      }
    }
  }
  debug9(printf("besti after continuous: %d\n",besti));

  if (besti < 0 || best_nmatches_intervening > best_nmatches_continuous) {
    /* Find best intervening solution */
    best_prob = GOOD_DISTAL_PROB;
    for (i = endi - 1; i >= 0; i--) {
      assert(this->medial_univdiagonal < this->univdiagonals[i]);
      if ((this->distal_trimpos[i] - this->splice_qpos) - this->distal_nmismatches[i] == best_nmatches_intervening &&
	  Splice_accept_p(this->splice_qpos,querylengthL,this->medial_prob,this->distal_probs[i],
			  /*univdiagonal_i*/this->medial_univdiagonal,
			  /*univdiagonal_j*/this->univdiagonals[i],
			  this->medial_support,this->medial_nmismatches,
			  this->distal_support[i],this->distal_nmismatches[i],
			  plusp,sense_forward_p,this->innerp) == true &&
	  (prob = this->distal_probs[i]) > best_prob) {
	besti = i;
	best_prob = prob;
      }
    }
  }

  this->nmismatches_known_p = true;

  if (besti < 0) {
    debug9(printf("Returning false\n"));
    return false;

  } else {
    debug9(printf("besti %d.  univdiagonal %u\n",besti,this->univdiagonals[besti]));

    /* Single intervening exon or same exon */
    *univdiagonal_L = this->univdiagonals[besti];
    *splice_qpos_L = this->splice_qpos;
    *distal_trimpos_L = this->distal_trimpos[besti];
    if (anchor_qpos_L != this->anchor_qpos) {
      /* Anchor qpos was changed, probably by combine_leftright_paths */
      *medial_nmismatches_L = -1;
    } else {
      *medial_nmismatches_L = this->medial_nmismatches;
    }
    if (this->distal_qpos[besti] != this->splice_qpos) {
      *distal_nmismatches_L = -1;
    } else {
      *distal_nmismatches_L = this->distal_nmismatches[besti];
    }
    *medial_prob_L = this->medial_prob;
    *distal_prob_L = this->distal_probs[besti];
    return true;
  }
}


/* For inner altsplices */
/* Diagonals are in order from closest to farthest from path, so need to check array from n-1 to 0 */
bool
Altsplice_resolve_qstart (Univcoord_T *univdiagonal_H, int *splice_qpos_H, int *distal_trimpos_H,
			  int *medial_nmismatches_H, int *distal_nmismatches_H,
			  double *medial_prob_H, double *distal_prob_H,
			  T this, int anchor_qpos_H, int querylengthL, int querylengthH,
			  Univcoord_T genomicendL, Compress_T query_compress,
			  bool plusp, int genestrand, bool sense_forward_p) {

  /* Candidates within max_insertlength are in startj..endj */
  int bestj = -1, startj, endj, j;
  int ncandidates;
  int insertlength;
  Univcoord_T genomicstartH;
  int best_nmatches_continuous, best_nmatches_intervening, nmatches, ref_nmismatches;
  double best_prob, prob;

#ifdef DEBUG9
  printf("Altsplice_resolve_qstart for univdiagonal %u\n",this->medial_univdiagonal);
  Altsplice_print(this);
  printf("\n");
#endif

#ifdef CHECK_ASSERTIONS
  check_descending(this->univdiagonals,/*order_medial_qpos*/NULL,this->nunivdiagonals);
#endif

  best_nmatches_continuous = 0;
  best_nmatches_intervening = 0;

  endj = this->nunivdiagonals - 1;
  /* Test for insertlength <= 0, or genomicstartH + querylengthL + querylengthH <= genomicendL */
  while (endj >= 0 && this->univdiagonals[endj] + querylengthL <= genomicendL) {
    debug9(genomicstartH = this->univdiagonals[endj] - querylengthH);
    debug9(insertlength = genomicstartH - genomicendL + querylengthL + querylengthH);
    debug9(printf(" Overreach: (insertlength %d, nmismatches %d+%d, prob %f+%f)\n",
		  insertlength,this->medial_nmismatches,this->distal_nmismatches[endj],
		  this->medial_prob,this->distal_probs[endj]));
    endj--;
  }

  
  j = startj = endj;
  while (j >= 0) {
    genomicstartH = this->univdiagonals[j] - querylengthH;
    insertlength = genomicstartH - genomicendL + querylengthL + querylengthH;
    if (insertlength <= max_insertlength) {
      debug9(printf(" Continuous: (insertlength %d, nmatches (%d-%d)-%d, prob %f+%f)\n",
		    insertlength,this->splice_qpos,this->distal_trimpos[j],this->distal_nmismatches[j],
		    this->medial_prob,this->distal_probs[j]));
      if (this->distal_nmismatches[j] < 0) {
	this->distal_nmismatches[j] =
	  Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress,
						/*univdiagonal*/this->univdiagonals[j],querylengthH,
						/*pos5*/this->distal_trimpos[j],/*pos3*/this->splice_qpos,
						plusp,genestrand);
      }
      if ((nmatches = (this->splice_qpos - this->distal_trimpos[j]) - this->distal_nmismatches[j]) > best_nmatches_continuous) {
	best_nmatches_continuous = nmatches;
      }
      startj = j;
    } else {
      debug9(printf(" Intervening: (insertlength %d, nmatches (%d-%d)-%d, prob %f+%f)\n",
		    insertlength,this->splice_qpos,this->distal_trimpos[j],this->distal_nmismatches[j],
		    this->medial_prob,this->distal_probs[j]));
      if (this->distal_nmismatches[j] < 0) {
	this->distal_nmismatches[j] =
	  Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress,
						/*univdiagonal*/this->univdiagonals[j],querylengthH,
						/*pos5*/this->distal_trimpos[j],/*pos3*/this->splice_qpos,
						plusp,genestrand);
      }
      if ((nmatches = (this->splice_qpos - this->distal_trimpos[j]) - this->distal_nmismatches[j]) > best_nmatches_intervening) {
	best_nmatches_intervening = nmatches;
      }
    }
    j--;
  }

  debug9(printf("best_nmatches_continuous %d\n",best_nmatches_continuous));
  debug9(printf("best_nmatches_intervening %d\n",best_nmatches_intervening));

  if (best_nmatches_continuous > 0) {
    /* Find shortest (first) good continuous solution */
    ncandidates = 0;
    for (j = endj; j >= startj; j--) {
      assert(this->univdiagonals[j] < this->medial_univdiagonal);
      if ((this->splice_qpos - this->distal_trimpos[j]) - this->distal_nmismatches[j] == best_nmatches_continuous &&
	  Splice_accept_p(this->splice_qpos,querylengthH,this->medial_prob,this->distal_probs[j],
			  /*univdiagonal_i*/this->univdiagonals[j],
			  /*univdiagonal_j*/this->medial_univdiagonal,
			  this->medial_support,this->medial_nmismatches,
			  this->distal_support[j],this->distal_nmismatches[j],
			  plusp,sense_forward_p,this->innerp) == true) {
	if (ncandidates++ == 0) {
	  bestj = j;
	}
      }
    }
  }
  debug9(printf("bestj after continuous: %d\n",bestj));

  if (bestj < 0 || best_nmatches_intervening > best_nmatches_continuous) {
    /* Find best intervening solution */
    best_prob = GOOD_DISTAL_PROB;
    for (j = endj - 1; j >= 0; j--) {
      if ((this->splice_qpos - this->distal_trimpos[j]) - this->distal_nmismatches[j] == best_nmatches_intervening &&
	  Splice_accept_p(this->splice_qpos,querylengthH,this->medial_prob,this->distal_probs[j],
			  /*univdiagonal_i*/this->univdiagonals[j],
			  /*univdiagonal_j*/this->medial_univdiagonal,
			  this->medial_support,this->medial_nmismatches,
			  this->distal_support[j],this->distal_nmismatches[j],
			  plusp,sense_forward_p,this->innerp) == true &&
	  (prob = this->distal_probs[j]) > best_prob) {
	bestj = j;
	best_prob = prob;
      }
    }
  }
	  
  this->nmismatches_known_p = true;

  if (bestj < 0) {
    debug9(printf("Returning false\n"));
    return false;

  } else {
    debug9(printf("bestj %d.  univdiagonal %u\n",bestj,this->univdiagonals[bestj]));
    /* Single intervening exon or same exon */
    *univdiagonal_H = this->univdiagonals[bestj];
    *splice_qpos_H = this->splice_qpos;
    *distal_trimpos_H = this->distal_trimpos[bestj];
    if (anchor_qpos_H != this->anchor_qpos) {
      /* Anchor qpos was changed, probably by combine_leftright_paths */
      *medial_nmismatches_H = -1;
    } else {
      *medial_nmismatches_H = this->medial_nmismatches;
    }
    if (this->distal_qpos[bestj] != this->splice_qpos) {
      *distal_nmismatches_H = -1;
    } else {
      *distal_nmismatches_H = this->distal_nmismatches[bestj];
    }
    *medial_prob_H = this->medial_prob;
    *distal_prob_H = this->distal_probs[bestj];
    return true;
  }
}


/* For inner altsplices */
bool
Altsplice_resolve_both (Univcoord_T *univdiagonal_L, int *splice_qpos_L, int *distal_trimpos_L,
			int *medial_nmismatches_L, int *distal_nmismatches_L,
			double *medial_prob_L, double *distal_prob_L,

			Univcoord_T *univdiagonal_H, int *splice_qpos_H, int *distal_trimpos_H,
			int *medial_nmismatches_H, int *distal_nmismatches_H,
			double *medial_prob_H, double *distal_prob_H,

			T thisL, int anchor_qpos_L, T thisH, int anchor_qpos_H,
			int querylengthL, int querylengthH,
			Compress_T query_compress_L, Compress_T query_compress_H,
			bool plusp, int genestrand, bool sense_forward_p) {
  int besti = -1, bestj = -1, besti_intervening = -1, bestj_intervening = -1, i, j;
  int ncandidates;
  int best_insertlength, insertlength;
  Univcoord_T genomicendL, genomicstartH;
  int best_nmatches_continuous, best_nmatches_intervening, nmatches, ref_nmismatches;
  double best_prob, prob;


#ifdef DEBUG9
  printf("Altsplice_resolve_both for univdiagonal_L %u and univdiagonal_H %u\n",
	 thisL->medial_univdiagonal,thisH->medial_univdiagonal);
  Altsplice_print(thisL);
  printf("\n");
  Altsplice_print(thisH);
  printf("\n");
#endif

  best_nmatches_continuous = best_nmatches_intervening = 0;

  for (i = 0; i < thisL->nunivdiagonals; i++) {
    genomicendL = thisL->univdiagonals[i];

    for (j = 0; j < thisH->nunivdiagonals; j++) {
      /* Test for insertlength > 0, or genomicstartH + querylengthL + querylengthH > genomicendL */
      if (thisH->univdiagonals[j] + querylengthL > genomicendL) {
	/* Valid insertlength */

	if (thisL->distal_nmismatches[i] < 0) {
	  thisL->distal_nmismatches[i] =
	    Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_L,
						  /*univdiagonal*/thisL->univdiagonals[i],querylengthL,
						  /*pos5*/thisL->splice_qpos,/*pos3*/thisL->distal_trimpos[i],
						  plusp,genestrand);
	}
	
	if (thisH->distal_nmismatches[j] < 0) {
	  thisH->distal_nmismatches[j] =
	    Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_H,
						  /*univdiagonal*/thisH->univdiagonals[j],querylengthH,
						  /*pos5*/thisH->distal_trimpos[j],/*pos3*/thisH->splice_qpos,
						  plusp,genestrand);
	}

	genomicstartH = thisH->univdiagonals[j] - querylengthH;
	insertlength = genomicstartH - genomicendL + querylengthL + querylengthH;
	debug9(printf(" (insertlength %u, nmatches (%d-%d)-%d + (%d-%d)-%d, prob %f+%f + %f+%f)\n",
		      insertlength,thisL->distal_trimpos[i],thisL->splice_qpos,thisL->distal_nmismatches[i],
		      thisH->splice_qpos,thisH->distal_trimpos[j],thisH->distal_nmismatches[j],
		      thisL->medial_prob,thisL->distal_probs[i],
		      thisH->medial_prob,thisH->distal_probs[j]));
	if (insertlength <= max_insertlength) {
	  if ((nmatches = (thisL->distal_trimpos[i] - thisL->splice_qpos) - thisL->distal_nmismatches[i] +
	       (thisH->splice_qpos - thisH->distal_trimpos[j]) - thisH->distal_nmismatches[j]) > best_nmatches_continuous) {
	    best_nmatches_continuous = nmatches;
	  }
	} else {
	  if ((nmatches = (thisL->distal_trimpos[i] - thisL->splice_qpos) - thisL->distal_nmismatches[i] +
	       (thisH->splice_qpos - thisH->distal_trimpos[j]) - thisH->distal_nmismatches[j]) > best_nmatches_intervening) {
	    best_nmatches_intervening = nmatches;
	  }
	}
      }
    }
  }

  debug9(printf("best_nmatches_continuous %d\n",best_nmatches_continuous));
  debug9(printf("best_nmatches_intervening %d\n",best_nmatches_intervening));

  if (best_nmatches_continuous > 0 || best_nmatches_intervening > 0) {
    debug9(printf("Finding shortest good continuous solution and best intervening solution\n"));

    /* Find shortest good continuous solution and best intervening solution */

    ncandidates = 0;
    best_insertlength = 0;	       /* Will get written with first candidate */
    best_prob = GOOD_DISTAL_PROB + GOOD_DISTAL_PROB;	/* threshold */

    for (i = 0; i < thisL->nunivdiagonals; i++) {
      genomicendL = thisL->univdiagonals[i];
      for (j = 0; j < thisH->nunivdiagonals; j++) {
	/* Test for insertlength > 0, or genomicstartH + querylengthL + querylengthH > genomicendL */
	if (thisH->univdiagonals[j] + querylengthL > genomicendL) {
	  /* Valid insertlength */
	  genomicstartH = thisH->univdiagonals[j] - querylengthH;
	  insertlength = genomicstartH - genomicendL + querylengthL + querylengthH;
	  if (insertlength <= max_insertlength) {
	    assert(thisL->medial_univdiagonal < thisL->univdiagonals[i]);
	    assert(thisH->univdiagonals[j] < thisH->medial_univdiagonal);
	    if ((thisL->distal_trimpos[i] - thisL->splice_qpos) - thisL->distal_nmismatches[i] +
		(thisH->splice_qpos - thisH->distal_trimpos[j]) - thisH->distal_nmismatches[j] == best_nmatches_continuous &&
		Splice_accept_p(thisL->splice_qpos,querylengthL,thisL->medial_prob,thisL->distal_probs[i],
				/*univdiagonal_i*/thisL->medial_univdiagonal,
				/*univdiagonal_j*/thisL->univdiagonals[i],
				thisL->medial_support,thisL->medial_nmismatches,
				thisL->distal_support[i],thisL->distal_nmismatches[i],
				plusp,sense_forward_p,thisL->innerp) == true &&
		Splice_accept_p(thisH->splice_qpos,querylengthH,thisH->medial_prob,thisH->distal_probs[j],
				/*univdiagonal_i*/thisH->univdiagonals[j],
				/*univdiagonal_j*/thisH->medial_univdiagonal,
				thisH->medial_support,thisH->medial_nmismatches,
				thisH->distal_support[j],thisH->distal_nmismatches[j],
				plusp,sense_forward_p,thisH->innerp) == true) {
	      if (ncandidates++ == 0 || insertlength < best_insertlength) {
		besti = i;
		bestj = j;
		best_insertlength = insertlength;
	      }
	    }

	  } else {
	    if ((thisL->distal_trimpos[i] - thisL->splice_qpos) - thisL->distal_nmismatches[i] +
		(thisH->splice_qpos - thisH->distal_trimpos[j]) - thisH->distal_nmismatches[j] == best_nmatches_intervening &&
		(prob = thisL->distal_probs[i] + thisH->distal_probs[j]) > best_prob) {
	      besti_intervening = i;
	      bestj_intervening = j;
	      best_prob = prob;
	    }
	  }
	}
      }
    }
    debug9(printf("best_insertlength %d\n",best_insertlength));
  }


  if (besti < 0 || best_nmatches_intervening > best_nmatches_continuous) {
    besti = besti_intervening;
    bestj = bestj_intervening;
  }

  thisL->nmismatches_known_p = true;
  thisH->nmismatches_known_p = true;

  if (besti < 0) {
    debug9(printf("Returning false\n"));
    return false;

  } else {
    debug9(printf("besti and bestj %d,%d.  univdiagonals %u and %u\n",
		  besti,bestj,thisL->univdiagonals[besti],thisH->univdiagonals[bestj]));

    *univdiagonal_L = thisL->univdiagonals[besti];
    *splice_qpos_L = thisL->splice_qpos;
    *distal_trimpos_L = thisL->distal_trimpos[besti];
    if (anchor_qpos_L != thisL->anchor_qpos) {
      /* Anchor qpos was changed, probably by combine_leftright_paths */
      *medial_nmismatches_L = -1;
    } else {
      *medial_nmismatches_L = thisL->medial_nmismatches;
    }
    if (thisL->distal_qpos[besti] != thisL->splice_qpos) {
      *distal_nmismatches_L = -1;
    } else {
      *distal_nmismatches_L = thisL->distal_nmismatches[besti];
    }
    *medial_prob_L = thisL->medial_prob;
    *distal_prob_L = thisL->distal_probs[besti];

    *univdiagonal_H = thisH->univdiagonals[bestj];
    *splice_qpos_H = thisH->splice_qpos;
    *distal_trimpos_H = thisH->distal_trimpos[bestj];
    if (anchor_qpos_H != thisH->anchor_qpos) {
      /* Anchor qpos was changed, probably by combine_leftright_paths */
      *medial_nmismatches_H = -1;
    } else {
      *medial_nmismatches_H = thisH->medial_nmismatches;
    }
    if (thisH->distal_qpos[bestj] != thisH->splice_qpos) {
      *distal_nmismatches_H = -1;
    } else {
      *distal_nmismatches_H = thisH->distal_nmismatches[bestj];
    }
    *medial_prob_H = thisH->medial_prob;
    *distal_prob_H = thisH->distal_probs[bestj];

    return true;
  }
}
  

/* For outer altsplices.  Returns a single best candidate, if there is one */
/* Diagonals are in order from closest to farthest from path, so need to check array from 0 to n-1 */
int
Altsplice_select_qend (T this, int querylength, Compress_T query_compress, bool plusp, int genestrand) {
  int besti = -1, i;
  int ncandidates;
  int best_nmatches, nmatches, ref_nmismatches;
  double prob;


#ifdef DEBUG8
  printf("Altsplice_select_qend for univdiagonal %u\n",this->medial_univdiagonal);
  Altsplice_print(this);
  printf("\n");
#endif

#ifdef CHECK_ASSERTIONS
  check_ascending(this->univdiagonals,/*order_medial_qpos*/NULL,this->nunivdiagonals);
#endif

  best_nmatches = 0;
  for (i = 0; i < this->nunivdiagonals; i++) {
    if (this->distal_nmismatches[i] < 0) {
      this->distal_nmismatches[i] =
	Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress,
					      /*univdiagonal*/this->univdiagonals[i],querylength,
					      /*pos5*/this->splice_qpos,/*pos3*/this->distal_trimpos[i],
					      plusp,genestrand);
    }
    debug8(printf(" Outer: (nmatches (%d-%d)-%d, prob %f+%f)\n",
		  this->distal_trimpos[i],this->splice_qpos,this->distal_nmismatches[i],
		  this->medial_prob,this->distal_probs[i]));
    if ((nmatches = (this->distal_trimpos[i] - this->splice_qpos) - this->distal_nmismatches[i]) > best_nmatches) {
      best_nmatches = nmatches;
    }
  }
  
  debug8(printf("best_nmatches %d\n",best_nmatches));

  this->nmismatches_known_p = true;

  ncandidates = 0;
  if (best_nmatches > 0) {
    /* Count good solutions */
    for (i = 0; i < this->nunivdiagonals; i++) {
      if ((this->distal_trimpos[i] - this->splice_qpos) - this->distal_nmismatches[i] == best_nmatches &&
#if 0
	  Splice_accept_p(this->splice_qpos,querylength,this->medial_prob,this->distal_probs[i],
			  /*univdiagonal_i*/this->medial_univdiagonal,
			  /*univdiagonal_j*/this->univdiagonals[i],
			  this->medial_support,this->medial_nmismatches,
			  this->distal_support[i],this->distal_nmismatches[i],
			  plusp_sense_forward_p,this->innerp) == true
#else
	  (prob = this->distal_probs[i]) >= GOOD_DISTAL_PROB
#endif	  
	  ) {
	debug8(printf("%i is a candidate\n",i));
	if (ncandidates++ == 0) {
	  /* Select shortest intron with a good distal prob */
	  besti = i;
	}
      }
    }
  }
  debug8(printf("ncandidates: %d\n",ncandidates));

  if (ncandidates == 0) {
    debug8(printf("No candidate found => returning false\n"));
    return -1;

  } else if (ncandidates > 1) {
    debug8(printf("Multiple candidates found => returning false\n"));
    return -1;

  } else {
    return besti;
  }
}



int
Altsplice_select_qstart (T this, int querylength, Compress_T query_compress, bool plusp, int genestrand) {
  int bestj = -1, j;
  int ncandidates;
  int best_nmatches, nmatches, ref_nmismatches;
  double prob;

#ifdef DEBUG8
  printf("Altsplice_select_qstart for univdiagonal %u\n",this->medial_univdiagonal);
  Altsplice_print(this);
  printf("\n");
#endif

#ifdef CHECK_ASSERTIONS
  check_descending(this->univdiagonals,/*order_medial_qpos*/NULL,this->nunivdiagonals);
#endif

  best_nmatches = 0;
  for (j = 0; j < this->nunivdiagonals; j++) {
    if (this->distal_nmismatches[j] < 0) {
      this->distal_nmismatches[j] =
	Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress,
					      /*univdiagonal*/this->univdiagonals[j],querylength,
					      /*pos5*/this->distal_trimpos[j],/*pos3*/this->splice_qpos,
					      plusp,genestrand);
    }
    debug8(printf(" Outer: (nmatches (%d-%d)-%d, prob %f+%f)\n",
		  this->splice_qpos,this->distal_trimpos[j],this->distal_nmismatches[j],
		  this->medial_prob,this->distal_probs[j]));
    if ((nmatches = (this->splice_qpos - this->distal_trimpos[j]) - this->distal_nmismatches[j]) > best_nmatches) {
	best_nmatches = nmatches;
    }
  }

  this->nmismatches_known_p = true;

  debug8(printf("best_nmatches %d\n",best_nmatches));

  ncandidates = 0;
  if (best_nmatches > 0) {
    /* Count good solutions */
    for (j = 0; j < this->nunivdiagonals; j++) {
      if ((this->splice_qpos - this->distal_trimpos[j]) - this->distal_nmismatches[j] == best_nmatches &&
#if 0
	  Splice_accept_p(this->splice_qpos,querylength,this->medial_prob,this->distal_probs[i],
			  /*univdiagonal_i*/this->univdiagonals[i],
			  /*univdiagonal_j*/this->medial_univdiagonal,
			  this->medial_support,this->medial_nmismatches,
			  this->distal_support[i],this->distal_nmismatches[i],
			  plusp,sense_forward_p,this->innerp) == true
#else
	  (prob = this->distal_probs[j]) >= GOOD_DISTAL_PROB
#endif
	  ) {
	debug8(printf("%d is a candidate\n",j));
	if (ncandidates++ == 0) {
	  /* Select shortest intron with a good distal prob */
	  bestj = j;
	}
      }
    }
  }
  debug8(printf("ncandidates: %d\n",ncandidates));

  if (ncandidates == 0) {
    debug8(printf("No candidate found => returning false\n"));
    return -1;

  } else if (ncandidates > 1) {
    debug8(printf("Multiple candidates found => returning false\n"));
    return -1;
	  
  } else {
    return bestj;
  }
}


void
Altsplice_setup (int max_insertlength_in, Genomebits_T genomebits_in,
		 Genomebits_T genomebits_alt_in) {

  max_insertlength = max_insertlength_in;
  genomebits = genomebits_in;
  genomebits_alt = genomebits_alt_in;

  return;
}
