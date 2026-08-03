/* Code modified from Sux package by Sebastiano Vigna,
   file bits/SimpleSelectHalf.hpp */

#include "select64-ones.h"

#include <stdio.h>
#include <stdlib.h>
#include "assert.h"
#include "mem.h"
#include "select64-common.h"

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#define log2_ones_per_inventory 10
#define ones_per_inventory (1 << log2_ones_per_inventory)
#define ones_per_inventory_mask (ones_per_inventory - 1)
#define log2_longwords_per_subinventory 2
#define longwords_per_subinventory (1 << log2_longwords_per_subinventory)
#define log2_ones_per_sub64 (log2_ones_per_inventory - log2_longwords_per_subinventory)
#define ones_per_sub64 (1 << log2_ones_per_sub64)
#define ones_per_sub64_mask (ones_per_sub64 - 1)
#define log2_ones_per_sub16 (log2_ones_per_sub64 - 2)
#define ones_per_sub16 (1 << log2_ones_per_sub16)
#define ones_per_sub16_mask (ones_per_sub16 - 1)


/** A simple Select implementation based on a two-level inventory,
 * and wired for approximately the same number of zeros and ones.
 *
 * The constructors of this class only store a reference
 * to a provided bit vector. Should the content of the
 * bit vector change, the results will be unpredictable.
 *
 * This implementation has been specifically developed to be used
 * with EliasFano.
 */


#define T Select64_ones_T
struct T {
  uint64_t *bits;
  int64_t *inventory;

  uint64_t nwords;
  uint64_t inventory_size;
  uint64_t num_ones;
};


void
Select64_ones_free (T *old) {
  /* FREE((*old)->bits); -- holds just a pointer */
  FREE((*old)->inventory);
  FREE(*old);
  return;
}


/* bits: a bit vector of 64-bit words */
/* nbits: the length (in bits) of the bit vector */

T
Select64_ones_new (uint64_t *const bits, const uint64_t nbits) {
  T new = (T) MALLOC(sizeof(*new));

  uint64_t c, i, d, j;
  uint64_t mask;
  uint64_t exact = 0, start, span, inventory_index;
  int offset;

  debug(printf("Select_ones_new called with nbits %d.  First word: %016llx\n",nbits,bits[0]));

  new->bits = bits;
  new->nwords = (nbits + 63) / 64;

  /* Init rank/select structure */
  c = 0;
  for (i = 0; i < new->nwords; i++) {
    c += __builtin_popcountll(bits[i]);
  }
  new->num_ones = c;
  assert(c <= nbits);

  new->inventory_size = (c + ones_per_inventory - 1) / ones_per_inventory;

#ifdef DEBUG
  printf("Number of bits: %llu.  Number of ones: %llu (%.2f%%)\n",
	 nbits, c, (c * 100.0) / nbits);
  printf("Ones per inventory: %d.  Ones per sub 64: %d.  sub 16: %d\n",
	 ones_per_inventory, ones_per_sub64, ones_per_sub16);
#endif

  new->inventory = (int64_t *) CALLOC((new->inventory_size * (longwords_per_subinventory + 1) + 1), sizeof(int64_t));

  d = 0;
  mask = ones_per_inventory - 1;

  /* First phase: we build an inventory for each one out of ones_per_inventory. */
  for (i = 0; i < new->nwords; i++) {
    for (j = 0; j < 64; j++) {
      if (i * 64 + j >= nbits) {
	break;
      }
      if (bits[i] & 1ULL << j) {
	if ((d & mask) == 0) {
	  new->inventory[(d >> log2_ones_per_inventory) * (longwords_per_subinventory + 1)] = i * 64 + j;
	}
	d++;
      }
    }
  }

  assert(c == d);
  new->inventory[new->inventory_size * (longwords_per_subinventory + 1)] = nbits;

#ifdef DEBUG
  printf("Inventory entries filled: %llu\n", new->inventory_size + 1);
#endif

  uint16_t *p16;
  int64_t *p64;

  d = 0;
  exact = 0;

  for (i = 0; i < new->nwords; i++) {
    for (j = 0; j < 64; j++) {
      if (i * 64 + j >= nbits) {
	break;
      }
      if (bits[i] & 1ULL << j) {
	if ((d & mask) == 0) {
	  inventory_index = (d >> log2_ones_per_inventory) * (longwords_per_subinventory + 1);
	  start = new->inventory[inventory_index];
	  span = new->inventory[inventory_index + longwords_per_subinventory + 1] - start;
	  if (span > (1 << 16)) {
	    new->inventory[inventory_index] = -new->inventory[inventory_index] - 1;
	  }
	  offset = 0;
	  p64 = &(new->inventory[inventory_index + 1]);
	  p16 = (uint16_t *) p64;
	}

	if (span < (1 << 16)) {
	  assert(i * 64 + j - start <= (1 << 16));
	  if ((d & ones_per_sub16_mask) == 0) {
	    assert(offset < longwords_per_subinventory * 4);
	    p16[offset++] = i * 64 + j - start;
	  }
	} else {
	  if ((d & ones_per_sub64_mask) == 0) {
	    assert(offset < longwords_per_subinventory);
	    p64[offset++] = i * 64 + j - start;
	    exact++;
	  }
	}

	d++;
      }
    }
  }

#ifdef DEBUG
  printf("Exact entries: %llu\n", exact);
  for (i = 0; i < new->inventory_size * (longwords_per_subinventory + 1) + 1; i++) {
    printf("Inventory %d %llu\n",i,new->inventory[i]);
  }
#endif

  return new;
}


uint64_t
Select64_ones_select (T this, const uint64_t rank) {
  debug(printf("Selecting %llu\n...", rank));

  debug(printf("Entering Select64_ones_select with rank, or %llu th one vs num_ones %llu\n...",
	       rank,this->num_ones));

  assert(rank < this->num_ones);

  const uint64_t inventory_index = rank >> log2_ones_per_inventory;
  assert(inventory_index <= this->inventory_size);
  const int64_t *inventory_start =
    &(this->inventory[(inventory_index << log2_longwords_per_subinventory) + inventory_index]);

  const int64_t inventory_rank = *inventory_start;
  const int subrank = rank & ones_per_inventory_mask;
#ifdef DEBUG
  printf("Rank: %llu  inventory index: %llu  inventory rank: %ld  subrank: %d\n",
	 rank, inventory_index, inventory_rank, subrank);
#endif

  uint64_t start;
  int residual;

  if (inventory_rank >= 0) {
    start = inventory_rank + ((uint16_t *)(inventory_start + 1))[subrank >> log2_ones_per_sub16];
    debug(printf("subrank %d >> log2_ones_per_sub16 %d is %d\n",
		 subrank,log2_ones_per_sub16,subrank >> log2_ones_per_sub16));
    residual = subrank & ones_per_sub16_mask;
  } else {
    assert((subrank >> log2_ones_per_sub64) < longwords_per_subinventory);
    start = -inventory_rank - 1 + *(inventory_start + 1 + (subrank >> log2_ones_per_sub64));
    residual = subrank & ones_per_sub64_mask;
  }

#ifdef DEBUG
  printf("Differential; start: %llu  residual: %d\n", start, residual);
  if (residual == 0) puts("No residual; returning start");
#endif

  if (residual == 0) return start;

  uint64_t word_index = start / 64;
  uint64_t word = this->bits[word_index] & -1ULL << start % 64;

  for (;;) {
    const int bit_count = __builtin_popcountll(word);
    if (residual < bit_count) {
      break;
    }
    word = this->bits[++word_index];
    residual -= bit_count;
  }

  return word_index * 64 + select64(word, residual);
}


uint64_t
Select64_ones_next (T this, const uint64_t rank, uint64_t *const next) {

  const uint64_t s = Select64_ones_select(this,rank);
  int curr = s / 64;

  uint64_t window = this->bits[curr] & (~0ULL << (s % 64));
  window &= window - 1;

  while (window == 0) {
    window = this->bits[++curr];
  }
  *next = curr * 64 + __builtin_ctzll(window);
  
  return s;
}

#if 0
/** Returns an estimate of the size (in bits) of this structure. */
size_t bitCount() const { return inventory.bitCount() - sizeof(inventory) * 8 + sizeof(*this) * 8; };
};
#endif

