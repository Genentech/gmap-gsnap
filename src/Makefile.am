
# Copied from ../util so "gmap_setup -B src" works during "make check"
CLEANFILES = fa_coords gmap_process

EXTRA_DIST = simde debug_groups.h

# Include master.c and master.h, so they get included in distribution
# EXTRA_DIST = mpidebug.c mpidebug.h master.c master.h

# -DUTILITYP=1 needed for iit_store, iit_dump, iit_get, gmapindex, and
# get-genome so they can handle both small and large genomes at run
# time

#lib_LTLIBRARIES = libgmap-@LIBGMAP_API_VERSION@.la
#if MAKE_LIB
#lib_LTLIBRARIES = libgmap.la
#endif


/* Obsolete: sarray-search */
/* Obsolete: genome128 genome128_hr genome128_consec */
/* Obsolete: spanningelt, indexdb_hr */
/* Obsolete: localdb, genome_hr, gdiag */
/* Obsolete: merge-records-heap, merge-records-simd */
/* Obsolete: uniqscan, uniqscanl */
/* Obsolete: cigar, mdprint, samprint */
/* Obsolete: stage3hr, substring, terminal */
/* Obsolete: ladder */
/* TO_FIX: distant-rna, distant-dna */


#if MAKE_LIB
#include_HEADERS = fopen.h bool.h types.h univcoord.h separator.h comp.h \
# except.h assert.h mem.h \
# intlist.h list.h \
# doublelist.h \
# littleendian.h bigendian.h \
# interval.h uintlist.h uint8list.h \
# iitdef.h iit-read.h iit-write.h parserange.h \
# univinterval.h iit-read-univ.h \
# table.h tableuint.h uinttable.h uinttableuint.h \
# stopwatch.h semaphore.h access.h \
# chrom.h filestring.h \
# md5.h complement.h bzip2.h sequence.h \
# genomicpos.h \
# bitpack64-read.h bitpack64-readtwo.h \
# maxent_hr.h \
# popcount.h genome128_hr.h genome128_consec.h \
# compress.h \
# mode.h chrnum.h genome.h samflags.h \
# getline.h datadir.h
#endif


# velocity-counts has been replaced by gcount-bulk and gcount-sc

# Not yet ready for release: gcount-bulk, gcount-sc, gexact, gfilter, compare2truth, sam_sort
bin_PROGRAMS = cpuid gmap gmapl get-genome gmapindex indexdb_cat \
               iit_store iit_get iit_dump \
               gsnap gsnapl \
               snpindex cmetindex atoiindex trindex

bin_PROGRAMS += gmap.nosimd
bin_PROGRAMS += gmapl.nosimd
bin_PROGRAMS += gsnap.nosimd
bin_PROGRAMS += gsnapl.nosimd
#bin_PROGRAMS += gexact.nosimd
#bin_PROGRAMS += gfilter.nosimd
if MAKE_ARM
  bin_PROGRAMS += gmap.arm
  bin_PROGRAMS += gmapl.arm
  bin_PROGRAMS += gsnap.arm
  bin_PROGRAMS += gsnapl.arm
#  bin_PROGRAMS += gexact.arm
#  bin_PROGRAMS += gfilter.arm
endif
if MAKE_SSE2
  bin_PROGRAMS += gmap.sse2
  bin_PROGRAMS += gmapl.sse2
  bin_PROGRAMS += gsnap.sse2
  bin_PROGRAMS += gsnapl.sse2
#  bin_PROGRAMS += gexact.sse2
#  bin_PROGRAMS += gfilter.sse2
endif
if MAKE_SSSE3
  bin_PROGRAMS += gmap.ssse3
  bin_PROGRAMS += gmapl.ssse3
  bin_PROGRAMS += gsnap.ssse3
  bin_PROGRAMS += gsnapl.ssse3
#  bin_PROGRAMS += gexact.ssse3
#  bin_PROGRAMS += gfilter.ssse3
endif
if MAKE_SSE4_1
  bin_PROGRAMS += gmap.sse41
  bin_PROGRAMS += gmapl.sse41
  bin_PROGRAMS += gsnap.sse41
  bin_PROGRAMS += gsnapl.sse41
#  bin_PROGRAMS += gexact.sse41
#  bin_PROGRAMS += gfilter.sse41
endif
if MAKE_SSE4_2
  bin_PROGRAMS += gmap.sse42
  bin_PROGRAMS += gmapl.sse42
  bin_PROGRAMS += gsnap.sse42
  bin_PROGRAMS += gsnapl.sse42
#  bin_PROGRAMS += gexact.sse42
#  bin_PROGRAMS += gfilter.sse42
endif
if MAKE_AVX2
  bin_PROGRAMS += gmap.avx2
  bin_PROGRAMS += gmapl.avx2
  bin_PROGRAMS += gsnap.avx2
  bin_PROGRAMS += gsnapl.avx2
#  bin_PROGRAMS += gexact.avx2
#  bin_PROGRAMS += gfilter.avx2
endif
if MAKE_AVX512
  bin_PROGRAMS += gmap.avx512
  bin_PROGRAMS += gmapl.avx512
  bin_PROGRAMS += gsnap.avx512
  bin_PROGRAMS += gsnapl.avx512
#  bin_PROGRAMS += gexact.avx512
#  bin_PROGRAMS += gfilter.avx512
endif
if MAKE_AVX512BW
  bin_PROGRAMS += gmap.avx512bw
  bin_PROGRAMS += gmapl.avx512bw
  bin_PROGRAMS += gsnap.avx512bw
  bin_PROGRAMS += gsnapl.avx512bw
#  bin_PROGRAMS += gexact.avx512bw
#  bin_PROGRAMS += gfilter.avx512bw
endif


#LIBGMAP_LA_FILES = bool.h types.h univcoord.h separator.h comp.h \
# except.c except.h assert.c assert.h mem.c mem.h \
# intlist.c intlist.h list.c list.h \
# doublelist.c doublelist.h \
# littleendian.c littleendian.h bigendian.c bigendian.h \
# interval.c interval.h uintlist.c uintlist.h uint8list.c uint8list.h \
# iit-read.c iit-read.h iit-write.c iit-write.h parserange.c parserange.h \
# univinterval.c univinterval.h iit-read-univ.c iit-read-univ.h \
# stopwatch.c stopwatch.h semaphore.c semaphore.h access.c access.h \
# table.c table.h tableuint.c tableuint.h uinttable.c uinttable.h uinttableuint.c uinttableuint.h\
# chrom.c chrom.h filestring.c filestring.h \
# md5.c md5.h complement.h bzip2.c bzip2.h fopen.c fopen.h sequence.c sequence.h \
# genomicpos.c genomicpos.h \
# bitpack64-read.c bitpack64-read.h bitpack64-readtwo.c bitpack64-readtwo.h \
# maxent_hr.c maxent_hr.h \
# popcount.c popcount.h genome128_hr.c genome128_hr.h genome128_consec.c genome128_consec.h \
# compress.c compress.h \
# mode.h chrnum.c chrnum.h genome.c genome.h \
# getline.c getline.h datadir.c datadir.h
#
#libgmap_CC = $(PTHREAD_CC)
#libgmap_la_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) $(POPCNT_CFLAGS) -fPIC -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\"
#libgmap_la_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) -version-info $(LIBGMAP_SO_VERSION)
#libgmap_la_LIBADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
#dist_libgmap_la_SOURCES = $(LIBGMAP_LA_FILES)



# EF64_T_FILES = types.h \
#                select64-common.c select64-common.h select64-ones.c select64-ones.h \
#                select64-zeroes.c select64-zeroes.h \
#                except.c except.h assert.c assert.h mem.c mem.h \
#                uintlist.c uintlist.h ef64.c ef64.h ef64.t.c
# ef64_t_CC = $(PTHREAD_CC)
# ef64_t_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS)
# ef64_t_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
# ef64_t_LDADD = $(PTHREAD_LIBS)
# dist_ef64_t_SOURCES = $(EF64_T_FILES)


CPUID_FILES = bool.h cpuid.c cpuid.h
cpuid_CC = $(CC)
cpuid_CFLAGS = $(AM_CFLAGS) -DMAIN=1
cpuid_LDFLAGS = $(AM_LDFLAGS)
cpuid_LDADD = 
dist_cpuid_SOURCES = $(CPUID_FILES)


GMAP_SELECT_FILES = bool.h cpuid.c cpuid.h gmap_select.c
gmap_CC = $(CC)
gmap_CFLAGS = $(AM_CFLAGS)
gmap_LDFLAGS = $(AM_LDFLAGS)
gmap_LDADD = 
dist_gmap_SOURCES = $(GMAP_SELECT_FILES)

GMAPL_SELECT_FILES = bool.h cpuid.c cpuid.h gmapl_select.c
gmapl_CC = $(CC)
gmapl_CFLAGS = $(AM_CFLAGS)
gmapl_LDFLAGS = $(AM_LDFLAGS)
gmapl_LDADD = 
dist_gmapl_SOURCES = $(GMAPL_SELECT_FILES)

GSNAP_SELECT_FILES = bool.h cpuid.c cpuid.h gsnap_select.c
gsnap_CC = $(CC)
gsnap_CFLAGS = $(AM_CFLAGS)
gsnap_LDFLAGS = $(AM_LDFLAGS)
gsnap_LDADD = 
dist_gsnap_SOURCES = $(GSNAP_SELECT_FILES)

GSNAPL_SELECT_FILES = bool.h cpuid.c cpuid.h gsnapl_select.c
gsnapl_CC = $(CC)
gsnapl_CFLAGS = $(AM_CFLAGS)
gsnapl_LDFLAGS = $(AM_LDFLAGS)
gsnapl_LDADD = 
dist_gsnapl_SOURCES = $(GSNAPL_SELECT_FILES)

GEXACT_SELECT_FILES = bool.h cpuid.c cpuid.h gexact_select.c
gexact_CC = $(CC)
gexact_CFLAGS = $(AM_CFLAGS)
gexact_LDFLAGS = $(AM_LDFLAGS)
gexact_LDADD = 
dist_gexact_SOURCES = $(GEXACT_SELECT_FILES)

GFILTER_SELECT_FILES = bool.h cpuid.c cpuid.h gfilter_select.c
gfilter_CC = $(CC)
gfilter_CFLAGS = $(AM_CFLAGS)
gfilter_LDFLAGS = $(AM_LDFLAGS)
gfilter_LDADD = 
dist_gfilter_SOURCES = $(GFILTER_SELECT_FILES)


GMAP_FILES = simd.h \
 bool.h types.h univcoord.h separator.h comp.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h uintlist.c uintlist.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 univinterval.c univinterval.h interval.c interval.h \
 stopwatch.c stopwatch.h semaphore.c semaphore.h access.c access.h \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iitdef.h iit-read.c iit-read.h \
 md5.c md5.h complement.h bzip2.c bzip2.h fopen.c fopen.h sequence.c sequence.h reader.c reader.h \
 genomicpos.c genomicpos.h compress.c compress.h compress-write.c compress-write.h \
 gbuffer.c gbuffer.h genome.c genome.h \
 popcount.c popcount.h dinucl_bits.c dinucl_bits.h genome_canonical.c genome_canonical.h \
 genome-write.c genome-write.h \
 bitpack64-read.c bitpack64-read.h bitpack64-readtwo.c bitpack64-readtwo.h \
 filesuffix.h signature.c signature.h indexdbdef.h indexdb.c indexdb.h \
 oligo.c oligo.h block.c block.h \
 chrom.c chrom.h segmentpos.c segmentpos.h \
 chrnum.c chrnum.h uinttable_rh.c uinttable_rh.h gregion.c gregion.h \
 matchdef.h match.c match.h matchpool.c matchpool.h \
 diagnostic.c diagnostic.h stage1.c stage1.h \
 diagdef.h diag.c diag.h diagpool.c diagpool.h \
 cmet.c cmet.h atoi.c atoi.h \
 orderstat.c orderstat.h oligoindex_hr.c oligoindex_hr.h \
 scores.h intron.c intron.h maxent.c maxent.h maxent_hr.c maxent_hr.h samflags.h pairdef.h pair.c pair.h \
 pairpool.c pairpool.h cellpool.c cellpool.h stage2.c stage2.h \
 doublelist.c doublelist.h smooth.c smooth.h \
 splicestringpool.c splicestringpool.h splicetrie_build.c splicetrie_build.h splicetrie.c splicetrie.h \
 boyer-moore.c boyer-moore.h \
 dynprog.c dynprog.h dynprog_simd.c dynprog_simd.h \
 dynprog_single.c dynprog_single.h dynprog_genome.c dynprog_genome.h dynprog_cdna.c dynprog_cdna.h dynprog_end.c dynprog_end.h \
 translation.c translation.h \
 pbinom.c pbinom.h changepoint.c changepoint.h sense.h fastlog.h stage3.c stage3.h \
 request.c request.h result.c result.h outputtype.h output.c output.h \
 inbuffer.c inbuffer.h samheader.c samheader.h printbuffer.c printbuffer.h outbuffer.c outbuffer.h \
 chimera.c chimera.h datadir.c datadir.h parserange.c parserange.h \
 getline.c getline.h getopt.c getopt1.c getopt.h gmap.c


# Note: dist_ commands get read by bootstrap, and don't follow the flags
# Previously had POPCNT_CFLAGS for nosimd versions, but not needed

gmap_nosimd_CC = $(PTHREAD_CC)
gmap_nosimd_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\"
gmap_nosimd_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gmap_nosimd_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gmap_nosimd_SOURCES = $(GMAP_FILES)

gmap_arm_CC = $(PTHREAD_CC)
gmap_arm_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DHAVE_ARM=1 -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 $(SIMD_SSE4_2_CFLAGS)
gmap_arm_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gmap_arm_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gmap_arm_SOURCES = $(GMAP_FILES)

gmap_sse2_CC = $(PTHREAD_CC)
gmap_sse2_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" $(POPCNT_CFLAGS) -DHAVE_SSE2=1 $(SIMD_SSE2_CFLAGS)
gmap_sse2_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gmap_sse2_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gmap_sse2_SOURCES = $(GMAP_FILES)

gmap_ssse3_CC = $(PTHREAD_CC)
gmap_ssse3_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 $(SIMD_SSSE3_CFLAGS)
gmap_ssse3_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gmap_ssse3_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gmap_ssse3_SOURCES = $(GMAP_FILES)

gmap_sse41_CC = $(PTHREAD_CC)
gmap_sse41_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 $(SIMD_SSE4_1_CFLAGS)
gmap_sse41_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gmap_sse41_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gmap_sse41_SOURCES = $(GMAP_FILES)

gmap_sse42_CC = $(PTHREAD_CC)
gmap_sse42_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 $(SIMD_SSE4_2_CFLAGS)
gmap_sse42_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gmap_sse42_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gmap_sse42_SOURCES = $(GMAP_FILES)

gmap_avx2_CC = $(PTHREAD_CC)
gmap_avx2_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 -DHAVE_AVX2=1 $(SIMD_AVX2_CFLAGS)
gmap_avx2_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gmap_avx2_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gmap_avx2_SOURCES = $(GMAP_FILES)

gmap_avx512_CC = $(PTHREAD_CC)
gmap_avx512_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 -DHAVE_AVX2=1 -DHAVE_AVX512=1 $(SIMD_AVX512_CFLAGS)
gmap_avx512_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gmap_avx512_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gmap_avx512_SOURCES = $(GMAP_FILES)

gmap_avx512bw_CC = $(PTHREAD_CC)
gmap_avx512bw_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 -DHAVE_AVX2=1 -DHAVE_AVX512=1 -DHAVE_AVX512BW=1 $(SIMD_AVX512_CFLAGS)
gmap_avx512bw_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gmap_avx512bw_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gmap_avx512bw_SOURCES = $(GMAP_FILES)



GMAPL_FILES = simd.h \
 bool.h types.h univcoord.h separator.h comp.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h uintlist.c uintlist.h uint8list.c uint8list.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 univinterval.c univinterval.h interval.c interval.h \
 stopwatch.c stopwatch.h semaphore.c semaphore.h access.c access.h \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iitdef.h iit-read.c iit-read.h \
 md5.c md5.h complement.h bzip2.c bzip2.h fopen.c fopen.h sequence.c sequence.h reader.c reader.h \
 genomicpos.c genomicpos.h compress.c compress.h compress-write.c compress-write.h \
 gbuffer.c gbuffer.h genome.c genome.h \
 popcount.c popcount.h dinucl_bits.c dinucl_bits.h genome_canonical.c genome_canonical.h \
 genome-write.c genome-write.h \
 bitpack64-read.c bitpack64-read.h bitpack64-readtwo.c bitpack64-readtwo.h \
 filesuffix.h signature.c signature.h indexdbdef.h indexdb.c indexdb.h \
 oligo.c oligo.h block.c block.h \
 chrom.c chrom.h segmentpos.c segmentpos.h \
 chrnum.c chrnum.h uinttable_rh.c uinttable_rh.h uint8table_rh.c uint8table_rh.h gregion.c gregion.h \
 matchdef.h match.c match.h matchpool.c matchpool.h \
 diagnostic.c diagnostic.h stage1.c stage1.h \
 diagdef.h diag.c diag.h diagpool.c diagpool.h \
 cmet.c cmet.h atoi.c atoi.h \
 orderstat.c orderstat.h oligoindex_hr.c oligoindex_hr.h \
 scores.h intron.c intron.h maxent.c maxent.h maxent_hr.c maxent_hr.h samflags.h pairdef.h pair.c pair.h \
 pairpool.c pairpool.h cellpool.c cellpool.h stage2.c stage2.h \
 doublelist.c doublelist.h smooth.c smooth.h \
 splicestringpool.c splicestringpool.h splicetrie_build.c splicetrie_build.h splicetrie.c splicetrie.h \
 boyer-moore.c boyer-moore.h \
 dynprog.c dynprog.h dynprog_simd.c dynprog_simd.h \
 dynprog_single.c dynprog_single.h dynprog_genome.c dynprog_genome.h dynprog_cdna.c dynprog_cdna.h dynprog_end.c dynprog_end.h \
 translation.c translation.h \
 pbinom.c pbinom.h changepoint.c changepoint.h sense.h fastlog.h stage3.c stage3.h \
 request.c request.h result.c result.h outputtype.h output.c output.h \
 inbuffer.c inbuffer.h samheader.c samheader.h printbuffer.c printbuffer.h outbuffer.c outbuffer.h \
 chimera.c chimera.h datadir.c datadir.h parserange.c parserange.h \
 getline.c getline.h getopt.c getopt1.c getopt.h gmap.c


# Note: dist_ commands get read by bootstrap, and don't follow the flags

gmapl_nosimd_CC = $(PTHREAD_CC)
gmapl_nosimd_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DLARGE_GENOMES=1
gmapl_nosimd_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gmapl_nosimd_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gmapl_nosimd_SOURCES = $(GMAPL_FILES)

gmapl_arm_CC = $(PTHREAD_CC)
gmapl_arm_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DLARGE_GENOMES=1 -DHAVE_ARM=1 -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 $(SIMD_SSE4_2_CFLAGS)
gmapl_arm_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gmapl_arm_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gmapl_arm_SOURCES = $(GMAPL_FILES)

gmapl_sse2_CC = $(PTHREAD_CC)
gmapl_sse2_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DLARGE_GENOMES=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 $(SIMD_SSE2_CFLAGS)
gmapl_sse2_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gmapl_sse2_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gmapl_sse2_SOURCES = $(GMAPL_FILES)

gmapl_ssse3_CC = $(PTHREAD_CC)
gmapl_ssse3_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DLARGE_GENOMES=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 $(SIMD_SSSE3_CFLAGS)
gmapl_ssse3_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gmapl_ssse3_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gmapl_ssse3_SOURCES = $(GMAPL_FILES)

gmapl_sse41_CC = $(PTHREAD_CC)
gmapl_sse41_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DLARGE_GENOMES=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 $(SIMD_SSE4_1_CFLAGS)
gmapl_sse41_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gmapl_sse41_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gmapl_sse41_SOURCES = $(GMAPL_FILES)

gmapl_sse42_CC = $(PTHREAD_CC)
gmapl_sse42_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DLARGE_GENOMES=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 $(SIMD_SSE4_2_CFLAGS)
gmapl_sse42_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gmapl_sse42_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gmapl_sse42_SOURCES = $(GMAPL_FILES)

gmapl_avx2_CC = $(PTHREAD_CC)
gmapl_avx2_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DLARGE_GENOMES=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 -DHAVE_AVX2=1 $(SIMD_AVX2_CFLAGS)
gmapl_avx2_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gmapl_avx2_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gmapl_avx2_SOURCES = $(GMAPL_FILES)

gmapl_avx512_CC = $(PTHREAD_CC)
gmapl_avx512_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DLARGE_GENOMES=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 -DHAVE_AVX2=1 -DHAVE_AVX512=1 $(SIMD_AVX512_CFLAGS)
gmapl_avx512_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gmapl_avx512_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gmapl_avx512_SOURCES = $(GMAPL_FILES)

gmapl_avx512bw_CC = $(PTHREAD_CC)
gmapl_avx512bw_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DLARGE_GENOMES=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 -DHAVE_AVX2=1 -DHAVE_AVX512=1 -DHAVE_AVX512BW=1 $(SIMD_AVX512_CFLAGS)
gmapl_avx512bw_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gmapl_avx512bw_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gmapl_avx512bw_SOURCES = $(GMAPL_FILES)


GSNAP_FILES = simd.h \
 bool.h types.h univcoord.h comp.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h uintlist.c uintlist.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 univinterval.c univinterval.h interval.c interval.h \
 stopwatch.c stopwatch.h semaphore.c semaphore.h access.c access.h \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iitdef.h iit-read.c iit-read.h \
 complement.h bzip2.c bzip2.h reader.c reader.h \
 genomicpos.c genomicpos.h compress.c compress.h \
 genome.c genome.h transcriptome.c transcriptome.h \
 popcount.c popcount.h \
 genomebits.c genomebits.h genomebits_consec.c genomebits_consec.h genomebits_count.c genomebits_count.h \
 genomebits_kmer.c genomebits_kmer.h genomebits_indel.c genomebits_indel.h \
 genomebits_decode.c genomebits_decode.h genomebits_mismatches.c genomebits_mismatches.h \
 genomebits_trim.c genomebits_trim.h \
 dinucl_bits.c dinucl_bits.h genome_sites.c genome_sites.h \
 bitpack64-read.c bitpack64-read.h bitpack64-readtwo.c bitpack64-readtwo.h \
 mergeinfo.c mergeinfo.h merge-method.h \
 merge-diagonals-heap-uint4.c merge-diagonals-heap-uint4.h \
 merge-uint4.c merge-uint4.h merge-diagonals-simd-uint4.c merge-diagonals-simd-uint4.h \
 record.h \
 filesuffix.h signature.c signature.h indexdbdef.h indexdb.c indexdb.h \
 localdb-read.c localdb-read.h \
 oligo.c oligo.h \
 chrnum.c chrnum.h \
 maxent_hr.c maxent_hr.h samflags.h \
 mapq.c mapq.h fopen.c fopen.h shortread.c shortread.h junction.c junction.h \
 concordance.c concordance.h \
 simplepair.c simplepair.h \
 orderstat.c orderstat.h \
 intlistpool.c intlistpool.h uintlistpool.c uintlistpool.h listpool.c listpool.h \
 hitlistpool.c hitlistpool.h pathpool.c pathpool.h \
 auxinfopool.c auxinfopool.h pathstorepool.c pathstorepool.h \
 intron.c intron.h \
 gbuffer.c gbuffer.h \
 sense.h fastlog.h \
 select64-common.c select64-common.h select64-ones.c select64-ones.h select64-zeroes.c select64-zeroes.h ef64.c ef64.h \
 nr-x.c nr-x.h knownsplicing.c knownsplicing.h knownindels.c knownindels.h \
 splice.c splice.h indel.c indel.h \
 spliceendsgen.c spliceendsgen.h spliceends.c spliceends.h \
 bitpack64-access.c bitpack64-access.h \
 trdiagdef.h trdiag.c trdiag.h trdiagpool.c trdiagpool.h tr-extension-search.c tr-extension-search.h \
 univdiagdef.h univdiag.c univdiag.h univdiagpool.c univdiagpool.h sedgesort.c sedgesort.h \
 auxinfo.c auxinfo.h pathstore.c pathstore.h \
 altsplice.c altsplice.h \
 path.c path.h pathpair.c pathpair.h path-solve.c path-solve.h path-fusion.c path-fusion.h \
 path-trim.c path-trim.h path-eval.c path-eval.h pathpair-eval.c pathpair-eval.h \
 path-print-alignment.c path-print-alignment.h path-print-m8.c path-print-m8.h path-print-sam.c path-print-sam.h \
 path-learn.c path-learn.h \
 method.c method.h \
 doublelist.c doublelist.h bitvector.c bitvector.h \
 transcriptpool.c transcriptpool.h exon.c exon.h transcript.c transcript.h \
 trpathpool.c trpathpool.h trpath.c trpath.h trpath-solve.c trpath-solve.h trpath-convert.c trpath-convert.h \
 trpath-eval.c trpath-eval.h \
 transcriptome-search.c transcriptome-search.h \
 vectorpool.c vectorpool.h repair.c repair.h transcript-remap.c transcript-remap.h \
 transcript-velocity.c transcript-velocity.h \
 intersect-small.c intersect-small.h intersect-uint2.c intersect-uint2.h \
 intersect-approx-uint4.c intersect-approx-uint4.h \
 intersect-approx-indices-uint4.c intersect-approx-indices-uint4.h \
 intersect-wdups-indices.c intersect-wdups-indices.h \
 intersect-lower-small.c intersect-lower-small.h \
 intersect-higher-small.c intersect-higher-small.h \
 intersect-indices-small.c intersect-indices-small.h \
 intersect-indices2-small.c intersect-indices2-small.h \
 intersect-concordance-uint4.c intersect-concordance-uint4.h \
 kmer-search.c kmer-search.h \
 extension-search.c extension-search.h \
 repetitive.c repetitive.h \
 uint8table.c uint8table.h uinttable.c uinttable.h \
 stage1hr.c stage1hr.h stage1hr-single.c stage1hr-single.h stage1hr-paired.c stage1hr-paired.h \
 uinttableuint.c uinttableuint.h \
 tableuint.c tableuint.h single-cell.c single-cell.h \
 request.c request.h resulthr.c resulthr.h outputtype.h output.c output.h \
 inbuffer.c inbuffer.h samheader.c samheader.h printbuffer.c printbuffer.h outbuffer.c outbuffer.h \
 uintpairtable.c uintpairtable.h splice-count.c splice-count.h support.c support.h \
 datadir.c datadir.h pass.h mode.h parserange.c parserange.h \
 getline.c getline.h getopt.c getopt1.c getopt.h gsnap.c


# Note: dist_ commands get read by bootstrap, and don't follow the flags
# Previously included -lrt for shm_open, but we are not calling that

gsnap_nosimd_CC = $(PTHREAD_CC)
gsnap_nosimd_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGSNAP=1
gsnap_nosimd_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gsnap_nosimd_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gsnap_nosimd_SOURCES = $(GSNAP_FILES)

gsnap_arm_CC = $(PTHREAD_CC)
gsnap_arm_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGSNAP=1 -DHAVE_ARM=1 -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 $(SIMD_SSE4_2_CFLAGS)
gsnap_arm_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gsnap_arm_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gsnap_arm_SOURCES = $(GSNAP_FILES)

gsnap_sse2_CC = $(PTHREAD_CC)
gsnap_sse2_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGSNAP=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 $(SIMD_SSE2_CFLAGS)
gsnap_sse2_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gsnap_sse2_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gsnap_sse2_SOURCES = $(GSNAP_FILES)

gsnap_ssse3_CC = $(PTHREAD_CC)
gsnap_ssse3_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGSNAP=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 $(SIMD_SSSE3_CFLAGS)
gsnap_ssse3_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gsnap_ssse3_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gsnap_ssse3_SOURCES = $(GSNAP_FILES)

gsnap_sse41_CC = $(PTHREAD_CC)
gsnap_sse41_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGSNAP=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 $(SIMD_SSE4_1_CFLAGS)
gsnap_sse41_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gsnap_sse41_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gsnap_sse41_SOURCES = $(GSNAP_FILES)

gsnap_sse42_CC = $(PTHREAD_CC)
gsnap_sse42_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGSNAP=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 $(SIMD_SSE4_2_CFLAGS)
gsnap_sse42_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gsnap_sse42_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gsnap_sse42_SOURCES = $(GSNAP_FILES)

gsnap_avx2_CC = $(PTHREAD_CC)
gsnap_avx2_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGSNAP=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 -DHAVE_AVX2=1 $(SIMD_AVX2_CFLAGS)
gsnap_avx2_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gsnap_avx2_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gsnap_avx2_SOURCES = $(GSNAP_FILES)

gsnap_avx512_CC = $(PTHREAD_CC)
gsnap_avx512_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGSNAP=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 -DHAVE_AVX2=1 -DHAVE_AVX512=1 $(SIMD_AVX512_CFLAGS)
gsnap_avx512_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gsnap_avx512_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gsnap_avx512_SOURCES = $(GSNAP_FILES)

gsnap_avx512bw_CC = $(PTHREAD_CC)
gsnap_avx512bw_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGSNAP=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 -DHAVE_AVX2=1 -DHAVE_AVX512=1 -DHAVE_AVX512BW=1 $(SIMD_AVX512_CFLAGS)
gsnap_avx512bw_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gsnap_avx512bw_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gsnap_avx512bw_SOURCES = $(GSNAP_FILES)




GSNAPL_FILES = simd.h \
 bool.h types.h univcoord.h comp.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h uintlist.c uintlist.h uint8list.c uint8list.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 univinterval.c univinterval.h interval.c interval.h \
 stopwatch.c stopwatch.h semaphore.c semaphore.h access.c access.h \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iitdef.h iit-read.c iit-read.h \
 complement.h bzip2.c bzip2.h reader.c reader.h \
 genomicpos.c genomicpos.h compress.c compress.h \
 genome.c genome.h transcriptome.c transcriptome.h \
 popcount.c popcount.h \
 genomebits.c genomebits.h genomebits_consec.c genomebits_consec.h genomebits_count.c genomebits_count.h \
 genomebits_kmer.c genomebits_kmer.h genomebits_indel.c genomebits_indel.h \
 genomebits_decode.c genomebits_decode.h genomebits_mismatches.c genomebits_mismatches.h \
 genomebits_trim.c genomebits_trim.h \
 dinucl_bits.c dinucl_bits.h genome_sites.c genome_sites.h \
 bitpack64-read.c bitpack64-read.h bitpack64-readtwo.c bitpack64-readtwo.h \
 mergeinfo.c mergeinfo.h merge-method.h \
 merge-diagonals-heap-uint8.c merge-diagonals-heap-uint8.h \
 merge-uint8.c merge-uint8.h merge-diagonals-simd-uint8.c merge-diagonals-simd-uint8.h \
 merge-diagonals-heap-uint4.c merge-diagonals-heap-uint4.h \
 merge-uint4.c merge-uint4.h merge-diagonals-simd-uint4.c merge-diagonals-simd-uint4.h \
 record.h \
 filesuffix.h signature.c signature.h indexdbdef.h indexdb.c indexdb.h \
 localdb-read.c localdb-read.h \
 oligo.c oligo.h \
 chrnum.c chrnum.h \
 maxent_hr.c maxent_hr.h samflags.h \
 mapq.c mapq.h fopen.c fopen.h shortread.c shortread.h junction.c junction.h \
 concordance.c concordance.h \
 simplepair.c simplepair.h \
 orderstat.c orderstat.h \
 intlistpool.c intlistpool.h uintlistpool.c uintlistpool.h uint8listpool.c uint8listpool.h \
 listpool.c listpool.h hitlistpool.c hitlistpool.h pathpool.c pathpool.h \
 auxinfopool.c auxinfopool.h pathstorepool.c pathstorepool.h \
 intron.c intron.h \
 gbuffer.c gbuffer.h \
 sense.h fastlog.h \
 select64-common.c select64-common.h select64-ones.c select64-ones.h select64-zeroes.c select64-zeroes.h ef64.c ef64.h \
 nr-x.c nr-x.h knownsplicing.c knownsplicing.h knownindels.c knownindels.h \
 splice.c splice.h indel.c indel.h \
 spliceendsgen.c spliceendsgen.h spliceends.c spliceends.h \
 bitpack64-access.c bitpack64-access.h \
 trdiagdef.h trdiag.c trdiag.h trdiagpool.c trdiagpool.h tr-extension-search.c tr-extension-search.h \
 univdiagdef.h univdiag.c univdiag.h univdiagpool.c univdiagpool.h sedgesort.c sedgesort.h \
 auxinfo.c auxinfo.h pathstore.c pathstore.h \
 altsplice.c altsplice.h \
 path.c path.h pathpair.c pathpair.h path-solve.c path-solve.h path-fusion.c path-fusion.h \
 path-trim.c path-trim.h path-eval.c path-eval.h pathpair-eval.c pathpair-eval.h \
 path-print-alignment.c path-print-alignment.h path-print-m8.c path-print-m8.h path-print-sam.c path-print-sam.h \
 path-learn.c path-learn.h \
 method.c method.h \
 doublelist.c doublelist.h bitvector.c bitvector.h \
 transcriptpool.c transcriptpool.h exon.c exon.h transcript.c transcript.h \
 trpathpool.c trpathpool.h trpath.c trpath.h trpath-solve.c trpath-solve.h trpath-convert.c trpath-convert.h \
 trpath-eval.c trpath-eval.h \
 transcriptome-search.c transcriptome-search.h \
 vectorpool.c vectorpool.h repair.c repair.h transcript-remap.c transcript-remap.h \
 transcript-velocity.c transcript-velocity.h \
 intersect-large.c intersect-large.h \
 intersect-small.c intersect-small.h intersect-uint2.c intersect-uint2.h \
 intersect-approx-uint4.c intersect-approx-uint4.h \
 intersect-approx-uint8.c intersect-approx-uint8.h \
 intersect-approx-indices-uint8.c intersect-approx-indices-uint8.h \
 intersect-wdups-indices.c intersect-wdups-indices.h \
 intersect-lower-large.c intersect-lower-large.h \
 intersect-higher-large.c intersect-higher-large.h \
 intersect-indices-small.c intersect-indices-small.h \
 intersect-indices2-large.c intersect-indices2-large.h \
 intersect-indices-uint8.c intersect-indices-uint8.h \
 intersect-concordance-uint8.c intersect-concordance-uint8.h \
 kmer-search.c kmer-search.h \
 extension-search.c extension-search.h \
 repetitive.c repetitive.h \
 uint8table.c uint8table.h \
 stage1hr.c stage1hr.h stage1hr-single.c stage1hr-single.h stage1hr-paired.c stage1hr-paired.h \
 uint8tableuint.c uint8tableuint.h \
 tableuint.c tableuint.h single-cell.c single-cell.h \
 request.c request.h resulthr.c resulthr.h outputtype.h output.c output.h \
 inbuffer.c inbuffer.h samheader.c samheader.h printbuffer.c printbuffer.h outbuffer.c outbuffer.h \
 uint8pairtable.c uint8pairtable.h splice-count.c splice-count.h support.c support.h \
 datadir.c datadir.h pass.h mode.h parserange.c parserange.h \
 getline.c getline.h getopt.c getopt1.c getopt.h gsnap.c


# Note: dist_ commands get read by bootstrap, and don't follow the flags

gsnapl_nosimd_CC = $(PTHREAD_CC)
gsnapl_nosimd_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGSNAP=1 -DLARGE_GENOMES=1
gsnapl_nosimd_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gsnapl_nosimd_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gsnapl_nosimd_SOURCES = $(GSNAPL_FILES)

gsnapl_arm_CC = $(PTHREAD_CC)
gsnapl_arm_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGSNAP=1 -DLARGE_GENOMES=1 -DHAVE_ARM=1 -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 $(SIMD_SSE4_2_CFLAGS)
gsnapl_arm_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gsnapl_arm_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gsnapl_arm_SOURCES = $(GSNAPL_FILES)

gsnapl_sse2_CC = $(PTHREAD_CC)
gsnapl_sse2_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGSNAP=1 -DLARGE_GENOMES=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 $(SIMD_SSE2_CFLAGS)
gsnapl_sse2_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gsnapl_sse2_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gsnapl_sse2_SOURCES = $(GSNAPL_FILES)

gsnapl_ssse3_CC = $(PTHREAD_CC)
gsnapl_ssse3_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGSNAP=1 -DLARGE_GENOMES=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 $(SIMD_SSSE3_CFLAGS)
gsnapl_ssse3_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gsnapl_ssse3_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gsnapl_ssse3_SOURCES = $(GSNAPL_FILES)

gsnapl_sse41_CC = $(PTHREAD_CC)
gsnapl_sse41_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGSNAP=1 -DLARGE_GENOMES=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 $(SIMD_SSE4_1_CFLAGS)
gsnapl_sse41_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gsnapl_sse41_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gsnapl_sse41_SOURCES = $(GSNAPL_FILES)

gsnapl_sse42_CC = $(PTHREAD_CC)
gsnapl_sse42_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGSNAP=1 -DLARGE_GENOMES=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 $(SIMD_SSE4_2_CFLAGS)
gsnapl_sse42_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gsnapl_sse42_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gsnapl_sse42_SOURCES = $(GSNAPL_FILES)

gsnapl_avx2_CC = $(PTHREAD_CC)
gsnapl_avx2_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGSNAP=1 -DLARGE_GENOMES=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 -DHAVE_AVX2=1 $(SIMD_AVX2_CFLAGS)
gsnapl_avx2_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gsnapl_avx2_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gsnapl_avx2_SOURCES = $(GSNAPL_FILES)

gsnapl_avx512_CC = $(PTHREAD_CC)
gsnapl_avx512_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGSNAP=1 -DLARGE_GENOMES=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 -DHAVE_AVX2=1 -DHAVE_AVX512=1 $(SIMD_AVX512_CFLAGS)
gsnapl_avx512_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gsnapl_avx512_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gsnapl_avx512_SOURCES = $(GSNAPL_FILES)

gsnapl_avx512bw_CC = $(PTHREAD_CC)
gsnapl_avx512bw_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGSNAP=1 -DLARGE_GENOMES=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 -DHAVE_AVX2=1 -DHAVE_AVX512=1 -DHAVE_AVX512BW=1 $(SIMD_AVX512_CFLAGS)
gsnapl_avx512bw_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gsnapl_avx512bw_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gsnapl_avx512bw_SOURCES = $(GSNAPL_FILES)


GEXACT_FILES = simd.h \
 bool.h types.h univcoord.h comp.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h uintlist.c uintlist.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 univinterval.c univinterval.h interval.c interval.h \
 stopwatch.c stopwatch.h semaphore.c semaphore.h access.c access.h \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iitdef.h iit-read.c iit-read.h \
 complement.h bzip2.c bzip2.h reader.c reader.h \
 genomicpos.c genomicpos.h compress.c compress.h \
 genome.c genome.h \
 popcount.c popcount.h \
 genomebits.c genomebits.h genomebits_consec.c genomebits_consec.h \
 bitpack64-read.c bitpack64-read.h bitpack64-readtwo.c bitpack64-readtwo.h \
 filesuffix.h \
 chrnum.c chrnum.h \
 samflags.h samprint-exact.c samprint-exact.h \
 fopen.c fopen.h shortread.c shortread.h \
 bitpack64-access.c bitpack64-access.h \
 bytecoding.c bytecoding.h sarray-read.c sarray-read.h \
 request.c request.h \
 inbuffer.c inbuffer.h samheader.c samheader.h printbuffer.c printbuffer.h outbuffer.c outbuffer.h \
 datadir.c datadir.h mode.h parserange.c parserange.h \
 getline.c getline.h getopt.c getopt1.c getopt.h gexact.c

gexact_nosimd_CC = $(PTHREAD_CC)
gexact_nosimd_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGEXACT=1
gexact_nosimd_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gexact_nosimd_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gexact_nosimd_SOURCES = $(GEXACT_FILES)

gexact_arm_CC = $(PTHREAD_CC)
gexact_arm_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGEXACT=1 -DHAVE_ARM=1 -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 $(SIMD_SSE4_2_CFLAGS)
gexact_arm_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gexact_arm_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gexact_arm_SOURCES = $(GEXACT_FILES)

gexact_sse2_CC = $(PTHREAD_CC)
gexact_sse2_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGEXACT=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 $(SIMD_SSE2_CFLAGS)
gexact_sse2_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gexact_sse2_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gexact_sse2_SOURCES = $(GEXACT_FILES)

gexact_ssse3_CC = $(PTHREAD_CC)
gexact_ssse3_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGEXACT=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 $(SIMD_SSSE3_CFLAGS)
gexact_ssse3_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gexact_ssse3_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gexact_ssse3_SOURCES = $(GEXACT_FILES)

gexact_sse41_CC = $(PTHREAD_CC)
gexact_sse41_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGEXACT=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 $(SIMD_SSE4_1_CFLAGS)
gexact_sse41_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gexact_sse41_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gexact_sse41_SOURCES = $(GEXACT_FILES)

gexact_sse42_CC = $(PTHREAD_CC)
gexact_sse42_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGEXACT=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 $(SIMD_SSE4_2_CFLAGS)
gexact_sse42_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gexact_sse42_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gexact_sse42_SOURCES = $(GEXACT_FILES)

gexact_avx2_CC = $(PTHREAD_CC)
gexact_avx2_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGEXACT=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 -DHAVE_AVX2=1 $(SIMD_AVX2_CFLAGS)
gexact_avx2_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gexact_avx2_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gexact_avx2_SOURCES = $(GEXACT_FILES)

gexact_avx512_CC = $(PTHREAD_CC)
gexact_avx512_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGEXACT=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 -DHAVE_AVX2=1 -DHAVE_AVX512=1 $(SIMD_AVX512_CFLAGS)
gexact_avx512_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gexact_avx512_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gexact_avx512_SOURCES = $(GEXACT_FILES)

gexact_avx512bw_CC = $(PTHREAD_CC)
gexact_avx512bw_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGEXACT=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 -DHAVE_AVX2=1 -DHAVE_AVX512=1 -DHAVE_AVX512BW=1 $(SIMD_AVX512_CFLAGS)
gexact_avx512bw_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gexact_avx512bw_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gexact_avx512bw_SOURCES = $(GEXACT_FILES)


GFILTER_FILES = simd.h \
 bool.h types.h univcoord.h comp.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h uintlist.c uintlist.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 stopwatch.c stopwatch.h semaphore.c semaphore.h access.c access.h \
 filestring.c filestring.h \
 complement.h bzip2.c bzip2.h reader.c reader.h oligo.c oligo.h \
 genomicpos.c genomicpos.h \
 popcount.c popcount.h \
 bitpack64-read.c bitpack64-read.h bitpack64-readtwo.c bitpack64-readtwo.h \
 filesuffix.h signature.c signature.h indexdbdef.h indexdb.c indexdb.h \
 intersectp-simd.h intersectp-simd.c \
 fopen.c fopen.h shortread.c shortread.h \
 request.c request.h \
 inbuffer.c inbuffer.h printbuffer.c printbuffer.h outbuffer.c outbuffer.h \
 datadir.c datadir.h \
 getline.c getline.h getopt.c getopt1.c getopt.h gfilter.c

gfilter_nosimd_CC = $(PTHREAD_CC)
gfilter_nosimd_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGFILTER=1
gfilter_nosimd_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gfilter_nosimd_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gfilter_nosimd_SOURCES = $(GFILTER_FILES)

gfilter_arm_CC = $(PTHREAD_CC)
gfilter_arm_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGFILTER=1 -DHAVE_ARM=1 -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 $(SIMD_SSE4_2_CFLAGS)
gfilter_arm_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gfilter_arm_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gfilter_arm_SOURCES = $(GFILTER_FILES)

gfilter_sse2_CC = $(PTHREAD_CC)
gfilter_sse2_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGFILTER=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 $(SIMD_SSE2_CFLAGS)
gfilter_sse2_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gfilter_sse2_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gfilter_sse2_SOURCES = $(GFILTER_FILES)

gfilter_ssse3_CC = $(PTHREAD_CC)
gfilter_ssse3_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGFILTER=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 $(SIMD_SSSE3_CFLAGS)
gfilter_ssse3_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gfilter_ssse3_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gfilter_ssse3_SOURCES = $(GFILTER_FILES)

gfilter_sse41_CC = $(PTHREAD_CC)
gfilter_sse41_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGFILTER=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 $(SIMD_SSE4_1_CFLAGS)
gfilter_sse41_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gfilter_sse41_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gfilter_sse41_SOURCES = $(GFILTER_FILES)

gfilter_sse42_CC = $(PTHREAD_CC)
gfilter_sse42_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGFILTER=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 $(SIMD_SSE4_2_CFLAGS)
gfilter_sse42_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gfilter_sse42_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gfilter_sse42_SOURCES = $(GFILTER_FILES)

gfilter_avx2_CC = $(PTHREAD_CC)
gfilter_avx2_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGFILTER=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 -DHAVE_AVX2=1 $(SIMD_AVX2_CFLAGS)
gfilter_avx2_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gfilter_avx2_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gfilter_avx2_SOURCES = $(GFILTER_FILES)

gfilter_avx512_CC = $(PTHREAD_CC)
gfilter_avx512_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGFILTER=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 -DHAVE_AVX2=1 -DHAVE_AVX512=1 $(SIMD_AVX512_CFLAGS)
gfilter_avx512_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gfilter_avx512_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gfilter_avx512_SOURCES = $(GFILTER_FILES)

gfilter_avx512bw_CC = $(PTHREAD_CC)
gfilter_avx512bw_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\" -DGFILTER=1 $(POPCNT_CFLAGS) -DHAVE_SSE2=1 -DHAVE_SSSE3=1 -DHAVE_SSE4_1=1 -DHAVE_SSE4_2=1 -DHAVE_AVX2=1 -DHAVE_AVX512=1 -DHAVE_AVX512BW=1 $(SIMD_AVX512_CFLAGS)
gfilter_avx512bw_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gfilter_avx512bw_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gfilter_avx512bw_SOURCES = $(GFILTER_FILES)


COMPARE2TRUTH_FILES = bool.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h uintlist.c uintlist.h list.c list.h \
 samread.c samread.h \
 getopt.c getopt1.c getopt.h compare2truth.c

compare2truth_CC = $(PTHREAD_CC)
compare2truth_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DTARGET=\"$(target)\" -DGMAPDB=\"$(GMAPDB)\"
compare2truth_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
compare2truth_LDADD = $(PTHREAD_LIBS)
dist_compare2truth_SOURCES = $(COMPARE2TRUTH_FILES)


# Build utility programs as non-SIMD programs
GMAPINDEX_FILES = simd.h \
 bool.h types.h univcoord.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 univinterval.c univinterval.h interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h semaphore.c semaphore.h access.c access.h \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iit-write-univ.c iit-write-univ.h \
 iitdef.h iit-read.c iit-read.h \
 md5.c md5.h complement.h bzip2.c bzip2.h fopen.c fopen.h sequence.c sequence.h genome.c genome.h \
 genomicpos.c genomicpos.h compress-write.c compress-write.h genome-write.c genome-write.h \
 compress.c compress.h popcount.c popcount.h \
 bitpack64-read.c bitpack64-read.h bitpack64-readtwo.c bitpack64-readtwo.h \
 bitpack64-access.c bitpack64-access.h bitpack64-incr.c bitpack64-incr.h bitpack64-write.c bitpack64-write.h \
 filesuffix.h signature.c signature.h indexdbdef.h indexdb.c indexdb.h indexdb-write.c indexdb-write.h \
 saca-k.c saca-k.h localdb-write.c localdb-write.h \
 table.c table.h tableuint.c tableuint.h tableuint8.c tableuint8.h tableint.c tableint.h \
 bytecoding.c bytecoding.h sarray-write.c sarray-write.h \
 chrom.c chrom.h segmentpos.c segmentpos.h \
 uint8list.c uint8list.h \
 parserange.c parserange.h getline.c getline.h \
 gmapindex.c

gmapindex_CC = $(PTHREAD_CC)
gmapindex_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DUTILITYP=1
gmapindex_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gmapindex_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gmapindex_SOURCES = $(GMAPINDEX_FILES)


INDEXDB_CAT_FILES = simd.h \
 bool.h types.h univcoord.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 genomicpos.c genomicpos.h \
 stopwatch.c stopwatch.h semaphore.c semaphore.h access.c access.h \
 univinterval.c univinterval.h interval.h interval.c \
 uintlist.c uintlist.h intlist.c intlist.h list.c list.h \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iitdef.h iit-read.h iit-read.c \
 popcount.c popcount.h \
 bitpack64-read.c bitpack64-read.h bitpack64-readtwo.c bitpack64-readtwo.h \
 bitpack64-access.c bitpack64-access.h bitpack64-incr.c bitpack64-incr.h bitpack64-write.c bitpack64-write.h \
 signature.c signature.h indexdbdef.h indexdb.c indexdb.h indexdb-write.c indexdb-write.h \
 complement.h compress-write.c compress-write.h bzip2.c bzip2.h fopen.c fopen.h \
 datadir.c datadir.h parserange.c parserange.h \
 getline.c getline.h getopt.c getopt1.c getopt.h indexdb-cat.c

indexdb_cat_CC = $(PTHREAD_CC)
indexdb_cat_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\" -DUTILITYP=1
indexdb_cat_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
indexdb_cat_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_indexdb_cat_SOURCES = $(INDEXDB_CAT_FILES)



SNPINDEX_FILES = simd.h \
 bool.h types.h univcoord.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 univinterval.c univinterval.h interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h semaphore.c semaphore.h access.c access.h \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iitdef.h iit-read.c iit-read.h \
 genomicpos.c genomicpos.h compress.c compress.h compress-write.c compress-write.h \
 popcount.c popcount.h \
 bitpack64-read.c bitpack64-read.h bitpack64-readtwo.c bitpack64-readtwo.h \
 bitpack64-access.c bitpack64-access.h bitpack64-incr.c bitpack64-incr.h bitpack64-write.c bitpack64-write.h \
 filesuffix.h signature.c signature.h indexdbdef.h indexdb.c indexdb.h indexdb-write.c indexdb-write.h \
 chrom.c chrom.h \
 complement.h md5.c md5.h bzip2.c bzip2.h fopen.c fopen.h sequence.c sequence.h genome.c genome.h \
 datadir.c datadir.h parserange.c parserange.h \
 getline.c getline.h getopt.c getopt1.c getopt.h snpindex.c

snpindex_CC = $(PTHREAD_CC)
snpindex_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\" -DUTILITYP=1
snpindex_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
snpindex_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_snpindex_SOURCES = $(SNPINDEX_FILES)


CMETINDEX_FILES = simd.h \
 bool.h types.h univcoord.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 genomicpos.c genomicpos.h \
 stopwatch.c stopwatch.h semaphore.c semaphore.h access.c access.h \
 univinterval.c univinterval.h interval.h interval.c \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iitdef.h iit-read.h iit-read.c \
 complement.h compress.c compress.h compress-write.c compress-write.h \
 popcount.c popcount.h \
 bitpack64-read.c bitpack64-read.h bitpack64-readtwo.c bitpack64-readtwo.h \
 bitpack64-access.c bitpack64-access.h bitpack64-incr.c bitpack64-incr.h bitpack64-write.c bitpack64-write.h \
 filesuffix.h signature.c signature.h indexdbdef.h indexdb.c indexdb.h indexdb-write.c indexdb-write.h \
 cmet.c cmet.h \
 complement.h md5.c md5.h bzip2.c bzip2.h fopen.c fopen.h sequence.c sequence.h genome.c genome.h \
 uintlist.c uintlist.h intlist.c intlist.h \
 uint8list.c uint8list.h \
 list.c list.h datadir.c datadir.h parserange.c parserange.h \
 getline.c getline.h getopt.c getopt1.c getopt.h cmetindex.c

cmetindex_CC = $(PTHREAD_CC)
cmetindex_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\" -DUTILITYP=1
cmetindex_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
cmetindex_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_cmetindex_SOURCES = $(CMETINDEX_FILES)


ATOIINDEX_FILES = simd.h \
 bool.h types.h univcoord.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 genomicpos.c genomicpos.h \
 stopwatch.c stopwatch.h semaphore.c semaphore.h access.c access.h \
 univinterval.c univinterval.h interval.h interval.c \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iitdef.h iit-read.h iit-read.c \
 complement.h compress.c compress.h compress-write.c compress-write.h \
 popcount.c popcount.h \
 bitpack64-read.c bitpack64-read.h bitpack64-readtwo.c bitpack64-readtwo.h \
 bitpack64-access.c bitpack64-access.h bitpack64-incr.c bitpack64-incr.h bitpack64-write.c bitpack64-write.h \
 filesuffix.h signature.c signature.h indexdbdef.h indexdb.c indexdb.h indexdb-write.c indexdb-write.h \
 atoi.c atoi.h \
 complement.h md5.c md5.h bzip2.c bzip2.h fopen.c fopen.h sequence.c sequence.h genome.c genome.h \
 uintlist.c uintlist.h intlist.c intlist.h \
 uint8list.c uint8list.h \
 list.c list.h datadir.c datadir.h parserange.c parserange.h \
 getline.c getline.h getopt.c getopt1.c getopt.h atoiindex.c

atoiindex_CC = $(PTHREAD_CC)
atoiindex_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\" -DUTILITYP=1
atoiindex_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
atoiindex_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_atoiindex_SOURCES = $(ATOIINDEX_FILES)


TRINDEX_FILES = simd.h \
 fopen.h bool.h types.h univcoord.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 genomicpos.c genomicpos.h \
 stopwatch.c stopwatch.h semaphore.c semaphore.h access.c access.h \
 interval.h interval.c \
 filestring.c filestring.h \
 iitdef.h iit-read.c iit-read.h \
 parserange.c parserange.h univinterval.c univinterval.h iit-read-univ.c iit-read-univ.h \
 popcount.c popcount.h \
 signature.c signature.h bitpack64-access.c bitpack64-access.h bitpack64-write.c bitpack64-write.h \
 uintlist.c uintlist.h intlist.c intlist.h list.c list.h \
 bitpack64-readtwo.c bitpack64-readtwo.h transcriptome.c transcriptome.h \
 datadir.c datadir.h \
 bzip2.c bzip2.h getline.c getline.h \
 getopt.c getopt1.c getopt.h trindex.c

trindex_CC = $(PTHREAD_CC)
trindex_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\" -DUTILITYP=1
trindex_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
trindex_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_trindex_SOURCES = $(TRINDEX_FILES)



GET_GENOME_FILES = bool.h types.h univcoord.h separator.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 univinterval.c univinterval.h interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h semaphore.c semaphore.h access.c access.h \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iitdef.h iit-read.c iit-read.h \
 md5.c md5.h complement.h bzip2.c bzip2.h fopen.c fopen.h sequence.c sequence.h \
 genome.c genome.h \
 genomicpos.c genomicpos.h chrom.c chrom.h \
 chrnum.c chrnum.h \
 datadir.c datadir.h parserange.c parserange.h \
 getline.c getline.h \
 getopt.c getopt1.c getopt.h get-genome.c

get_genome_CC = $(PTHREAD_CC)
get_genome_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\" -DUTILITYP=1
get_genome_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
get_genome_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_get_genome_SOURCES = $(GET_GENOME_FILES)


IIT_STORE_FILES = fopen.h bool.h types.h univcoord.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 univinterval.c univinterval.h interval.c interval.h \
 uintlist.c uintlist.h \
 stopwatch.c stopwatch.h semaphore.c semaphore.h access.c access.h \
 doublelist.c doublelist.h \
 iit-write-univ.c iit-write-univ.h iitdef.h iit-write.c iit-write.h \
 tableint.c tableint.h table.c table.h chrom.c chrom.h \
 bzip2.c bzip2.h getline.c getline.h \
 getopt.c getopt1.c getopt.h iit_store.c

iit_store_CC = $(PTHREAD_CC)
iit_store_CFLAGS = $(AM_CFLAGS)  $(PTHREAD_CFLAGS) -DUTILITYP=1
iit_store_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
iit_store_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_iit_store_SOURCES = $(IIT_STORE_FILES)


IIT_GET_FILES = fopen.h bool.h types.h univcoord.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h list.c list.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 univinterval.c univinterval.h interval.c interval.h \
 uintlist.c uintlist.h \
 stopwatch.c stopwatch.h semaphore.c semaphore.h access.c access.h \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iitdef.h iit-read.c iit-read.h \
 complement.h parserange.c parserange.h \
 bzip2.c bzip2.h getline.c getline.h \
 getopt.c getopt1.c getopt.h iit_get.c

iit_get_CC = $(PTHREAD_CC)
iit_get_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DUTILITYP=1
iit_get_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
iit_get_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_iit_get_SOURCES = $(IIT_GET_FILES)

IIT_DUMP_FILES = fopen.h bool.h types.h univcoord.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 intlist.c intlist.h list.c list.h \
 univinterval.c univinterval.h interval.c interval.h uintlist.c uintlist.h \
 stopwatch.c stopwatch.h semaphore.c semaphore.h access.c access.h \
 filestring.c filestring.h \
 iit-read-univ.c iit-read-univ.h iitdef.h iit-read.c iit-read.h parserange.c parserange.h \
 getopt.c getopt1.c getopt.h iit_dump.c

iit_dump_CC = $(PTHREAD_CC)
iit_dump_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DUTILITYP=1
iit_dump_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
iit_dump_LDADD = $(PTHREAD_LIBS)
dist_iit_dump_SOURCES = $(IIT_DUMP_FILES)


IIT_PILEUP_FILES = fopen.h bool.h types.h univcoord.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 semaphore.c semaphore.h stopwatch.c stopwatch.h access.c access.h \
 list.c list.h filestring.c filestring.h \
 genomicpos.h interval.c interval.h \
 intlist.c intlist.h uintlist.c uintlist.h iitdef.h iit-read.c iit-read.h \
 univinterval.c univinterval.h iit-read-univ.c iit-read-univ.h parserange.c parserange.h \
 getopt.c getopt1.c getopt.h iit_pileup.c

iit_pileup_CC = $(PTHREAD_CC)
iit_pileup_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DUTILITYP=1
iit_pileup_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
iit_pileup_LDADD = $(PTHREAD_LIBS)
dist_iit_pileup_SOURCES = $(IIT_PILEUP_FILES)



SAM_SORT_FILES = bool.h types.h univcoord.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 littleendian.c littleendian.h bigendian.c bigendian.h \
 genomicpos.h chrnum.h \
 intlist.c intlist.h list.c list.h \
 stopwatch.c stopwatch.h semaphore.c semaphore.h access.c access.h \
 univinterval.c univinterval.h interval.c interval.h \
 uintlist.c uintlist.h \
 filestring.c filestring.h \
 parserange.c parserange.h iit-read-univ.c iit-read-univ.h iitdef.h iit-read.c iit-read.h \
 samflags.h samheader.c samheader.h samread.c samread.h \
 bzip2.c bzip2.h getline.c getline.h datadir.c datadir.h \
 getopt.c getopt1.c getopt.h sam_sort.c

sam_sort_CC = $(PTHREAD_CC)
sam_sort_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\" -DUTILITYP=1
sam_sort_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
sam_sort_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_sam_sort_SOURCES = $(SAM_SORT_FILES)


UMI_CORRECT_FILES = bool.h types.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h uintlist.c uintlist.h filestring.c filestring.h interval.c interval.h \
 iitdef.h iit-read.c iit-read.h \
 samread.c samread.h \
 table.c table.h tableint.c tableint.h \
 bzip2.c bzip2.h getline.c getline.h \
 stopwatch.c stopwatch.h semaphore.c semaphore.h access.c access.h \
 list.c list.h datadir.c datadir.h \
 getopt.c getopt1.c getopt.h umi-correct.c

umi_correct_CC = $(PTHREAD_CC)
umi_correct_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\" -DUTILITYP=1
umi_correct_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
umi_correct_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_umi_correct_SOURCES = $(UMI_CORRECT_FILES)


VELOCITY_COUNTS_FILES = bool.h \
 except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h uintlist.c uintlist.h filestring.c filestring.h interval.c interval.h \
 iitdef.h iit-read.c iit-read.h \
 samread.c samread.h \
 table.c table.h tableint.c tableint.h \
 popcount.c popcount.h \
 bzip2.c bzip2.h getline.c getline.h \
 stopwatch.c stopwatch.h semaphore.c semaphore.h access.c access.h \
 list.c list.h datadir.c datadir.h \
 getopt.c getopt1.c getopt.h velocity-counts.c

velocity_counts_CC = $(PTHREAD_CC)
velocity_counts_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\" -DUTILITYP=1 $(POPCNT_CFLAGS)
velocity_counts_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
velocity_counts_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_velocity_counts_SOURCES = $(VELOCITY_COUNTS_FILES)


GCOUNT_BULK_FILES = simd.h \
 bool.h except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h uintlist.c uintlist.h filestring.c filestring.h interval.c interval.h \
 iitdef.h iit-read.c iit-read.h \
 parserange.c parserange.h univinterval.c univinterval.h iit-read-univ.c iit-read-univ.h chrnum.c chrnum.h \
 genomicpos.c genomicpos.h bitpack64-readtwo.c bitpack64-readtwo.h transcriptome.c transcriptome.h \
 samread.c samread.h \
 table.c table.h tableint.c tableint.h tableuint.c tableuint.h \
 uinttableuint.c uinttableuint.h uinttable.c uinttable.h \
 doublelist.c doublelist.h \
 sedgesort.c sedgesort.h equiv-class-bulk.c equiv-class-bulk.h \
 orderstat.c orderstat.h \
 popcount.c popcount.h bitvector.c bitvector.h \
 bzip2.c bzip2.h getline.c getline.h \
 stopwatch.c stopwatch.h semaphore.c semaphore.h access.c access.h \
 list.c list.h datadir.c datadir.h \
 getopt.c getopt1.c getopt.h gcount-bulk.c

gcount_bulk_CC = $(PTHREAD_CC)
gcount_bulk_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\" -DUTILITYP=1 $(POPCNT_CFLAGS)
gcount_bulk_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gcount_bulk_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gcount_bulk_SOURCES = $(GCOUNT_BULK_FILES)


GCOUNT_SC_FILES = simd.h \
 bool.h except.c except.h assert.c assert.h mem.c mem.h \
 intlist.c intlist.h uintlist.c uintlist.h filestring.c filestring.h interval.c interval.h \
 iitdef.h iit-read.c iit-read.h \
 parserange.c parserange.h univinterval.c univinterval.h iit-read-univ.c iit-read-univ.h chrnum.c chrnum.h \
 genomicpos.c genomicpos.h bitpack64-readtwo.c bitpack64-readtwo.h transcriptome.c transcriptome.h \
 samread.c samread.h \
 table.c table.h tableint.c tableint.h tableuint.c tableuint.h \
 uinttableuint.c uinttableuint.h \
 doublelist.c doublelist.h \
 transcript-splice.c transcript-splice.h \
 sedgesort.c sedgesort.h equiv-class-sc.c equiv-class-sc.h \
 orderstat.c orderstat.h \
 popcount.c popcount.h bitvector.c bitvector.h \
 bzip2.c bzip2.h getline.c getline.h \
 stopwatch.c stopwatch.h semaphore.c semaphore.h access.c access.h \
 list.c list.h datadir.c datadir.h \
 getopt.c getopt1.c getopt.h gcount-sc.c

gcount_sc_CC = $(PTHREAD_CC)
gcount_sc_CFLAGS = $(AM_CFLAGS) $(PTHREAD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\" -DUTILITYP=1 $(POPCNT_CFLAGS)
gcount_sc_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG)
gcount_sc_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
dist_gcount_sc_SOURCES = $(GCOUNT_SC_FILES)


#SPLICING_SCORE_FILES = bool.h types.h univcoord.h \
# except.c except.h assert.c assert.h mem.c mem.h \
# intlist.c intlist.h list.c list.h \
# univinterval.c univinterval.h interval.c interval.h \
# uintlist.c uintlist.h \
# chrom.c chrom.h stopwatch.c stopwatch.h semaphore.c semaphore.h access.c access.h \
# iit-read-univ.c iit-read-univ.h iitdef.h iit-read.c iit-read.h \
# filestring.c filestring.h \
# md5.c md5.h complement.h bzip2.c bzip2.h fopen.c fopen.h sequence.c sequence.h \
# genome.c genome.h \
# genomicpos.c genomicpos.h \
# chrnum.c chrnum.h \
# maxent.c maxent.h maxent_hr.c maxent_hr.h \
# branchpoint.c branchpoint.h \
# parserange.c parserange.h datadir.c datadir.h getopt.c getopt1.c getopt.h splicing-score.c
#
#splicing_score_CFLAGS = $(AM_CFLAGS)  $(PTHREAD_CFLAGS) -DGMAPDB=\"$(GMAPDB)\"
#splicing_score_LDFLAGS = $(AM_LDFLAGS) $(STATIC_LDFLAG) $(PTHREAD_CFLAGS)
#splicing_score_LDADD = $(PTHREAD_LIBS) $(ZLIB_LIBS) $(BZLIB_LIBS)
#dist_splicing_score_SOURCES =
#nodist_splicing_score_SOURCES = $(SPLICING_SCORE_FILES)
