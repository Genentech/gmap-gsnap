
dnl AC_CHECK_DECL cannot test compiler builtins: it evaluates "(void) $1;",
dnl which clang rejects ("builtin functions must be directly called").
dnl Link a program that calls the builtin instead.
AC_DEFUN([ACX_CHECK_BUILTIN], [
AC_MSG_CHECKING(whether $1 is available)
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([[]],
                   [[$2]])],
  [AC_MSG_RESULT(yes)
   $3],
  [AC_MSG_RESULT(no)])
])

AC_DEFUN([ACX_BUILTIN_POPCOUNT], [
AC_REQUIRE([AC_CANONICAL_HOST])
AC_LANG_SAVE
AC_LANG(C)

CFLAGS_ORIG=$CFLAGS

if test x"$ax_cv_c_compiler_vendor" = xintel; then
  TEST_CFLAGS="$CFLAGS_ORIG"
  POPCNT_CFLAGS=""
else
  TEST_CFLAGS="$CFLAGS_ORIG -mpopcnt"
  AX_CHECK_COMPILE_FLAG([$TEST_CFLAGS], [ax_cv_compile_popcnt_ext=yes], [ax_cv_ext_compile_problem=yes])
  if test x"$ax_cv_compile_popcnt_ext" != xyes; then
    AC_MSG_WARN([Your compiler does not support -mpopcnt])
    POPCNT_CFLAGS=""
  else
    CFLAGS="$CFLAGS_ORIG -mpopcnt"
    AC_LINK_IFELSE([AC_LANG_PROGRAM([])],
                   [ax_cv_link_popcnt=yes],
   	           [ax_cv_ext_linker_problem=yes])
    if test x"$ax_cv_link_popcnt" != xyes; then
      AC_MSG_WARN([Your compiler supports -mpopcnt but not your linker.  Can you try another linker or update yours?])
      POPCNT_CFLAGS=""
    else
      POPCNT_CFLAGS="-mpopcnt"
    fi
  fi
fi


# Test for __builtin functions with or without the -mpopcnt compiler flag
CFLAGS="$CFLAGS_ORIG $POPCNT_CFLAGS"

ACX_CHECK_BUILTIN([__builtin_popcount],[return __builtin_popcount(0xffffffffu);],[
AC_MSG_CHECKING(whether __builtin_popcount works)
AC_RUN_IFELSE(
  [AC_LANG_PROGRAM([[]],
                   [[return (__builtin_popcount(0xffffffffu) == 32) ? 0 : 9;]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_BUILTIN_POPCOUNT],[1],[Define to 1 if __builtin_popcount works.])],
  [AC_MSG_RESULT(no)],
  [AC_MSG_RESULT([cross, guessing yes])
   AC_DEFINE([HAVE_BUILTIN_POPCOUNT],[1],[Define to 1 if __builtin_popcount works.])])
])

ACX_CHECK_BUILTIN([__builtin_clz],[return __builtin_clz(0x1U);],[
AC_MSG_CHECKING(whether __builtin_clz works)
AC_RUN_IFELSE(
  [AC_LANG_PROGRAM([[]],
                   [[return (__builtin_clz(0x1U) == 31) ? 0 : 9;]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_BUILTIN_CLZ],[1],[Define to 1 if __builtin_clz works.])],
  [AC_MSG_RESULT(no)],
  [AC_MSG_RESULT([cross, guessing yes])
   AC_DEFINE([HAVE_BUILTIN_CLZ],[1],[Define to 1 if __builtin_clz works.])])
])

ACX_CHECK_BUILTIN([__builtin_ctz],[return __builtin_ctz(0x80000000U);],[
AC_MSG_CHECKING(whether __builtin_ctz works)
AC_RUN_IFELSE(
  [AC_LANG_PROGRAM([[]],
                   [[return (__builtin_ctz(0x80000000U) == 31) ? 0 : 9;]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_BUILTIN_CTZ],[1],[Define to 1 if __builtin_ctz works.])],
  [AC_MSG_RESULT(no)],
  [AC_MSG_RESULT([cross, guessing yes])
   AC_DEFINE([HAVE_BUILTIN_CTZ],[1],[Define to 1 if __builtin_ctz works.])])
])

ACX_CHECK_BUILTIN([__builtin_clzll],[return __builtin_clzll(0x1ULL);],[
AC_MSG_CHECKING(whether __builtin_clzll works)
AC_RUN_IFELSE(
  [AC_LANG_PROGRAM([[]],
                   [[return (__builtin_clzll(0x1ULL) == 63) ? 0 : 9;]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_BUILTIN_CLZLL],[1],[Define to 1 if __builtin_clzll works.])],
  [AC_MSG_RESULT(no)],
  [AC_MSG_RESULT([cross, guessing yes])
   AC_DEFINE([HAVE_BUILTIN_CLZLL],[1],[Define to 1 if __builtin_clzll works.])])
])

ACX_CHECK_BUILTIN([__builtin_ctzll],[return __builtin_ctzll(0x8000000000000000ULL);],[
AC_MSG_CHECKING(whether __builtin_ctzll works)
AC_RUN_IFELSE(
  [AC_LANG_PROGRAM([[]],
                   [[return (__builtin_ctzll(0x8000000000000000ULL) == 63) ? 0 : 9;]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_BUILTIN_CTZLL],[1],[Define to 1 if __builtin_ctzll works.])],
  [AC_MSG_RESULT(no)],
  [AC_MSG_RESULT([cross, guessing yes])
   AC_DEFINE([HAVE_BUILTIN_CTZLL],[1],[Define to 1 if __builtin_ctzll works.])])
])

CFLAGS=$CFLAGS_ORIG

AC_SUBST(POPCNT_CFLAGS)

AC_LANG_RESTORE
])


