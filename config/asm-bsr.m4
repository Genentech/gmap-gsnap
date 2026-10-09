
AC_DEFUN([ACX_ASM_BSR], [
AC_LANG_SAVE
AC_LANG(C)

AC_MSG_CHECKING(for bsr instruction in assembly)
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([[#include <stdlib.h>]],
                   [[int msb; unsigned int x = rand(); asm volatile("bsr %1,%0" : "=r"(msb) : "r"(x));]])],
  [AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_ASM_BSR],[1],[Define to 1 if bsr command is available in assembly])],
  [AC_MSG_RESULT(no)])

AC_LANG_RESTORE
])

