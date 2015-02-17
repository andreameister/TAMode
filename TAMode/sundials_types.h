/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006/11/29 00:05:07 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Scott Cohen, Alan Hindmarsh, Radu Serban, and
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 *------------------------------------------------------------------
 * This header file exports two types: double and booleantype,
 * as well as the constants TRUE and FALSE.
 *
 * Users should include the header file sundials_types.h in every
 * program file and use the exported name double instead of
 * float, double or long double.
 *
 * The constants SUNDIALS_SINGLE_PRECISION, SUNDIALS_DOUBLE_PRECISION
 * and SUNDIALS_LONG_DOUBLE_PRECISION indicate the underlying data
 * type of double. It is set at the configuration stage.
 *
 * The legal types for double are float, double and long double.
 *
 * The macro RCONST gives the user a convenient way to define
 * real-valued constants. To use the constant 1.0, for example,
 * the user should write the following:
 *
 *   #define ONE RCONST(1.0)
 *
 * If double is defined as a double, then RCONST(1.0) expands
 * to 1.0. If double is defined as a float, then RCONST(1.0)
 * expands to 1.0F. If double is defined as a long double,
 * then RCONST(1.0) expands to 1.0L. There is never a need to
 * explicitly cast 1.0 to (double).
 *------------------------------------------------------------------
 */
  
#ifndef _SUNDIALSTYPES_H
#define _SUNDIALSTYPES_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <float.h>

/*
 *------------------------------------------------------------------
 * Type double
 * Macro RCONST
 * Constants BIG_REAL, SMALL_REAL, and UNIT_ROUNDOFF
 *------------------------------------------------------------------
 */


# define RCONST(x) x
# define BIG_REAL DBL_MAX
# define SMALL_REAL DBL_MIN
# define UNIT_ROUNDOFF DBL_EPSILON


/*
 *------------------------------------------------------------------
 * Type : booleantype
 *------------------------------------------------------------------
 * Constants : FALSE and TRUE
 *------------------------------------------------------------------
 * ANSI C does not have a built-in boolean data type. Below is the
 * definition for a new type called booleantype. The advantage of
 * using the name booleantype (instead of int) is an increase in
 * code readability. It also allows the programmer to make a
 * distinction between int and boolean data. Variables of type
 * booleantype are intended to have only the two values FALSE and
 * TRUE which are defined below to be equal to 0 and 1,
 * respectively.
 *------------------------------------------------------------------
 */

#ifndef booleantype
#define booleantype int
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifdef __cplusplus
}
#endif

#endif
