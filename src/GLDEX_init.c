#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void check_gld(void *, void *, void *, void *, void *, void *);
extern void dgl(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void mult_check_gld(void *, void *, void *, void *, void *, void *, void *);
extern void optim_fun3(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void optim_fun3_v(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void pgl(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void q_fmkl_gld_minmax_check(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void q_rs_gld_minmax_check(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Fortran calls */
extern void F77_NAME(halton)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(sobol)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"check_gld",               (DL_FUNC) &check_gld,                6},
    {"dgl",                     (DL_FUNC) &dgl,                     11},
    {"mult_check_gld",          (DL_FUNC) &mult_check_gld,           7},
    {"optim_fun3",              (DL_FUNC) &optim_fun3,              11},
    {"optim_fun3_v",            (DL_FUNC) &optim_fun3_v,            12},
    {"pgl",                     (DL_FUNC) &pgl,                     11},
    {"q_fmkl_gld_minmax_check", (DL_FUNC) &q_fmkl_gld_minmax_check, 11},
    {"q_rs_gld_minmax_check",   (DL_FUNC) &q_rs_gld_minmax_check,   11},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"halton", (DL_FUNC) &F77_NAME(halton),  7},
    {"sobol",  (DL_FUNC) &F77_NAME(sobol),  11},
    {NULL, NULL, 0}
};

void R_init_GLDEX(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
