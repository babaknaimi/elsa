#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP C_elsa_local_krige(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_elsa_local_krige_v2(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,
                                  SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,
                                  SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern SEXP C_categorize(SEXP, SEXP);
extern SEXP dist(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dnn(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_elsa(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP elsa_cell(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP elsa_cell_test(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP elsa_test(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP elsa_vector(SEXP, SEXP, SEXP);
extern SEXP elsa_vector_test(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP elsac(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP elsac_cell(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP elsac_cell_test(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP elsac_test(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP elsac_vector(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP elsac_vector_test(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_geary(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP geary_vector(SEXP, SEXP);
extern SEXP GG(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GG_vector(SEXP, SEXP);
extern SEXP localgeary(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP localgeary_vector(SEXP, SEXP);
extern SEXP localmoran(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP localmoran_vector(SEXP, SEXP);
extern SEXP Melsa(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_moran(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP moran_vector(SEXP, SEXP);
extern SEXP poly_loop2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP semivar(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP semivar_vector(SEXP, SEXP);
extern SEXP v_elsa(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP v_elsa_cell(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP v_elsa_cell_Ea(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP v_elsa_cell_Ec(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP v_elsa_vector(SEXP, SEXP, SEXP);
extern SEXP v_elsa_vector_Ea(SEXP, SEXP, SEXP);
extern SEXP v_elsa_vector_Ec(SEXP, SEXP, SEXP);
extern SEXP v_elsac(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP v_elsac_cell(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP v_elsac_cell_Ea(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP v_elsac_cell_Ec(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP v_elsac_vector(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP v_elsac_vector_Ea(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP v_elsac_vector_Ec(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"C_elsa_local_krige", (DL_FUNC) &C_elsa_local_krige, 17},
    {"C_elsa_local_krige_v2", (DL_FUNC) &C_elsa_local_krige_v2, 24},
    {"C_categorize",         (DL_FUNC) &C_categorize,      2},
    {"dist",               (DL_FUNC) &dist,                5},
    {"dnn",                (DL_FUNC) &dnn,                 5},
    {"C_elsa",               (DL_FUNC) &C_elsa,            6},
    {"elsa_cell",          (DL_FUNC) &elsa_cell,           7},
    {"elsa_cell_test",     (DL_FUNC) &elsa_cell_test,     10},
    {"elsa_test",          (DL_FUNC) &elsa_test,           9},
    {"elsa_vector",        (DL_FUNC) &elsa_vector,         3},
    {"elsa_vector_test",   (DL_FUNC) &elsa_vector_test,    6},
    {"elsac",              (DL_FUNC) &elsac,               8},
    {"elsac_cell",         (DL_FUNC) &elsac_cell,          9},
    {"elsac_cell_test",    (DL_FUNC) &elsac_cell_test,    12},
    {"elsac_test",         (DL_FUNC) &elsac_test,         11},
    {"elsac_vector",       (DL_FUNC) &elsac_vector,        5},
    {"elsac_vector_test",  (DL_FUNC) &elsac_vector_test,   8},
    {"C_geary",              (DL_FUNC) &C_geary,           5},
    {"geary_vector",       (DL_FUNC) &geary_vector,        2},
    {"GG",                 (DL_FUNC) &GG,                  5},
    {"GG_vector",          (DL_FUNC) &GG_vector,           2},
    {"localgeary",         (DL_FUNC) &localgeary,          5},
    {"localgeary_vector",  (DL_FUNC) &localgeary_vector,   2},
    {"localmoran",         (DL_FUNC) &localmoran,          5},
    {"localmoran_vector",  (DL_FUNC) &localmoran_vector,   2},
    {"Melsa",              (DL_FUNC) &Melsa,               6},
    {"C_moran",              (DL_FUNC) &C_moran,           5},
    {"moran_vector",       (DL_FUNC) &moran_vector,        2},
    {"poly_loop2",         (DL_FUNC) &poly_loop2,          8},
    {"semivar",            (DL_FUNC) &semivar,             5},
    {"semivar_vector",     (DL_FUNC) &semivar_vector,      2},
    {"v_elsa",             (DL_FUNC) &v_elsa,              6},
    {"v_elsa_cell",        (DL_FUNC) &v_elsa_cell,         7},
    {"v_elsa_cell_Ea",     (DL_FUNC) &v_elsa_cell_Ea,      7},
    {"v_elsa_cell_Ec",     (DL_FUNC) &v_elsa_cell_Ec,      7},
    {"v_elsa_vector",      (DL_FUNC) &v_elsa_vector,       3},
    {"v_elsa_vector_Ea",   (DL_FUNC) &v_elsa_vector_Ea,    3},
    {"v_elsa_vector_Ec",   (DL_FUNC) &v_elsa_vector_Ec,    3},
    {"v_elsac",            (DL_FUNC) &v_elsac,             8},
    {"v_elsac_cell",       (DL_FUNC) &v_elsac_cell,        9},
    {"v_elsac_cell_Ea",    (DL_FUNC) &v_elsac_cell_Ea,     9},
    {"v_elsac_cell_Ec",    (DL_FUNC) &v_elsac_cell_Ec,     9},
    {"v_elsac_vector",     (DL_FUNC) &v_elsac_vector,      5},
    {"v_elsac_vector_Ea",  (DL_FUNC) &v_elsac_vector_Ea,   5},
    {"v_elsac_vector_Ec",  (DL_FUNC) &v_elsac_vector_Ec,   5},
    {NULL, NULL, 0}
};

void R_init_elsa(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
