// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// assenteAcesso
List assenteAcesso(List B, List C, List D, int g, double r, double eps);
RcppExport SEXP _menina_assenteAcesso(SEXP BSEXP, SEXP CSEXP, SEXP DSEXP, SEXP gSEXP, SEXP rSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type B(BSEXP);
    Rcpp::traits::input_parameter< List >::type C(CSEXP);
    Rcpp::traits::input_parameter< List >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type g(gSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(assenteAcesso(B, C, D, g, r, eps));
    return rcpp_result_gen;
END_RCPP
}
// avalieAcesso
List avalieAcesso(List A, List D, int g, double r, double eps);
RcppExport SEXP _menina_avalieAcesso(SEXP ASEXP, SEXP DSEXP, SEXP gSEXP, SEXP rSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type A(ASEXP);
    Rcpp::traits::input_parameter< List >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type g(gSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(avalieAcesso(A, D, g, r, eps));
    return rcpp_result_gen;
END_RCPP
}
// busqueBase
List busqueBase(List C, List D, int g, double r, double eps);
RcppExport SEXP _menina_busqueBase(SEXP CSEXP, SEXP DSEXP, SEXP gSEXP, SEXP rSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type C(CSEXP);
    Rcpp::traits::input_parameter< List >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type g(gSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(busqueBase(C, D, g, r, eps));
    return rcpp_result_gen;
END_RCPP
}
// computeComponente
List computeComponente(List D, int g, double r, double eps);
RcppExport SEXP _menina_computeComponente(SEXP DSEXP, SEXP gSEXP, SEXP rSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type g(gSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(computeComponente(D, g, r, eps));
    return rcpp_result_gen;
END_RCPP
}
// definaDispositivo
List definaDispositivo(int g, int mu);
RcppExport SEXP _menina_definaDispositivo(SEXP gSEXP, SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type g(gSEXP);
    Rcpp::traits::input_parameter< int >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(definaDispositivo(g, mu));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_menina_assenteAcesso", (DL_FUNC) &_menina_assenteAcesso, 6},
    {"_menina_avalieAcesso", (DL_FUNC) &_menina_avalieAcesso, 5},
    {"_menina_busqueBase", (DL_FUNC) &_menina_busqueBase, 5},
    {"_menina_computeComponente", (DL_FUNC) &_menina_computeComponente, 4},
    {"_menina_definaDispositivo", (DL_FUNC) &_menina_definaDispositivo, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_menina(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
