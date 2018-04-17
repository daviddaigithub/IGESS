#define ARMA_DONT_USE_WRAPPER
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
#include <RcppArmadillo.h>
#include <R.h>
#include <Rinternals.h>
#include "Rcpp_aux.hpp"
#include "IGESS_aux.hpp"
#include "IGESS.hpp"
#include "readPlink.hpp"

using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
RcppExport SEXP IGESS(SEXP Xs, arma::vec& y, SEXP  SS = R_NilValue,  SEXP opts = R_NilValue,
                                std::string logfile="screen", double lbPval = 1e-12, bool verbose = false){
  streambuf* coutBuf = NULL;
  ofstream of(logfile.c_str(),std::ios::app);
  if(logfile != "screen"){
    coutBuf = cout.rdbuf();
    streambuf* fileBuf = of.rdbuf();
    cout.rdbuf(fileBuf);
  }
  Options* lp_opt = read_opts(opts);
  mat summary;
  int N, P, type;
  void* ptr;
  List ret;
  preprocess_summary(Xs, SS, lbPval, type, N, P, summary, ptr);
  IGESSfit* fit = NULL;
  if(type == type_int){
    int* data = static_cast<Mat<int> *>(ptr) ->memptr();
    fit = iGess(data, y, P, &summary,lp_opt, type, false, verbose);
  }else{
    double* data = static_cast<Mat<double> *>(ptr) ->memptr();
    fit = iGess(data , y, P, &summary,lp_opt, type, false, verbose);

  }
  if(logfile != "screen"){
    cout.rdbuf(coutBuf);
    cout << "Write Personal Information over..." << endl;
  }
  of.flush();
  of.close();
  ret = wrap_fit(fit);
  return ret;
}


// [[Rcpp::export]]
RcppExport SEXP IGESSCV(SEXP Xs, arma::vec& y, SEXP  SS = R_NilValue,  SEXP opts = R_NilValue, std::string logfile="screen", double lbPval = 1e-12, Rcpp::String measure = "mse"){
  streambuf* coutBuf = NULL;
  ofstream of(logfile.c_str(),std::ios::app);
  if(logfile != "screen"){
    coutBuf = cout.rdbuf();
    streambuf* fileBuf = of.rdbuf();
    cout.rdbuf(fileBuf);
  }
  Options* lp_opt = read_opts(opts);
  mat summary;
  int N, P, type;
  void* ptr;
  List ret;
  preprocess_summary(Xs, SS, lbPval, type, N, P, summary, ptr);
  double predict = 0;
  if(type == type_int){
    // cout << static_cast<Mat<int> *>(ptr) -> n_rows << " " << static_cast<Mat<int> *>(ptr) -> n_cols << endl;
    int* data = static_cast<Mat<int> *>(ptr) ->memptr();
    predict = iGessCV(data, y,P,  &summary, lp_opt, type, false, measure);
  }else{
    // cout << static_cast<Mat<double> *>(ptr) -> n_rows << " " << static_cast<Mat<double> *>(ptr) -> n_cols << endl;
    double* data = static_cast<Mat<double> *>(ptr) ->memptr();
    predict = iGessCV(data, y,P,  &summary, lp_opt, type, false, measure);
  }
  if(logfile != "screen"){
    cout.rdbuf(coutBuf);
    cout << "Write Personal Information over..." << endl;
  }
  of.flush();
  of.close();
  ret[measure] = wrap(predict);
  return ret;
}


// [[Rcpp::export]]
RcppExport SEXP IGESS_Plink(Rcpp::String genoplinkfile, SEXP  SS = R_NilValue, SEXP opts = R_NilValue, std::string logfile="screen", double lbPval = 1e-12){
  streambuf* coutBuf = NULL;
  ofstream of(logfile.c_str(),std::ios::app);
  if(logfile != "screen"){
    coutBuf = cout.rdbuf();
    streambuf* fileBuf = of.rdbuf();
    cout.rdbuf(fileBuf);
  }
  Options* lp_opt = read_opts(opts);
  mat summary;
  int N, P;
  int type = type_int;
  List ret;
  GenoInfo geno(genoplinkfile);
  summary = combine_summary_X(geno.X, geno.chroms.snps, SS, lbPval);
  N = geno.X.n_rows;
  P = geno.X.n_cols;
  IGESSfit* fit = iGess(geno.X.memptr(), geno.y, P, &summary,lp_opt, type, false);
  ret = wrap_fit(fit);
  if(logfile != "screen"){
    cout.rdbuf(coutBuf);
    cout << "Write Personal Information over..." << endl;
  }
  of.flush();
  of.close();
  ret["X"]  = geno.X;
  ret["y"]  = geno.y;
  return ret;
}


// [[Rcpp::export]]
RcppExport SEXP IGESSCV_Plink(Rcpp::String genoplinkfile, SEXP  SS = R_NilValue,  SEXP opts = R_NilValue, std::string logfile="screen", double lbPval = 1e-12,  Rcpp::String measure = "mse"){
  streambuf* coutBuf = NULL;
  ofstream of(logfile.c_str(),std::ios::app);
  if(logfile != "screen"){
    coutBuf = cout.rdbuf();
    streambuf* fileBuf = of.rdbuf();
    cout.rdbuf(fileBuf);
  }
  Options* lp_opt = read_opts(opts);
  mat summary;
  int N, P;
  int type = type_int;
  List ret;
  GenoInfo geno(genoplinkfile);
  summary = combine_summary_X(geno.X, geno.chroms.snps, SS, lbPval);
  N = geno.X.n_rows;
  P = geno.X.n_cols;
  double predict = iGessCV(geno.X.memptr(), geno.y, P,  &summary, lp_opt, type, false, measure);
  ret[measure] = wrap(predict);
  if(logfile != "screen"){
    cout.rdbuf(coutBuf);
    cout << "Write Personal Information over..." << endl;
  }
  of.flush();
  of.close();
  return ret;
}



// [[Rcpp::export]]
RcppExport SEXP IGESS_Predict(SEXP fit_,  arma::mat& X){
  vec ypred(X.n_rows);
  IGESSfit* fit = NULL;
  if(!Rf_isNull(fit_)){
    Rcpp::List fitList(fit_);
    fit = new IGESSfit(fitList["gammas"], fitList["mu"], fitList["cov"]);
    ypred = fit -> predict(&X);
  }else{
    cout << "Invalid input of iGESS fit!" << endl;
    return Rcpp::wrap(ypred);
  }
  return Rcpp::wrap(ypred);
}


