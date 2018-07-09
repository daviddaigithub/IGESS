//
//  IGESS.hpp
//  IGESSArma
//
//  Created by DaiMingwei on 16/11/8.
//  Copyright © 2016年 daviddai. All rights reserved.
//

#ifndef IGESS_hpp
#define IGESS_hpp
#include <stdio.h>
#include <math.h>
#include <string>
#include <stdlib.h>
#include "IGESS_aux.hpp"

#define type_int 0
#define type_double 1

struct PairCORAUC{
  double cor;
  double auc;
};

template<typename T>
void test_template(T* t);


IGESSfit* iGess(void* lpfX, vec y, int P, mat* lpsummaryinfo = NULL, Options* opt = NULL, int type = type_int, bool efficient = true, bool verbose = false);


double iGessCV(void* lpfX, vec y, int P, mat* lpsummaryinfo = NULL, Options* opt = NULL,int type = type_int, bool efficient = true, std::string measure = "auc");
arma::Col<uword> cross_valind(arma::uword N, arma::uword nfold);

template<typename T>
void igess_update(T* x_j, double* gamma, double* mu, double d, double s, double logPi, double sigma2beta, double sigma2e, int N, double xy, double* ytilde_pt, double* lpsummay = NULL, vec* lpparam = NULL, double xmean_j = 0, double xsd_j = 1, int type = type_int, bool befloat = true);//,

template<class T>
void addX (double* y, double a, T* x, int n);

template <typename T>
double dotX (T* x, double* y, int n);
#endif /* IGESS_hpp */
