//
//  IGESS.cpp
//  IGESSArma
//
//  Created by DaiMingwei on 16/11/8.
//  Copyright © 2016年 daviddai,eeyang. All rights reserved.
//

#include "IGESS.hpp"



template <typename T>
double dotX (T* x, double* y, int n) {
  double sum = 0;
  //  #pragma omp parallel for reduction(+:sum)
  for (int i = 0; i < n; i++)
    sum += x[i] * y[i];
  return sum;
}



template<class T>
void addX (double* y, double a, T* x, int n) {
  //   #pragma omp parallel for num_threads(4)
  for (int i = 0; i < n; i++)
    y[i] += a * x[i];
}


template<typename T>
void igess_update(T* x_j, double* gamma, double* mu, double d, double s,  double logPi, double sigma2beta, double sigma2e, int N, double xy, double* ytilde_pt, double* lpsummay, vec* lpparam, double xmean_j, bool befloat){//

  double r = (*gamma) * (*mu);
  double numerator = xy + d * r -  dotX(x_j, ytilde_pt, N);
  mat y_tilde(ytilde_pt, N, 1, false);

  if(befloat == false){
     double sum_y_tilde = as_scalar(sum(y_tilde));
     numerator += xmean_j * sum_y_tilde;
  }
  (*mu) = s / sigma2e * numerator;

  double SSR_logratio(0);

  SSR_logratio = (*mu) * (*mu) / (2 * s) + 0.5 * log( s / sigma2beta);

  double ai = 0; //additional information provided by the summary statistisc
  if (lpsummay != NULL && lpparam != NULL) {
    uword K = lpparam -> size();
    double* lppara = lpparam -> memptr();
    for (uword i = 0; i < K; i++) {
      ai += log(lppara[i]) + (lppara[i] - 1) * log(lpsummay[i]);
    }
    SSR_logratio += ai;
  }


  (*gamma) = 1/(1+exp(-(logPi + SSR_logratio)));

  double rnew = (*gamma) * (*mu);

  addX (ytilde_pt, rnew - r, x_j, N);
  if(befloat == false){
      y_tilde -= (rnew - r) * xmean_j;
  }
}



IGESSfit* iGess(void* lpfX, vec y, int P, mat* lpsummaryinfo, Options* opt, int type, bool efficient, bool verbose)
{

  uword N = y.n_rows;
  void* X_mat;
  Mat<double> Mat_f;
  bool befloat;

  get_matrix_pointer(lpfX, N, P, type, efficient, Mat_f, X_mat, befloat);
  mat SZX; //store column means of X, return by function 'centering'
  double mean_y;// = as_scalar(mean(y)); //mean of y, return by function 'centering'

  centering(X_mat, y, mean_y, SZX, befloat, N, P);
  uword K = lpsummaryinfo != NULL ? (lpsummaryinfo -> n_rows) : 0;
  opt = opt != NULL ? opt : new Options();
  int max_iter = opt -> max_iter;
  int display_gap = opt -> display_gap;
  mat xty;
  vec diagXTX = cal_diagXTX(y, X_mat, befloat, SZX, N, xty);

  double pi_p = 0.01; //pi for prior proportion
  double mu0 = 0;
  double alpha0 = 0.5; //initial parameters of beta distribtuion

  Vardist vardist(P, mu0, pi_p);

  double sigma2e = var(y) / 2;
  double sigma2beta = sigma2e;

  vec beta = vardist.gamma % vardist.mu;
  vec ytilde = vec(N);//(*lpfX) * beta;
  ytilde.zeros();

  Col<double>* lpparams = NULL;
  if ( K > 0 ){
    lpparams = new Col<double>(K);
    lpparams -> fill(alpha0); //parameters for beta distribtuions
  }

  double L0 = -INFINITY;
  double L = 0;
  double* lpgamma = vardist.gamma.memptr();
  double* lpmu = vardist.mu.memptr();
  double* lpd = diagXTX.memptr();
  double* lpytilde = ytilde.memptr();
  double* lpxy = xty.memptr();
  double* lpsummary = lpsummaryinfo != NULL ? lpsummaryinfo -> memptr() : NULL;
  uword iter = 0;
  double* mean_x = SZX.memptr();


  for (iter = 0; iter < max_iter; iter ++) {
    clock_t t1 = clock();
    if(iter == 0)  cout <<"Begin Iterations" << endl;
    double logPi = log(pi_p / (1 - pi_p));
    double sigma2e_Sigma2beta = sigma2e / sigma2beta;
    vec xxsigma = diagXTX + sigma2e_Sigma2beta;
    vardist.sigma2beta = sigma2e / xxsigma;
    double* S = vardist.sigma2beta.memptr();
    double gamma_sum = 0;
    if(befloat){
      Mat<double> * mat_f = static_cast<Mat<double> *>(X_mat);
      double* lp_Xf = mat_f -> memptr();
      for (int j = 0; j < P; j++) {
        igess_update(lp_Xf + N*j, lpgamma + j, lpmu + j, lpd[j], S[j], logPi, sigma2beta, sigma2e, (int)N,  lpxy[j], lpytilde, lpsummary + K * j, lpparams, mean_x[j], befloat);
        gamma_sum += *(lpgamma + j);
      }
    }else{
      Mat<int> * mat_i = static_cast<Mat<int> *>(X_mat);
      int* lp_Xi = mat_i -> memptr();
      for (int j = 0; j < P; j++) {
        igess_update(lp_Xi + N*j, lpgamma + j, lpmu + j, lpd[j], S[j], logPi, sigma2beta, sigma2e, (int)N,  lpxy[j], lpytilde, lpsummary + K * j, lpparams, mean_x[j], befloat);
        gamma_sum += *(lpgamma + j);
      }
    }

    update_betaparam(P, K, gamma_sum, lpsummary, lpgamma, lpparams);

    update_param(N, P, vardist, sum(square(y-ytilde)), diagXTX, sigma2e, sigma2beta, pi_p);

    L = lb_linear(ytilde, diagXTX,y, sigma2e, vardist) + lb_gamma(vardist.gamma, logPi) +lb_klbeta(vardist, sigma2beta)
           + lb_pvalue(P, K, lpsummary, lpgamma, lpparams);

    if(iter % display_gap == 0 && verbose){
      cout <<iter <<"th iteration L=" << L << ";sigma2e = " << sigma2e << "; sigma2beta = " << sigma2beta << " ;pi = " <<
        pi_p <<" time = " << (clock() - t1)*1.0/CLOCKS_PER_SEC << endl;
    }

    if(L < L0){
      cout << "Lowerbound decreasing, Error at iteration "<<iter <<"th iteration, diff = " << L-L0 << endl;
      break;
    }else if(L - L0 < 1e-5)//if(fabs((L - L0)/L0) < 1e-8) //
    {
      if(verbose){
        cout<<"Converge at " << iter << "th iteration, L = "<< L << endl;
      }
      break;
    }
    L0 = L;

  }
  mat v = vardist.gamma % vardist.mu;
  double cov = mean_y - as_scalar(SZX * conv_to<mat>::from(vardist.gamma % vardist.mu));
  IGESSfit* fit = new IGESSfit(N, P,  K, iter, L,  sigma2e, sigma2beta, pi_p, vardist.gamma, vardist.mu
                                 , vardist.sigma2beta, lpparams,  ytilde, cov);
  return fit;
}

/******************************************************************
 Function:       IGESSCV
Description:    Calculate the correlation and auc of the cross validation of the input data
Calls:
Input:
mat* lpfX,   //(*lpfX) is a N by P matrix of float,
N denotes the number of samples,
P is the number of SNPs
vec y,       //y is N by 1 vector of float corresponding to phenotype of each individual
mat* Z,      // (*Z) is a N by M matrix of float,  M is the number of covariates for each individual
//default value : NULL
mat* lpsummaryinfo,   //(*lpsummaryinfo) is a P by K matrix, K is the number of GWAS with Summary Statistics
//default value : NULL
Options* opt          // the options for the functions
//default value : NULL
Return:
 PairCORAUC      // A Struct of the value of AUC and Correlation Value
******************************************************************/
double iGessCV(void* lpfX, vec y, int P, mat* lpsummaryinfo, Options* opt, int type, bool efficient, std::string measure){

 //   cout << "Start of the CV function....";
    uword N = y.n_rows;
    void* X_mat;
    Mat<double> Mat_f;
    bool befloat;
    get_matrix_pointer(lpfX, N, P, type, efficient, Mat_f, X_mat, befloat);
  //  cout <<"Outer:" <<Mat_f.at(0) << "  " << Mat_f.at(1) << endl;
    opt = opt != NULL ? opt : new Options();
    uword nfold = opt -> n_fold;
    Col<uword> indices = cross_valind(N, nfold);
    vec ylabel = y;
    vec predY(N);
//    Mat<float> Xf = Mat_f;

  //  vec diagXTX = conv_to<vec>::from( sum(Xf % Xf, 0) );
//    cout <<"diagXTX=" << diagXTX[0] <<"  " << diagXTX[1] << " " << endl;

    cout << "Start " << nfold <<" cross validations!" << endl;
    for (uword i = 1; i <= nfold; i++) {
        cout <<"."<< endl;
        Col<uword> train_idx = find(indices != i);
        Col<uword> test_idx = find(indices == i);
        if( befloat ){
          Mat<double> * mat_f = static_cast<Mat<double> *>(X_mat);
          Mat<double> trainM = mat_f -> rows(train_idx);
          vec ytrain = y(train_idx);
          Mat<double> testM = mat_f -> rows(test_idx);
          vec ytest = y(test_idx);
          double* trainX = trainM.memptr();
          IGESSfit* f = iGess(trainX, ytrain, P,  lpsummaryinfo, opt, type, efficient);
          vec predy = f -> predict(&testM);
          predY.elem(test_idx) = predy;
       //   double accuracy = calauc(conv_to<vec>::from(ylabel(test_idx)), conv_to<vec>::from(predy));
          delete f;
        }else{
          Mat<int> * mat_i = static_cast<Mat<int> *>(X_mat);
          Mat<int> trainM = mat_i -> rows(train_idx);
          vec ytrain = y(train_idx);
          Mat<double> testM = conv_to<mat>::from(mat_i -> rows(test_idx));
          vec ytest = y(test_idx);
          int* trainX = trainM.memptr();
          IGESSfit* f = iGess(trainX, ytrain, P,  lpsummaryinfo, opt, type, efficient);
          vec predy = f -> predict(&testM);
          predY.elem(test_idx) = predy;
          delete f;
        }
    }

    double predict = 0;
    if(measure.compare("auc") == 0 ){
       predict = calauc(conv_to<vec>::from(ylabel), conv_to<vec>::from(predY));
    }else if(measure.compare("mse") == 0 ){
       vec diff = y - predY;
       predict = mean(diff % diff);//as_scalar(cor(y, predY));
    }else if(measure.compare("cor") == 0 ){
       predict = as_scalar(cor(y, predY));
    }
    return predict;
}

/**shuffle the index for cross validation*/
arma::Col<uword> cross_valind(arma::uword N, arma::uword nfold){
    arma::Col<uword> indices(N);
    arma_rng::set_seed_random();
    arma::Col<uword> vec_n = arma::shuffle(arma::linspace < arma::Col <uword> >(1, N, N));

    indices.fill(nfold);
    for(uword n = 1; n <= nfold-1; n++){
        arma::Col<uword> in = vec_n.rows((n-1)*N/nfold,n*N/nfold - 1);
        indices.elem(in - 1 ).fill(n);
    }
    return indices;
}
