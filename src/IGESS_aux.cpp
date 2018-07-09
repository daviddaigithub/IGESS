//
//  IGESS_aux.cpp
//  IGESSArma
//
//  Created by DaiMingwei on 16/11/8.
//  Copyright © 2016年 daviddai. All rights reserved.
//

#include "IGESS_aux.hpp"

double calauc(arma::vec label, arma::vec pred){
    double auc = 0;
    double m = mean(label);
    vec label2 = label;
    label2(find(label >= m)).fill(1);
    label2(find(label <= m)).fill(0);
    label = label2;
    uword N = pred.size();
    uword N_pos = sum(label);
    uvec  idx = sort_index(pred,"descend");
    vec above = linspace<vec>(1, N, N) - cumsum(label(idx));
    auc = (1 - sum(above % label(idx)) / (N_pos * (N-N_pos)));
    auc = auc > 0.5?auc:(1-auc);
    return auc;
}


/*update the parameters of the gamma distributions*/
void update_betaparam(uword P, uword K, double gamma_sum, double * lpsummary, double* lpgamma, vec* lpparams)
{
    if(K == 0) return;
    vec alphalogpvec(K);
    alphalogpvec.fill(0);
    double* lpalphalogp = alphalogpvec.memptr();
    for(uword k=0; k < K; k++)
    {
        for (uword j = 0; j < P; j++)
        {
            double pvalue = *(lpsummary + K * j + k);
            lpalphalogp[k] += (*(lpgamma + j))*(-log(pvalue));
        }
        (*lpparams)[k] = gamma_sum / lpalphalogp[k];
        if( (*lpparams)[k] > 1){
           (*lpparams)[k] = 1;
        }
    }

}


vec IGESSfit::predict(mat* X)
{
  //  uword N = X -> n_rows;
  //  mat Z0(N, 1, fill::ones);
    vec yhat = this->cov + (*X) * conv_to<vec>::from(this->gammas % this->mu);
    return yhat;
}

void IGESSfit::cal_powerfdr(DataModel* model, double threshold, PerformanceObject* po)
{
    vec gFDR = fdr2FDR(this->gammas);
    uvec ufound = find(gFDR < threshold);
    uword nerr = sum((*model -> labels)(ufound) == 0);
    uword nfound = ufound.size();

    double FDR = (nfound != 0) ? nerr * 1.0 / nfound : 0;

    uword nCausal = sum((*model ->labels) != 0);
    double power =  (ufound.size() - nerr) * 1.0  / nCausal ;

    po -> FDR = FDR;
    po -> power = power;

}

double IGESSfit::cal_auc(DataModel* model){
    vec label = conv_to<vec>::from((*model -> labels));
    double auc = calauc(label, this -> gammas);
    return auc;
}




double lb_pvalue(uword P, uword K, double * lpsummary, double* lpgamma, vec* lpparams){
    double lb = 0;
    if(K == 0) return lb;
    double* lpparam = lpparams -> memptr();
    double pvalue = 0;
    for (uword k = 0; k < K; k++) {
        double param_k = *(lpparam + k);
        for (uword j = 0; j < P; j++){
            pvalue =  *(lpsummary + j * K + k);
            lb += (*(lpgamma + j)) * (log(param_k) + (param_k - 1)*log(pvalue));
        }
    }
    return lb;

}

void update_param(uword N, uword P, Vardist vardist, double sumyytilde, vec diagXTX, double& sigma2e, double& sigma2beta, double& pi_p){
    vec term1 = vardist.gamma % (vardist.sigma2beta + square(vardist.mu));
    double term2 = sum((term1 - square(vardist.gamma % vardist.mu)) % diagXTX);

    sigma2e = (sumyytilde + term2) / N;
    double sum_vardist_gamma = sum(vardist.gamma);
    pi_p = sum_vardist_gamma / P;
    sigma2beta = sum(term1) / sum_vardist_gamma;
}

double dotXX (double* x, float* y, uword n) {
    double z = 0;
    for (uword i = 0; i < n; i++)
        z += x[i] * y[i];
    return z;
}

double dotXX (float* x, float* y, uword n) {
    double z = 0;
    for (uword i = 0; i < n; i++)
        z += x[i] * y[i];
    return z;
}


//lower bound for the linear part
double lb_linear(vec ytilde, vec diagXTX, vec y, double sigma2e, Vardist vardist){
    uword n = y.n_elem;
    double lb1 = - (0.5*n)*log(2 * M_PI * sigma2e);
    double lb2 = - 0.5 * sum(square(y-ytilde))/sigma2e;
    double lb3 =
    -0.5*sum( (vardist.gamma % (vardist.sigma2beta + square(vardist.mu)) -
               square(vardist.gamma % vardist.mu)) % diagXTX )/sigma2e;
    return lb1 + lb2 + lb3;

};

vec logpexp(vec x) {
    vec y = x;
    uvec idx = (x < 16);
    y(idx) = log(1 + exp(x(idx)));
    return y;
}

double logpexp(double x) {
    double y = log(1 + exp(x));
    return y;
}

//lower bound for gamma part
double lb_gamma(vec gamma, double log_pi){
    return sum((gamma - 1) * log_pi + (-logpexp(-log_pi)));
};

//lower bound for the klbeta part
double lb_klbeta(Vardist vardist, double sigma2beta){
    double lb = -sum(vardist.gamma % log(vardist.gamma+(vardist.gamma==0)) + (1-vardist.gamma) % log(1-vardist.gamma+(vardist.gamma == 1))) + 0.5*sum(vardist.gamma % (1+log(vardist.sigma2beta / sigma2beta)- (square(vardist.mu) + vardist.sigma2beta)/ sigma2beta ));
    return lb;
}






//convert local fdr to Global FDR
vec fdr2FDR(vec fdr){
    uword M = fdr.size();
    uvec indices = sort_index(fdr);
    vec sort_fdr = fdr(indices);
    vec FDR = cumsum(sort_fdr) / linspace<vec>(1, M, M);
    FDR.elem(find(FDR  > 1)).ones();
    FDR(indices) = FDR;
    return FDR;
}

void get_matrix_pointer(void* lpfX, int N, int P, int type, bool efficient,Mat<double> & Mat_f, void* & X_mat, bool & befloat){
  if(type == type_double){
    clock_t t1 = clock();
    Mat<double>* X_p = new Mat<double>( static_cast<double *>(lpfX) , N, P, false);
   //  if(type == type_double){
   // //   Mat_f  = *X_p;
   // //   delete X_p;
   //  }
   //  // else{
   //  //   Mat<int>* X_p = new Mat<int>( static_cast<int *>(lpfX) , N, P, false);
   //  //   Mat_f  = conv_to<fmat>::from(*X_p);
   //  //   X_p -> clear();
   //  //   delete X_p;
   //  // }
    X_mat = X_p;
    befloat = true;
    // cout <<" Time Elapsed " << 1.0*(clock() - t1)/CLOCKS_PER_SEC << " seconds"<<endl;
  }else if(type == type_int){
    Mat<int>* X_p = new Mat<int>( static_cast<int *>(lpfX) , N, P, false);
    X_mat  = X_p;
    befloat = false;
  }
}

Mat<double> cal_means(void* X_mat, bool befloat, int N){
  mat sumX;
  if( befloat ){
    sumX = conv_to<mat>::from(sum(*static_cast<Mat<double> *>(X_mat), 0));
  }
  else{
    sumX = conv_to<mat>::from(sum(*static_cast<Mat<int> *>(X_mat), 0));
  }

  mat SZX =  sumX * 1.0 / N ;//conv_to<mat>::from(mean(X, 0)); //column means of X
  return SZX;
}


void centering(void* X_mat, vec& y, double& mean_y, bool befloat, int N, int P, double* x_mean, double* x_sd){
  cal_means(X_mat, befloat, N, x_mean, x_sd);
  mean_y = as_scalar(mean(y)); //mean of y
  y -= mean_y;
  if ( befloat ){
    Mat<double>* m = static_cast<Mat<double> *>(X_mat);
    for(int i = 0; i < P; i++){
      m -> col(i) -= x_mean[i];
      m -> col(i) /= x_sd[i];
    }
  }
}

template <typename T>
void calMeanSd(T * data, int N, int M, double* mean, double* sd){
  int value;
  for(int j = 0; j < M; j++){
    mean[j] = 0;
    sd[j] = 0;
    for (int i = 0; i < N; i++) {
      value = data[j * N + i];
      mean[j] += value;
      sd[j] += value * value;
    }
    sd[j] = sqrt(sd[j] / N - (mean[j] / N) * (mean[j] / N));
    mean[j] = mean[j] / N;
  }
}

void cal_means(void* X_mat, bool befloat, int N, double* mean, double* sd){
  mat sumX;
  if( befloat ){
    // sumX = conv_to<mat>::from(sum(*static_cast<Mat<double> *>(X_mat), 0));
    Mat<double>* m = static_cast<Mat<double> *>(X_mat);
    calMeanSd(m->memptr(), N, m->n_cols, mean, sd);
  }
  else{
    Mat<char>* m = static_cast<Mat<char> *>(X_mat);
    calMeanSd(m->memptr(), N, m->n_cols, mean, sd);

    // sumX = conv_to<mat>::from(sum(*static_cast<Mat<int> *>(X_mat), 0));
  }
  // mat SZX =  sumX * 1.0 / N ;//conv_to<mat>::from(mean(X, 0)); //column means of X
  // return SZX;
}


vec cal_diagXTX(const vec& y, void* X_mat, bool befloat, const mat& SZX, const mat& SDX, int N, mat& xty){
  double sum_y = sum(y);
  //mat xty;
  vec diagXTX;
  // t1 = clock();
  if( befloat ){
    Mat<double> * mat_f = static_cast<Mat<double> *>(X_mat);
    xty = y.t( ) * (* mat_f);
    diagXTX = conv_to<vec>::from( sum(*mat_f % (*mat_f), 0) );
  }else{
    Mat<char> * mat_i = static_cast<Mat<char> *>(X_mat);
    xty = y.t() * (* mat_i) - sum_y * SZX;
    xty /= SDX;
    mat diag(xty.n_elem,1);
    diag.fill(N);
    diagXTX = diag;
  }
  return diagXTX;

}




