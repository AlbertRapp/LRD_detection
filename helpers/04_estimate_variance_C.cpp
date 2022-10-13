#include <Rcpp.h>
#include <iostream>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector estimate_variance_C(NumericVector Y, int kMax) {
  int TMax = Y.size();
  NumericVector sample_var(kMax);
  
  for(int k = 0; k < kMax; ++k) {
    int n_blocks = TMax - k;
    
    NumericVector block_means(n_blocks);
    for (int i = 0; i < n_blocks; ++i) {
      NumericVector subvector(k + 1);
      for (int j = 0; j <= k; ++j) {
        subvector[j] = Y[i + j];
      }
      block_means[i] = mean(subvector);
    }
    
    double all_mean = mean(block_means);
    double tmp_sum = 0;
    for (int j = 0; j < n_blocks; ++j) {
      tmp_sum += pow(block_means[j] - all_mean, 2);
    }
    sample_var[k] = tmp_sum / n_blocks;
  }
  return sample_var;
}