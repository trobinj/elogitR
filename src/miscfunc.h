#ifndef miscfunc_h
#define miscfunc_h

arma::mat uniquerows(arma::mat x, double tol = 1e-08);
int rdiscrete(arma::vec wght);
void heapPermutations(arma::vec &x, int n, arma::mat &y);
arma::mat heapPermutations(arma::vec x);
arma::mat heapPermutations(int a, int b);
bool ordered(arma::vec x);
arma::uvec obs_samp(arma::vec eta);

#endif