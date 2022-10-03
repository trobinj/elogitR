#include <RcppArmadillo.h>
#include <roptim.h>

#include "miscfunc.h"

using namespace roptim;

// [[Rcpp::depends(RcppArmadillo,roptim)]]

namespace mlogitspc {

class elresponses 
{
private:

  int K;
  arma::mat z;

public:

  elresponses(int K) : K(K)
  {
    arma::mat y = heapPermutations(1, K);
    int n = y.n_rows;
    arma::vec u;
    for (int i = 0; i < n; ++i) {
      if (ordered((y.row(i)).t())) {
        u = arma::conv_to<arma::vec>::from(sort_index(y.row(i)));
        z = join_vert(z, u.t() + 1);
      }
    }
  }

  arma::mat elset() 
  {
    return z; 
  }

  arma::mat elset(int Kmax) 
  {
    int m = z.n_cols;
    arma::mat u(z);
    for (int i = Kmax; i < m; ++i) {
      u.col(i).fill(0.0);
    }
    return uniquerows(u);
  }

  arma::mat flset() 
  {
    int m = z.n_cols;
    arma::mat u(z);
    for (int i = 1; i < m - 1; ++i) {
      u.col(i).fill(0.0);
    }
    return uniquerows(u);
  }
};

class elobs 
{
private:

  arma::rowvec y;
  arma::rowvec x;
  arma::mat znum;
  arma::mat zden;
  int Kmax, K, p;

public:

  elobs(arma::rowvec y, arma::rowvec x, int Kmax) : y(y), x(x), Kmax(Kmax)
  {
    K = y.n_cols;
    p = x.n_cols;

    zden.zeros(K, K);
    znum.zeros(K, K);
    zden.row(0) = arma::ones(K).t();
    int a = y(0) - 2;
    int b = y(0);
    znum(0, y(0) - 1) = 1;

    for (int k = 1; k < K; ++k) {
      if (k > Kmax - 1) {
        znum = znum.rows(0, k - 1);
        zden = zden.rows(0, k - 1);
        break;
      }
      if (a >= 0) {
        zden(k, a) = 1;
      } else {
        znum = znum.rows(0, k - 1);
        zden = zden.rows(0, k - 1);
        break;
      }
      if (b <= K - 1) {
        zden(k, b) = 1;
      } else {
        znum = znum.rows(0, k - 1);
        zden = zden.rows(0, k - 1);
        break;
      }
      if (y(k) < y(k - 1)) {
        a = a - 1;
      } else {
        b = b + 1;
      }
      znum(k, y(k) - 1) = 1;
    }
  }

  double logl(const arma::vec &theta) 
  {
    arma::mat beta(theta);
    beta.reshape(K - 1, p);

    arma::vec eta(K);
    for (int i = 1; i < K; ++i) {
      eta(i) = dot(x, beta.row(i - 1));
    }

    return sum(znum * eta) - sum(log(zden * exp(eta)));
  }
};

class flobs 
{
private:

  std::vector<elobs> data;
  int n;

public:

  flobs(arma::rowvec y, arma::rowvec x, arma::mat yset) 
  {
    int K = y.n_cols;

    arma::mat z = yset.rows(find(yset.col(0) == y(0) && 
      yset.col(K - 1) == y(K - 1)));

    n = z.n_rows;
    data.reserve(n);
    for (int i = 0; i < n; ++i) {
      data.emplace_back(z.row(i), x, K);
    }
  }

  double logl(const arma::vec &theta) 
  {
    double prob = 0.0;
    for (int i = 0; i < n; ++i) {
      prob = prob + exp(data[i].logl(theta));
    }
    return log(prob);
  }
};

template <typename obstype> 
class eldata : public Functor 
{
private:

  std::vector<obstype> data;
  int n, p, K;

public:

  eldata(arma::mat y, arma::mat x, int Kmax) 
  {
    n = y.n_rows;
    p = x.n_cols;
    K = y.n_cols;

    data.reserve(n);
    for (int i = 0; i < n; ++i) {
      data.emplace_back(y.row(i), x.row(i), Kmax);
    }
  }

  eldata(arma::mat y, arma::mat x) 
  {
    n = y.n_rows;
    p = x.n_cols;
    K = y.n_cols;

    elresponses responses(K);
    arma::mat yset = responses.elset();

    data.reserve(n);
    for (int i = 0; i < n; ++i) {
      data.emplace_back(y.row(i), x.row(i), yset);
    }
  }

  double operator()(const arma::vec &theta) override 
  {
    double logl = 0.0;
    for (int i = 0; i < n; ++i) {
      logl = logl + data[i].logl(theta);
    }
    return -logl;
  }
};

}

// [[Rcpp::export]]
arma::umat mnl_samp(arma::vec theta, arma::mat x, int K)
{
  int n = x.n_rows;
  int p = x.n_cols;

  arma::umat y(n, K);
  arma::vec eta(K);
  arma::mat beta(theta);
  beta.reshape(K - 1, p); 

  for (int i = 0; i < n; ++i) {
    for (int j = 1; j < K; ++j) {
      eta(j) = dot(x.row(i), beta.row(j - 1));
    }
    y.row(i) = obs_samp(eta).t(); 
  }

  return y;
}

// [[Rcpp::export]]
Rcpp::List elmnl(arma::vec theta, arma::mat y, arma::mat x, int Kmax, bool hessian) 
{
  using namespace roptim;
  using namespace mlogitspc;
  using namespace Rcpp;

  eldata<elobs> data(y, x, Kmax);

  Roptim<eldata<elobs>> opt("BFGS");
  opt.set_hessian(hessian);
  opt.minimize(data, theta);

  if (opt.convergence() > 0) {
    Rcpp::Rcout << "convergence warning/error: " << opt.convergence() << "\n";
    Rcpp::Rcout << opt.message() << "\n";
  }

  arma::vec stderror(theta.n_elem);
  if (hessian) {
    stderror = sqrt(diagvec(inv(opt.hessian())));
  }

  return List::create(
    Named("estimate") = wrap(opt.par()),
    Named("stderror") = wrap(stderror),
    Named("aic") = 2 * opt.value() + 2 * theta.n_elem
  );
}

// [[Rcpp::export]]
Rcpp::List flmnl(arma::vec theta, arma::mat y, arma::mat x, bool hessian) 
{
  using namespace roptim;
  using namespace mlogitspc;
  using namespace Rcpp;

  eldata<flobs> data(y, x);

  Roptim<eldata<flobs>> opt("BFGS");
  opt.set_hessian(hessian);
  opt.minimize(data, theta);

  if (opt.convergence() > 0) {
    Rcpp::Rcout << "convergence warning/error: " << opt.convergence() << "\n";
    Rcpp::Rcout << opt.message() << "\n";
  }

  arma::vec stderror(theta.n_elem);
  if (hessian) {
    stderror = sqrt(diagvec(inv(opt.hessian())));
  }

  return List::create(
    Named("estimate") = wrap(opt.par()),
    Named("stderror") = wrap(stderror),
    Named("aic") = 2 * opt.value() + 2 * theta.n_elem
  );
}
