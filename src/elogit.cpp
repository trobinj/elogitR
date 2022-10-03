#include <RcppArmadillo.h>
#include <roptim.h>

#include "miscfunc.h"

using namespace roptim;

// [[Rcpp::depends(RcppArmadillo,roptim)]]

class ghquadNormal 
{

private:

  int n;
  arma::vec node;
  arma::vec wght;

public:

  ghquadNormal(int n) : n(n) 
  {
    const double sqrt2 = arma::datum::sqrt2;

    arma::mat J(n, n, arma::fill::zeros);
    for (int i = 0; i < n - 1; ++i) {
      J(i + 1, i) = sqrt((i + 1) / 2.0);
      J(i, i + 1) = J(i + 1, i);
    }
    arma::vec eigenval(n);
    arma::mat eigenvec(n, n);
    arma::eig_sym(eigenval, eigenvec, J);

    node.set_size(n);
    wght.set_size(n);
    
    for (int i = 0; i < n; ++i) {
      node(i) = eigenval(i) * sqrt2;
      wght(i) = pow(eigenvec(0, i), 2) / 
        pow(arma::norm(eigenvec.col(i)), 2);
    }
  } 

  arma::vec wghts() { return wght; }
  arma::vec nodes() { return node; }
  arma::vec nodes(double m, double s) { return node * s + m; }
};

void elparms(const arma::vec& theta, arma::vec& alph, arma::vec& gamm, arma::vec& beta) 
{
  int K = alph.n_elem;
  int p = beta.n_elem;

  alph(0) = 0.0;
  alph.subvec(1, K - 1) = theta.subvec(0, K - 2);
  gamm(0) = 0.0;
  gamm.subvec(1, K - 2) = theta.subvec(K - 1, 2 * K - 4);
  gamm(K - 1) = 1.0;
  beta = theta.subvec(2 * K - 3, 2 * K + p - 4);
}

void rmparms(const arma::vec& theta, arma::vec& alph, arma::vec& gamm, arma::vec& beta) 
{
  int K = alph.n_elem;
  int p = beta.n_elem;

  alph(0) = 0.0;
  alph(1) = 0.0;
  alph.subvec(2, K - 1) = -cumsum(theta.subvec(0, K - 3));
  gamm = arma::regspace(0, K - 1);
  beta = theta.subvec(K - 2, K + p - 3);
}

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

class elobs {

 private:

  arma::rowvec y;
  arma::rowvec x;
  arma::mat znum;
  arma::mat zden;
  int Kmax;

 public:

  elobs(arma::rowvec y, arma::rowvec x, int Kmax) : y(y), x(x), Kmax(Kmax) 
  {
    int K = y.n_cols;
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

  double logl(arma::vec& alph, arma::vec& gamm, arma::vec& beta) 
	{
    arma::vec eta = alph + gamm * x * beta;
    return sum(znum * eta) - sum(log(zden * exp(eta)));
  }

  void setnode(double z) { x(0) = z; }

  arma::rowvec gety() { return y; }
  arma::rowvec getx() { return x; }
};

class flobs {

 private:

  std::vector<elobs> data;
  int n;

 public:

  flobs(arma::rowvec y, arma::rowvec x, arma::mat yset) {
    int K = y.n_cols;

    arma::mat z =
        yset.rows(find(yset.col(0) == y(0) && yset.col(K - 1) == y(K - 1)));

    n = z.n_rows;
    data.reserve(n);
    for (int i = 0; i < n; ++i) {
      data.emplace_back(z.row(i), x, K);
    }
  }

  double logl(arma::vec& alph, arma::vec& gamm, arma::vec& beta) {
    double prob = 0.0;
    for (int i = 0; i < n; ++i) {
      prob = prob + exp(data[i].logl(alph, gamm, beta));
    }
    return log(prob);
  }

  void setnode(double z) {
    for (auto& datum : data) {
      datum.setnode(z);
    }
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

  eldata(arma::mat y, arma::mat x) {
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

  double operator()(const arma::vec& theta) override 
  {
    arma::vec alph(K);
    arma::vec gamm(K);
    arma::vec beta(p);
    elparms(theta, alph, gamm, beta);

    double logl = 0.0;
    for (int i = 0; i < n; ++i) {
      logl = logl + data[i].logl(alph, gamm, beta);
    }
    return -logl;
  }
};

template <typename obstype>
class elblock 
{
private:

  std::vector<obstype> data;
  arma::vec node, wght;
  int n, pnts;

public:

  elblock(arma::mat y, int Kmax, ghquadNormal ghdata) 
  {
    n = y.n_rows;

    arma::mat x(n, n + 1);
    x = join_rows(arma::ones(n, 1), -arma::eye(n, n));

    data.reserve(n);
    for (int i = 0; i < n; ++i) {
      data.emplace_back(y.row(i), x.row(i), Kmax);
    }

    node = ghdata.nodes();
    wght = ghdata.wghts();
    pnts = wght.n_elem;
  }

  elblock(arma::mat y, arma::mat yset, ghquadNormal ghdata) 
  {
    n = y.n_rows;

    arma::mat x(n, n + 1);
    x = join_rows(arma::ones(n, 1), -arma::eye(n, n));

    data.reserve(n);
    for (int i = 0; i < n; ++i) {
      data.emplace_back(y.row(i), x.row(i), yset);
    }

    node = ghdata.nodes();
    wght = ghdata.wghts();
    pnts = wght.n_elem;
  }

  double logl(arma::vec& alph, arma::vec& gamm, arma::vec& beta) 
  {
    double logl;
    double prob = 0.0;
    for (int q = 0; q < pnts; ++q) {
      logl = 0.0;
      for (int i = 0; i < n; ++i) {
        data[i].setnode(node(q));
        logl = logl + data[i].logl(alph, gamm, beta);
      }
      prob = prob + wght(q) * exp(logl);
    }
    return log(prob);
  }
};

template <typename obstype>
class rmdata : public Functor 
{

private:

  std::vector<elblock<obstype>> data;
  int n, m, K;

public:

  rmdata(arma::mat y, int m, int p, int Kmax) : m(m) 
  {
    n = y.n_rows / m;
    K = y.n_cols;

    ghquadNormal ghdata(p);

    data.reserve(n);
    for (int i = 0; i < n; ++i) {
      data.emplace_back(y.rows(i * m, i * m + m - 1), Kmax, ghdata);
    }
  }

  rmdata(arma::mat y, int m, int p) : m(m) 
  {
    n = y.n_rows / m;
    K = y.n_cols;

    ghquadNormal ghdata(p);

    elresponses responses(K);
    arma::mat yset = responses.elset();

    data.reserve(n);
    for (int i = 0; i < n; ++i) {
      data.emplace_back(y.rows(i * m, i * m + m - 1), yset, ghdata);
    }
  }

  double operator()(const arma::vec& theta) override 
	{
    arma::vec alph(K);
    arma::vec gamm(K);
    arma::vec beta(m + 1);
    rmparms(theta, alph, gamm, beta);

    double logl = 0.0;
    for (int i = 0; i < n; ++i) {
      logl = logl + data[i].logl(alph, gamm, beta);
    }
    return -logl;
  }
};

// [[Rcpp::export]]
arma::umat stm_samp(arma::vec theta, arma::mat x, int K) 
{
  int n = x.n_rows;
  int p = x.n_cols;

  arma::vec alph = theta.subvec(0, K - 1);
  arma::vec gamm = theta.subvec(K, 2 * K - 1);
  arma::vec beta = theta.subvec(2 * K, 2 * K + p - 1);

  arma::umat y(n, K);
  for (int i = 0; i < n; ++i) {
    y.row(i) = obs_samp(alph + gamm * (x.row(i) * beta)).t();
  }

  return y;
}

// [[Rcpp::export]]
arma::mat elpatterns(int K, int Kmax) 
{
  elresponses data(K);
  return data.elset(Kmax);
}

// [[Rcpp::export]]
arma::mat flpatterns(int K) 
{
  elresponses data(K);
  return data.flset();
}

// [[Rcpp::export]]
double elprofile(double z, arma::vec theta, arma::rowvec y, arma::rowvec x, int Kmax) 
{
  int p = x.n_elem;
  int K = y.n_elem;

  arma::vec alph = theta.subvec(0, K - 1);
  arma::vec gamm = theta.subvec(K, 2 * K - 1);
  arma::vec beta = theta.subvec(2 * K, 2 * K + p - 1);

  elobs data(y, x, Kmax);
  data.setnode(z);

  return data.logl(alph, gamm, beta);
}

// [[Rcpp::export]]
double flprofile(double z, arma::vec theta, arma::rowvec y, arma::rowvec x) 
{
  int p = x.n_elem;
  int K = y.n_elem;

  arma::vec alph = theta.subvec(0, K - 1);
  arma::vec gamm = theta.subvec(K, 2 * K - 1);
  arma::vec beta = theta.subvec(2 * K, 2 * K + p - 1);

  elresponses responses(K);
  arma::mat yset = responses.elset();

  flobs data(y, x, yset);
  data.setnode(z);

  return data.logl(alph, gamm, beta);
}

// [[Rcpp::export]]
Rcpp::List elstm(arma::vec theta, arma::mat y, arma::mat x, int Kmax, bool hessian) 
{
  using namespace roptim;
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
Rcpp::List flstm(arma::vec theta, arma::mat y, arma::mat x, bool hessian) 
{
  using namespace roptim;
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

// [[Rcpp::export]]
Rcpp::List elrsm(arma::vec theta, arma::mat y, int m, int p, int Kmax, bool hessian) 
{
  using namespace Rcpp;
  using namespace roptim;

  rmdata<elobs> data(y, m, p, Kmax);

  Roptim<rmdata<elobs>> opt("BFGS");
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
    Named("stderror") = wrap(stderror)
  );
}

// [[Rcpp::export]]
Rcpp::List flrsm(arma::vec theta, arma::mat y, int m, int p, bool hessian) 
{
  using namespace Rcpp;

  rmdata<flobs> data(y, m, p);

  Roptim<rmdata<flobs>> opt("BFGS");
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
    Named("stderror") = wrap(stderror)
  );
}

