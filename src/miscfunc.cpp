#include <RcppArmadillo.h>

arma::mat uniquerows(arma::mat x, double tol = 1e-8)
{
	int n = x.n_rows;
	arma::uvec indx(n, arma::fill::ones);

	for (int i = 0; i < n - 1; ++i) {
		for (int j = i + 1; j < n; ++j) {
			if (approx_equal(x.row(i), x.row(j), "absdiff", tol)) {
				indx(j) = 0;
			}
		}
	}

	return x.rows(find(indx == 1));
}

int rdiscrete(arma::vec wght)
{
	int n = wght.n_elem;
	arma::vec prob = wght / accu(wght);
	double u = R::runif(0.0, 1.0);
	double cprb = 0.0;
	for (int y = 0; y < n - 1; ++y) {
		cprb = cprb + prob(y);
		if (u < cprb)
		{
			return y;
		}
	}
	return n - 1;
}

void heapPermutations(arma::vec &x, int n, arma::mat &y)
{
	if (n == 1) {
		y = join_vert(y, x.t());
	}
	else {
		for (int i = 0; i < n; ++i) {
			heapPermutations(x, n - 1, y);
			if (n % 2 == 1) {
				std::swap(x(0), x(n - 1));
			}
			else {
				std::swap(x(i), x(n - 1));
			}
		}
	}
}
arma::mat heapPermutations(arma::vec x)
{
	arma::mat y;
	arma::vec z(x);
	int n = z.n_elem;
	heapPermutations(z, n, y);
	return y;
}
arma::mat heapPermutations(int a, int b)
{
	arma::mat y;
	arma::vec z = arma::regspace(a, b);
	int n = z.n_elem;
	heapPermutations(z, n, y);
	return y;
}

bool ordered(arma::vec x) {
  if (x.is_sorted("ascend")) {
    return true;
  }
  if (x.is_sorted("descend")) {
    return true;
  }
  int n = x.n_elem;
  arma::vec left, rght;
  for (int i = 1; i < n - 1; ++i) {
    if (x(i) == 1) {
      left = x.head(i);
      rght = x.tail(n - i - 1);
      if (left.is_sorted("descend") && rght.is_sorted("ascend")) {
        return true;
      } else {
        return false;
      }
    }
  }
  return false;
}

arma::uvec obs_samp(arma::vec eta) 
{
  int K = eta.n_elem;

  arma::uvec y(K);
  arma::vec delt(K);
  arma::vec eeta = exp(eta);

  y(0) = rdiscrete(eeta) + 1;
  int a = y(0) - 2;
  int b = y(0);

  for (int k = 1; k < K; ++k) {
    delt.zeros();
    if (a >= 0) {
      delt(a) = 1;
    }
    if (b <= K - 1) {
      delt(b) = 1;
    }
    y(k) = rdiscrete(delt % eeta) + 1;
    if (y(k) < y(k - 1)) {
      a = a - 1;
    } else {
      b = b + 1;
    }
  }

  return y;
}