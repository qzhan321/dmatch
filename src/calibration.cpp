#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace std;


const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec dmvnrmArma(arma::mat x,
                     arma::rowvec mean,
                     arma::mat sigma,
                     bool logd = true) {
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;

  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;
    out(i) = constants - 0.5 * arma::sum(z%z) + rootisum;
  }

  if (logd == false) {
    out = exp(out);
  }
  return(out);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double likelihood_para_cpp(arma::mat X, arma::mat Y, arma::rowvec mu_x, arma::rowvec mu_y, arma::mat A,
                       arma::rowvec d, arma::mat D1, arma::mat D2, bool logd = true){

  double likelihood = 0;

  uint32_t p = Y.n_cols;
  arma::mat moveX = X*A;
  arma::mat moveY = Y*A;

  for(int j = 0; j < p-1; ++j){
    moveX.col(j) -= d(j);
    moveY.col(j) -= d(j);
  }

  likelihood += sum(dmvnrmArma(moveX, mu_x, D2, logd));
  likelihood += sum(dmvnrmArma(moveY, mu_y, D1, logd));

  return likelihood;

}

/*
 *  Internal C++ function for Mahalanobis distance
 */
arma::vec mahaInt(arma::mat & X,
                  arma::vec & mu,
                  arma::mat & sigma,
                  const bool isChol = false)
{
  using namespace arma;

  // Some sanity checks

  if(mu.n_elem != sigma.n_cols) Rcpp::stop("The mean vector has a different dimensions from the covariance matrix.");
  if(X.n_cols != sigma.n_cols)  Rcpp::stop("The number of columns of X is different from the dimension of the covariance matrix.");

  // Calculate transposed cholesky dec. unless sigma is alread a cholesky dec.
  mat cholDec;
  if( isChol == false ) {
    cholDec = trimatl(chol(sigma).t());
  }
  else{
    cholDec = trimatl(sigma.t());
    if(any(cholDec.diag() <= 0.0))  Rcpp::stop("The supplied cholesky decomposition has values <= 0.0 on the main diagonal.");
  }

  vec D = cholDec.diag();

  vec out(X.n_rows);

  // Declaring some private variables
  uint32_t d = X.n_cols;
  uint32_t n = X.n_rows;

  vec tmp(d);

  double acc;
  uint32_t icol, irow, ii;

  // For each of the "n" random vectors, forwardsolve the corresponding linear system.
  // Forwardsolve because I'm using the lower triangle Cholesky.

  for(icol = 0; icol < n; icol++)
  {

    for(irow = 0; irow < d; irow++)
    {
      acc = 0.0;

      for(ii = 0; ii < irow; ii++) acc += tmp.at(ii) * cholDec.at(irow, ii);

      tmp.at(irow) = ( X.at(icol, irow) - mu.at(irow) - acc ) / D.at(irow);
    }

    out.at(icol) = sum(square(tmp));
  }

  return out;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec dmvtArma(arma::mat X, arma::vec mu, arma::mat sigma, double df, bool log){

  unsigned int d = X.n_cols;
  using namespace arma;

  mat cholDec = arma::chol(sigma);

  vec out = mahaInt(X, mu, cholDec, true);

  uint32_t ii;
  uint32_t n = X.n_rows;
  double logDet = sum(arma::log(cholDec.diag()));
  double c = lgamma((d + df)/2.0) - (lgamma(df/2.0) + logDet + d/2.0 * std::log(M_PI * df));

  for(ii = 0; ii < n; ii++){
    out.at(ii) = c - 0.5 * (df + d) * log1p(out.at(ii)/df);
  }

  if (log == false) out = exp(out);

  return( out );
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat callmtCov(arma::mat X, int nu){
  Environment myEnv("package:SinglecellCalibration");
  Function mycov = myEnv["callmtCov_t"];
  Rcpp::NumericMatrix mycovRes = wrap(mycov(X, nu));
  return as<arma::mat>(mycovRes);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double likelihood_para_t_cpp(arma::mat X, arma::mat Y, arma::rowvec mu_x, arma::rowvec mu_y, arma::mat A,
                           arma::rowvec d, arma::mat Sigma1, arma::mat Sigma2, double df, bool logd = true){

  double likelihood = 0;

  uint32_t p = Y.n_cols;
  arma::mat moveX = X*A;
  arma::mat moveY = Y*A;

  for(int j = 0; j < p-1; ++j){
    moveX.col(j) -= d(j);
    moveY.col(j) -= d(j);
  }

  likelihood += sum(dmvtArma(moveX, mu_x, Sigma2, df, logd));
  likelihood += sum(dmvtArma(moveY, mu_y, Sigma1, df, logd));

  return likelihood;

}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat rotation_matrix(double theta){
  using namespace arma;
  double a1 = sin(theta); // * M_PI);
  double a2 = cos(theta); // * M_PI);
  arma::mat rot = arma::zeros<arma::mat>(2, 2);
  rot(0, 0) = rot(1, 1) = a2;
  rot(0, 1) = - a1;
  rot(1, 0) = a1;
  return rot;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat rotation_matrix_dev(double theta){
  using namespace arma;
  double a1 = sin(theta); // * M_PI
  double a2 = cos(theta); // * M_PI
  arma::mat rot = arma::zeros<arma::mat>(2, 2);
  rot(0, 0) = rot(1, 1) = - a1;
  rot(0, 1) = - a2;
  rot(1, 0) = a2;
  return rot;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List run_simple_est_cpp(arma::mat X, arma::mat Y, arma::mat x1, arma::mat y1, int steps = 1000){

  uint32_t n = X.n_rows;
  uint32_t m = Y.n_rows;
  uint32_t p = X.n_cols;

  arma::rowvec mu_x = mean(x1, 0);
  arma::rowvec mu_y = mean(y1, 0);
  arma::mat D1 = cov(y1);
  arma::mat D2 = cov(x1);

  arma::rowvec estd = (mean(X, 0) - mu_x + mean(Y, 0) - mu_y)/2;

  arma::mat X_t = X.t();
  arma::mat Y_t = Y.t();
  arma::mat XX = X_t*X;
  arma::mat YY = Y_t*Y;

  arma::mat D1_inv = inv(D1);
  arma::mat D2_inv = inv(D2);
  arma::mat D12_inv_sum = D1_inv + D2_inv;

  arma::mat omega1 = arma::zeros<arma::mat>(m, p);

  for(int j = 0; j < m-1; ++j){
    omega1.row(j) = mu_y + estd;
  }

  arma::mat omega2 = arma::zeros<arma::mat>(n, p);
  for(int j = 0; j < n-1; ++j){
    omega2.row(j) = mu_x + estd;
  }

  arma::mat beta1 = inv(YY)*Y_t*omega1;
  arma::mat beta2 = inv(XX)*X_t*omega2;

  arma::mat estA = (beta1 + beta2)/2;
  arma::mat C1 = YY*estA*D1_inv/m + XX*estA*D2_inv/n;
  arma::mat C2 = Y_t*omega1*D1_inv/m + X_t*omega2*D2_inv/n;
  arma::mat A_gradient = C1 - C2;
  arma::mat pre = arma::ceil(log10(abs(A_gradient)));
  arma::rowvec pre2 = estd;
  arma::rowvec d_gradient = estd;
  double start = 0;
  for(int kk = 0; kk < steps-1; ++kk){

    for(int j = 0; j < m-1; ++j){
      omega1.row(j) = mu_y + estd;
    }

    for(int j = 0; j < n-1; ++j){
      omega2.row(j) = mu_x + estd;
    }

    C1 = YY*estA*D1_inv/m + XX*estA*D2_inv/n;
    C2 = Y_t*omega1*D1_inv/m + X_t*omega2*D2_inv/n;
    A_gradient = C1 - C2;

    arma::mat XestA = X*estA;
    arma::mat YestA = Y*estA;
    arma::rowvec errorX = arma::zeros<arma::rowvec>(p);
    arma::rowvec errorY = arma::zeros<arma::rowvec>(p);
    for(int j = 0; j < n-1; ++j){
      errorX = errorX + XestA.row(j);
    }
    for(int j = 0; j < m-1; ++j){
      errorY = errorY + YestA.row(j);
    }
    arma::rowvec ave_errorX = errorX/n - mu_x;
    arma::rowvec ave_errorY = errorY/m - mu_y;
    arma::rowvec dmat_part = ave_errorX*D2_inv + ave_errorY*D1_inv;
    d_gradient = estd*D12_inv_sum - dmat_part;

    if(kk==1){
      pre = arma::ceil(log10(abs(A_gradient)));
      pre2 = arma::ceil(log10(abs(d_gradient)));
      if(pre2.max() > pre.max()){
        start = pre2.max();
      } else {
        start = pre.max();
      }
      start = pow(10, -start);
    //  Rcpp::Rcout << "start" << start << std::endl;
    }
   // pre = arma::ceil(log10(abs(A_gradient)));

   // pre2 = arma::ceil(log10(abs(d_gradient)));

    if(start >= 0.01){
      estA = estA - start*A_gradient/(kk+1);
      estd = estd - start*d_gradient/(kk+1);
    } else {
      estA = estA - start*A_gradient;
      estd = estd - start*d_gradient;
    }

  }

  return Rcpp::List::create(Rcpp::Named("A") = estA,
                            Rcpp::Named("d") = estd);

}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List run_simple_est_theta_cpp(arma::mat X, arma::mat Y, arma::mat x1, arma::mat y1, int steps = 500){

  uint32_t n = X.n_rows;
  uint32_t m = Y.n_rows;
  uint32_t p = X.n_cols;

  arma::rowvec mu_x = mean(x1, 0);
  arma::rowvec mu_y = mean(y1, 0);
  arma::mat D1 = cov(y1);
  arma::mat D2 = cov(x1);

  arma::rowvec estd = -(mean(X, 0) - mu_x + mean(Y, 0) - mu_y)/2;

  arma::mat X_t = X.t();
  arma::mat Y_t = Y.t();

  arma::mat D1_inv = inv(D1);
  arma::mat D2_inv = inv(D2);
  arma::mat D12_inv_sum = D1_inv + D2_inv;

  arma::rowvec pre2 = estd;
  arma::rowvec d_gradient = estd;
  arma::mat Theta_gradient = arma::zeros<arma::mat>(1, 1);
  double start = 0;
  double estTheta = 0;
  arma::mat estA = rotation_matrix(estTheta);

  for(int kk = 0; kk < steps-1; ++kk){

    arma::mat XestA = X*estA;
    arma::mat YestA = Y*estA;
    arma::rowvec errorX = arma::zeros<arma::rowvec>(p);
    arma::rowvec errorY = arma::zeros<arma::rowvec>(p);
    for(int j = 0; j < n-1; ++j){
      errorX = errorX + XestA.row(j);
    }
    for(int j = 0; j < m-1; ++j){
      errorY = errorY + YestA.row(j);
    }
    arma::rowvec ave_errorX = errorX/n - mu_x;
    arma::rowvec ave_errorY = errorY/m - mu_y;
    arma::rowvec dmat_part = ave_errorX*D2_inv + ave_errorY*D1_inv;
    d_gradient = estd*D12_inv_sum - dmat_part;

    arma::mat A_T = estA.t();
    arma::mat A_dev = rotation_matrix_dev(estTheta);
    arma::mat A_dev_T = A_dev.t();

    arma::mat omega1 = A_T*D1_inv*A_dev;
    arma::mat omega2 = A_T*D2_inv*A_dev;


    arma::rowvec mu_y_est_d = mu_y + estd;
    arma::rowvec mu_x_est_d = mu_x + estd;
    arma::mat omega3 = A_dev_T*mu_y_est_d.t();
    arma::mat omega4 = A_dev_T*mu_x_est_d.t();

    arma::mat Theta_gradient1 = arma::zeros<arma::mat>(1, 1);
    for(int j = 0; j < m-1; ++j){
      Theta_gradient1 +=  Y.row(j)*omega1*Y_t.col(j) - Y.row(j)*omega3;
    }

    arma::mat Theta_gradient2 = arma::zeros<arma::mat>(1, 1);
    for(int j = 0; j < n-1; ++j){
      Theta_gradient2 +=  X.row(j)*omega2*X_t.col(j) - X.row(j)*omega4;
    }
    Theta_gradient = Theta_gradient1/m + Theta_gradient2/n;

    arma::mat pre = arma::ceil(log10(abs(Theta_gradient)));
    pre2 = arma::ceil(log10(abs(d_gradient)));
    if(pre2.max() > pre.max()){
      start = pre2.max();
    } else {
      start = pre.max();
    }
    start = pow(10, -start);

    Rcpp::Rcout << "A_gradient " << Theta_gradient(0, 0) << std::endl;
    Rcpp::Rcout << "d_gradient " << d_gradient << std::endl;
//  if(start >= 0.01){
//      estTheta = estTheta - start*Theta_gradient(0, 0)/(kk+1);
//     estd = estd - start*d_gradient/(kk+1);
//    } else {
      estTheta = estTheta - start*Theta_gradient(0, 0)/(kk+1);
      estd = estd - start*d_gradient/(kk+1);
//    }

    estA = rotation_matrix(estTheta);

  }



  return Rcpp::List::create(Rcpp::Named("Theta") = estTheta,
                            Rcpp::Named("A") = estA,
                            Rcpp::Named("d") = estd);

}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List run_simple_est_theta_one_cluster_cpp(arma::mat X, arma::mat x1, int steps = 500){

  uint32_t n = X.n_rows;
  uint32_t p = X.n_cols;

  arma::rowvec mu_x = mean(x1, 0);
  arma::mat D2 = cov(x1);

  arma::rowvec estd = -(mean(X, 0) - mu_x)/2;

  arma::mat X_t = X.t();

  arma::mat D2_inv = inv(D2);

  arma::rowvec pre2 = estd;
  arma::rowvec d_gradient = estd;
  arma::mat Theta_gradient = arma::zeros<arma::mat>(1, 1);
  double start = 0;
  double estTheta = 0;
  arma::mat estA = rotation_matrix(estTheta);

  for(int kk = 0; kk < steps-1; ++kk){

    arma::mat XestA = X*estA;
    arma::rowvec errorX = arma::zeros<arma::rowvec>(p);
    for(int j = 0; j < n-1; ++j){
      errorX = errorX + XestA.row(j);
    }

    arma::rowvec ave_errorX = errorX/n - mu_x;
    arma::rowvec dmat_part = ave_errorX*D2_inv;
    d_gradient = estd*D2_inv - dmat_part;

    arma::mat A_T = estA.t();
    arma::mat A_dev = rotation_matrix_dev(estTheta);
    arma::mat A_dev_T = A_dev.t();

    arma::mat omega2 = A_T*D2_inv*A_dev;

    arma::rowvec mu_x_est_d = mu_x + estd;
    arma::mat omega4 = A_dev_T*mu_x_est_d.t();

    arma::mat Theta_gradient2 = arma::zeros<arma::mat>(1, 1);
    for(int j = 0; j < n-1; ++j){
      Theta_gradient2 +=  X.row(j)*omega2*X_t.col(j) - X.row(j)*omega4;
    }
    Theta_gradient = Theta_gradient2/n;

    if(kk==0){
      arma::mat pre = arma::ceil(log10(abs(Theta_gradient)));
      pre2 = arma::ceil(log10(abs(d_gradient)));
      if(pre2.max() > pre.max()){
        start = pre2.max();
      } else {
        start = pre.max();
      }
      start = pow(10, -start);
      //  Rcpp::Rcout << "start" << start << std::endl;
    }


    if(start >= 0.01){
      estTheta = estTheta - start*Theta_gradient(0, 0)/(kk+1);
      estd = estd - start*d_gradient/(kk+1);
    } else {
      estTheta = estTheta - start*Theta_gradient(0, 0);
      estd = estd - start*d_gradient;
    }

    estA = rotation_matrix(estTheta);

  }



  return Rcpp::List::create(Rcpp::Named("Theta") = estTheta,
                            Rcpp::Named("A") = estA,
                            Rcpp::Named("d") = estd);

}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List run_simple_est_theta_coord_cpp(arma::mat X, arma::mat Y, arma::mat x1, arma::mat y1, int steps = 500){

  uint32_t n = X.n_rows;
  uint32_t m = Y.n_rows;
  uint32_t p = X.n_cols;

  arma::rowvec mu_x = mean(x1, 0);
  arma::rowvec mu_y = mean(y1, 0);
  arma::mat D1 = cov(y1);
  arma::mat D2 = cov(x1);

  arma::rowvec estd = -(mean(X, 0) - mu_x + mean(Y, 0) - mu_y)/2;

  arma::mat X_t = X.t();
  arma::mat Y_t = Y.t();

  arma::mat D1_inv = inv(D1);
  arma::mat D2_inv = inv(D2);
  arma::mat D12_inv_sum = D1_inv + D2_inv;

  arma::rowvec pre2 = estd;
  arma::rowvec d_gradient = estd;
  arma::mat Theta_gradient = arma::zeros<arma::mat>(1, 1);
  double start1 = 0;
  double start2 = 0;
  double estTheta = 0;
  arma::mat estA = rotation_matrix(estTheta);

  for(int kk = 0; kk < steps-1; ++kk){

    arma::mat XestA = X*estA;
    arma::mat YestA = Y*estA;
    arma::rowvec errorX = arma::zeros<arma::rowvec>(p);
    arma::rowvec errorY = arma::zeros<arma::rowvec>(p);
    for(int j = 0; j < n-1; ++j){
      errorX = errorX + XestA.row(j);
    }
    for(int j = 0; j < m-1; ++j){
      errorY = errorY + YestA.row(j);
    }
    arma::rowvec ave_errorX = errorX/n - mu_x;
    arma::rowvec ave_errorY = errorY/m - mu_y;
    arma::rowvec dmat_part = ave_errorX*D2_inv + ave_errorY*D1_inv;
    d_gradient = estd*D12_inv_sum - dmat_part;

    if(kk==0){
      pre2 = arma::ceil(log10(abs(d_gradient)));
      start1 = pre2.max();
      start1 = pow(10, -start1);
      Rcpp::Rcout << "start1" << start1 << std::endl;
    }

    if(start1 >= 0.01){
      estd = estd - start1*d_gradient/(kk+1);
    } else {
      estd = estd - start1*d_gradient;
    }

    arma::mat A_T = estA.t();
    arma::mat A_dev = rotation_matrix_dev(estTheta);
    arma::mat A_dev_T = A_dev.t();

    arma::mat omega1 = A_T*D1_inv*A_dev;
    arma::mat omega2 = A_T*D2_inv*A_dev;

    arma::rowvec mu_y_est_d = mu_y + estd;
    arma::rowvec mu_x_est_d = mu_x + estd;
    arma::mat omega3 = A_dev_T*mu_y_est_d.t();
    arma::mat omega4 = A_dev_T*mu_x_est_d.t();

    arma::mat Theta_gradient1 = arma::zeros<arma::mat>(1, 1);
    for(int j = 0; j < m-1; ++j){
      Theta_gradient1 +=  Y.row(j)*omega1*Y_t.col(j) - Y.row(j)*omega3;
    }

    arma::mat Theta_gradient2 = arma::zeros<arma::mat>(1, 1);
    for(int j = 0; j < n-1; ++j){
      Theta_gradient2 +=  X.row(j)*omega2*X_t.col(j) - X.row(j)*omega4;
    }
    Theta_gradient = Theta_gradient1/m + Theta_gradient2/n;

    if(kk==0){
      arma::mat pre = arma::ceil(log10(abs(Theta_gradient)));
      start2 = pre.max();
      start2 = pow(10, -start2);
      Rcpp::Rcout << "start2" << start2 << std::endl;
    }

    if(start2 >= 0.01){
      estTheta = estTheta - start2*Theta_gradient(0, 0)/(kk+1);
    } else {
      estTheta = estTheta - start2*Theta_gradient(0, 0);
    }

    estA = rotation_matrix(estTheta);

  }



  return Rcpp::List::create(Rcpp::Named("Theta") = estTheta,
                            Rcpp::Named("A") = estA,
                            Rcpp::Named("d") = estd);

}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List run_simple_est_theta_three_cluster_cpp(arma::mat X, arma::mat Y, arma::mat Z, arma::mat x1, arma::mat y1, arma::mat z1, int steps = 500){

  uint32_t n = X.n_rows;
  uint32_t m = Y.n_rows;
  uint32_t l = Z.n_rows;
  uint32_t p = X.n_cols;

  arma::rowvec mu_x = mean(x1, 0);
  arma::rowvec mu_y = mean(y1, 0);
  arma::rowvec mu_z = mean(z1, 0);
  arma::mat D1 = cov(y1);
  arma::mat D2 = cov(x1);
  arma::mat D3 = cov(z1);

  arma::rowvec estd = -(mean(X, 0) - mu_x + mean(Y, 0) - mu_y + mean(Z, 0) - mu_z)/2;

  arma::mat X_t = X.t();
  arma::mat Y_t = Y.t();
  arma::mat Z_t = Z.t();

  arma::mat D1_inv = inv(D1);
  arma::mat D2_inv = inv(D2);
  arma::mat D3_inv = inv(D3);
  arma::mat D123_inv_sum = D1_inv + D2_inv + D3_inv;

  arma::rowvec pre2 = estd;
  arma::rowvec d_gradient = estd;
  arma::mat Theta_gradient = arma::zeros<arma::mat>(1, 1);
  double start = 0;
  double estTheta = 0;
  arma::mat estA = rotation_matrix(estTheta);

  for(int kk = 0; kk < steps-1; ++kk){

    arma::mat XestA = X*estA;
    arma::mat YestA = Y*estA;
    arma::mat ZestA = Z*estA;
    arma::rowvec errorX = arma::zeros<arma::rowvec>(p);
    arma::rowvec errorY = arma::zeros<arma::rowvec>(p);
    arma::rowvec errorZ = arma::zeros<arma::rowvec>(p);
    for(int j = 0; j < n-1; ++j){
      errorX = errorX + XestA.row(j);
    }
    for(int j = 0; j < m-1; ++j){
      errorY = errorY + YestA.row(j);
    }
    for(int j = 0; j < l-1; ++j){
      errorZ = errorZ + ZestA.row(j);
    }
    arma::rowvec ave_errorX = errorX/n - mu_x;
    arma::rowvec ave_errorY = errorY/m - mu_y;
    arma::rowvec ave_errorZ = errorZ/m - mu_z;
    arma::rowvec dmat_part = ave_errorX*D2_inv + ave_errorY*D1_inv + ave_errorZ*D3_inv;
    d_gradient = estd*D123_inv_sum - dmat_part;

    arma::mat A_T = estA.t();
    arma::mat A_dev = rotation_matrix_dev(estTheta);
    arma::mat A_dev_T = A_dev.t();

    arma::mat omega1 = A_T*D1_inv*A_dev;
    arma::mat omega2 = A_T*D2_inv*A_dev;
    arma::mat omega3 = A_T*D3_inv*A_dev;


    arma::rowvec mu_y_est_d = mu_y + estd;
    arma::rowvec mu_x_est_d = mu_x + estd;
    arma::rowvec mu_z_est_d = mu_z + estd;
    arma::mat omega4 = A_dev_T*mu_y_est_d.t();
    arma::mat omega5 = A_dev_T*mu_x_est_d.t();
    arma::mat omega6 = A_dev_T*mu_z_est_d.t();

    arma::mat Theta_gradient1 = arma::zeros<arma::mat>(1, 1);
    for(int j = 0; j < m-1; ++j){
      Theta_gradient1 +=  Y.row(j)*omega1*Y_t.col(j) - Y.row(j)*omega4;
    }

    arma::mat Theta_gradient2 = arma::zeros<arma::mat>(1, 1);
    for(int j = 0; j < n-1; ++j){
      Theta_gradient2 +=  X.row(j)*omega2*X_t.col(j) - X.row(j)*omega6;
    }
    arma::mat Theta_gradient3 = arma::zeros<arma::mat>(1, 1);
    for(int j = 0; j < l-1; ++j){
      Theta_gradient3 +=  Z.row(j)*omega3*Z_t.col(j) - Z.row(j)*omega6;
    }
    Theta_gradient = Theta_gradient1/m + Theta_gradient2/n + Theta_gradient3/l;

    if(kk==0){
      arma::mat pre = arma::ceil(log10(abs(Theta_gradient)));
      pre2 = arma::ceil(log10(abs(d_gradient)));
      if(pre2.max() > pre.max()){
        start = pre2.max();
      } else {
        start = pre.max();
      }
      start = pow(10, -start);
      //  Rcpp::Rcout << "start" << start << std::endl;
    }


    if(start >= 0.01){
      estTheta = estTheta - start*Theta_gradient(0, 0)/(kk+1);
      estd = estd - start*d_gradient/(kk+1);
    } else {
      estTheta = estTheta - start*Theta_gradient(0, 0);
      estd = estd - start*d_gradient;
    }

    estA = rotation_matrix(estTheta);

  }



  return Rcpp::List::create(Rcpp::Named("Theta") = estTheta,
                            Rcpp::Named("A") = estA,
                            Rcpp::Named("d") = estd);

}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List run_simple_est_coord_cpp(arma::mat X, arma::mat Y, arma::mat x1, arma::mat y1, int steps = 1000){

  uint32_t n = X.n_rows;
  uint32_t m = Y.n_rows;
  uint32_t p = X.n_cols;

  arma::rowvec mu_x = mean(x1, 0);
  arma::rowvec mu_y = mean(y1, 0);
  arma::mat D1 = cov(y1);
  arma::mat D2 = cov(x1);

  arma::rowvec estd = (mean(X, 0) - mu_x + mean(Y, 0) - mu_y)/2;

  arma::mat X_t = X.t();
  arma::mat Y_t = Y.t();
  arma::mat XX = X_t*X;
  arma::mat YY = Y_t*Y;

  arma::mat D1_inv = inv(D1);
  arma::mat D2_inv = inv(D2);
  arma::mat D12_inv_sum = D1_inv + D2_inv;

  arma::mat omega1 = arma::zeros<arma::mat>(m, p);

  for(int j = 0; j < m-1; ++j){
    omega1.row(j) = mu_y + estd;
  }

  arma::mat omega2 = arma::zeros<arma::mat>(n, p);
  for(int j = 0; j < n-1; ++j){
    omega2.row(j) = mu_x + estd;
  }

  arma::mat beta1 = inv(YY)*Y_t*omega1;
  arma::mat beta2 = inv(XX)*X_t*omega2;

  arma::mat estA = (beta1 + beta2)/2;
  arma::mat C1 = YY*estA*D1_inv/m + XX*estA*D2_inv/n;
  arma::mat C2 = Y_t*omega1*D1_inv/m + X_t*omega2*D2_inv/n;
  arma::mat A_gradient = C1 - C2;
  arma::mat pre = arma::ceil(log10(abs(A_gradient)));
  arma::rowvec pre2 = estd;
  arma::rowvec d_gradient = estd;
  double start1 = 0;
  double start2 = 0;
  for(int kk = 0; kk < steps-1; ++kk){

    for(int j = 0; j < m-1; ++j){
      omega1.row(j) = mu_y + estd;
    }

    for(int j = 0; j < n-1; ++j){
      omega2.row(j) = mu_x + estd;
    }

    C1 = YY*estA*D1_inv/m + XX*estA*D2_inv/n;
    C2 = Y_t*omega1*D1_inv/m + X_t*omega2*D2_inv/n;
    A_gradient = C1 - C2;

    if(kk==1){
      pre = arma::ceil(log10(abs(A_gradient)));
      start1 = pre.max();
      start1 = pow(10, -start1);
    //  Rcpp::Rcout << "start1" << start1 << std::endl;
    }

    if(start1 >= 0.01){
      estA = estA - start1*A_gradient/(kk+1);
    } else {
      estA = estA - start1*A_gradient;
    }


    arma::mat XestA = X*estA;
    arma::mat YestA = Y*estA;
    arma::rowvec errorX = arma::zeros<arma::rowvec>(p);
    arma::rowvec errorY = arma::zeros<arma::rowvec>(p);
    for(int j = 0; j < n-1; ++j){
      errorX = errorX + XestA.row(j);
    }
    for(int j = 0; j < m-1; ++j){
      errorY = errorY + YestA.row(j);
    }
    arma::rowvec ave_errorX = errorX/n - mu_x;
    arma::rowvec ave_errorY = errorY/m - mu_y;
    arma::rowvec dmat_part = ave_errorX*D2_inv + ave_errorY*D1_inv;
    d_gradient = estd*D12_inv_sum - dmat_part;

    if(kk==1){
      pre2 = arma::ceil(log10(abs(d_gradient)));
      start2 = pre2.max();
      start2 = pow(10, -start2);
   //   Rcpp::Rcout << "start2" << start2 << std::endl;
    }

    if(start1 >= 0.01){
      estd = estd - start2*d_gradient/(kk+1);
    } else {
      estd = estd - start2*d_gradient;
    }


  }

 // Rcpp::Rcout << "d" << d_gradient  << std::endl;
 // Rcpp::Rcout << "A" << A_gradient  << std::endl;

  return Rcpp::List::create(Rcpp::Named("A") = estA,
                            Rcpp::Named("d") = estd);

}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List run_simple_est_coord_three_clusters_cpp(arma::mat X, arma::mat Y, arma::mat Z, arma::mat x1, arma::mat y1, arma::mat z1, int steps = 1000){

  uint32_t n = X.n_rows;
  uint32_t m = Y.n_rows;
  uint32_t l = Z.n_rows;
  uint32_t p = X.n_cols;

  arma::rowvec mu_x = mean(x1, 0);
  arma::rowvec mu_y = mean(y1, 0);
  arma::rowvec mu_z = mean(z1, 0);
  arma::mat D1 = cov(y1);
  arma::mat D2 = cov(x1);
  arma::mat D3 = cov(z1);

  arma::rowvec estd = (mean(X, 0) - mu_x + mean(Y, 0) - mu_y + mean(Z, 0) - mu_z)/2;

  arma::mat X_t = X.t();
  arma::mat Y_t = Y.t();
  arma::mat Z_t = Z.t();
  arma::mat XX = X_t*X;
  arma::mat YY = Y_t*Y;
  arma::mat ZZ = Z_t*Z;

  arma::mat D1_inv = inv(D1);
  arma::mat D2_inv = inv(D2);
  arma::mat D3_inv = inv(D3);
  arma::mat D123_inv_sum = D1_inv + D2_inv + D3_inv;

  arma::mat omega1 = arma::zeros<arma::mat>(m, p);

  for(int j = 0; j < m-1; ++j){
    omega1.row(j) = mu_y + estd;
  }

  arma::mat omega2 = arma::zeros<arma::mat>(n, p);
  for(int j = 0; j < n-1; ++j){
    omega2.row(j) = mu_x + estd;
  }

  arma::mat omega3 = arma::zeros<arma::mat>(l, p);
  for(int j = 0; j < n-1; ++j){
    omega3.row(j) = mu_z + estd;
  }

  arma::mat beta1 = inv(YY)*Y_t*omega1;
  arma::mat beta2 = inv(XX)*X_t*omega2;
  arma::mat beta3 = inv(ZZ)*Z_t*omega2;

  arma::mat estA = (beta1 + beta2 + beta3)/2;
  arma::mat C1 = YY*estA*D1_inv/m + XX*estA*D2_inv/n + ZZ*estA*D3_inv/l;
  arma::mat C2 = Y_t*omega1*D1_inv/m + X_t*omega2*D2_inv/n + Z_t*omega3*D3_inv/n;
  arma::mat A_gradient = C1 - C2;
  arma::mat pre = arma::ceil(log10(abs(A_gradient)));
  arma::rowvec pre2 = estd;
  arma::rowvec d_gradient = estd;
  double start1 = 0;
  double start2 = 0;
  for(int kk = 0; kk < steps-1; ++kk){

    for(int j = 0; j < m-1; ++j){
      omega1.row(j) = mu_y + estd;
    }

    for(int j = 0; j < n-1; ++j){
      omega2.row(j) = mu_x + estd;
    }

    for(int j = 0; j < l-1; ++j){
      omega3.row(j) = mu_z + estd;
    }

    C1 = YY*estA*D1_inv/m + XX*estA*D2_inv/n + ZZ*estA*D3_inv/l;
    C2 = Y_t*omega1*D1_inv/m + X_t*omega2*D2_inv/n + Z_t*omega3*D3_inv/n;
    A_gradient = C1 - C2;

    if(kk==1){
      pre = arma::ceil(log10(abs(A_gradient)));
      start1 = pre.max();
      start1 = pow(10, -start1);
      //  Rcpp::Rcout << "start1" << start1 << std::endl;
    }

    if(start1 >= 0.01){
      estA = estA - start1*A_gradient/(kk+1);
    } else {
      estA = estA - start1*A_gradient;
    }


    arma::mat XestA = X*estA;
    arma::mat YestA = Y*estA;
    arma::mat ZestA = Z*estA;
    arma::rowvec errorX = arma::zeros<arma::rowvec>(p);
    arma::rowvec errorY = arma::zeros<arma::rowvec>(p);
    arma::rowvec errorZ = arma::zeros<arma::rowvec>(p);
    for(int j = 0; j < n-1; ++j){
      errorX = errorX + XestA.row(j);
    }
    for(int j = 0; j < m-1; ++j){
      errorY = errorY + YestA.row(j);
    }
    for(int j = 0; j < m-1; ++j){
      errorZ = errorZ + ZestA.row(j);
    }

    arma::rowvec ave_errorX = errorX/n - mu_x;
    arma::rowvec ave_errorY = errorY/m - mu_y;
    arma::rowvec ave_errorZ = errorZ/m - mu_z;
    arma::rowvec dmat_part = ave_errorX*D2_inv + ave_errorY*D1_inv + ave_errorZ*D3_inv;
    d_gradient = estd*D123_inv_sum - dmat_part;

    if(kk==1){
      pre2 = arma::ceil(log10(abs(d_gradient)));
      start2 = pre2.max();
      start2 = pow(10, -start2);
      //   Rcpp::Rcout << "start2" << start2 << std::endl;
    }

    if(start1 >= 0.01){
      estd = estd - start2*d_gradient/(kk+1);
    } else {
      estd = estd - start2*d_gradient;
    }

  }

  // Rcpp::Rcout << "d" << d_gradient  << std::endl;
  // Rcpp::Rcout << "A" << A_gradient  << std::endl;

  return Rcpp::List::create(Rcpp::Named("A") = estA,
                            Rcpp::Named("d") = estd);

}
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List run_simple_t_distribution_est_cpp(arma::mat X, arma::mat Y, arma::mat x1, arma::mat y1,
                                             int steps = 1000, double nu = 5){


  double nu_factor = nu/(nu-2);
  uint32_t n = X.n_rows;
  uint32_t m = Y.n_rows;
  uint32_t p = X.n_cols;

  arma::rowvec mu_x = mean(x1, 0);
  arma::rowvec mu_y = mean(y1, 0);

//  Rcpp::Rcout << "0  "    << std::endl;
  arma::mat Sigma1 = callmtCov(y1, nu)/nu_factor;
  arma::mat Sigma2 = callmtCov(y1, nu)/nu_factor;


  arma::rowvec estd = (mean(X, 0) - mu_x + mean(Y, 0) - mu_y)/2;

  arma::mat X_t = X.t();
  arma::mat Y_t = Y.t();

  arma::mat Sigma1_inv = inv(Sigma1);
  arma::mat Sigma2_inv = inv(Sigma2);

  arma::mat estA = arma::zeros<arma::mat>(p, p);
  estA.eye();
  arma::mat A = estA.t();
  arma::mat pre = estA;
  arma::rowvec pre2 = estd;
  double start = 0;
  arma::mat A_gradient = arma::zeros<arma::mat>(2, 2);
  arma::rowvec d_gradient = arma::zeros<arma::rowvec>(2);
  for(int kk = 0; kk < steps-1; ++kk){

    Rcpp::Rcout << "estd" << estd << std::endl;
    Rcpp::Rcout << "estA" <<  estA << std::endl;

    A = estA.t();

    arma::rowvec Omega1 = mu_y + estd;
    arma::rowvec Omega2 = mu_x + estd;


    arma::mat A_gradient_part1 = arma::zeros<arma::mat>(2, 2);
    arma::mat A_gradient_part2 = arma::zeros<arma::mat>(2, 2);

    arma::rowvec d_gradient_part1 = arma::zeros<arma::rowvec>(2);
    arma::rowvec d_gradient_part2 = arma::zeros<arma::rowvec>(2);

    for(int j = 0; j < m-1; ++j){
      arma::mat ay_d_u = A*Y_t.col(j) - Omega1.t();
      arma::mat S_ay_d_u = Sigma1_inv*ay_d_u;
      arma::mat C_1u = Y_t.col(j)*Y.row(j)*estA*Sigma1_inv - Y_t.col(j)*Omega1*Sigma1_inv;
      arma::mat C_1d = ay_d_u.t()*S_ay_d_u;
      A_gradient_part1 = A_gradient_part1 + C_1u/(1+C_1d(0, 0));
      arma::mat S_ay_d_u_t = S_ay_d_u.t();
      d_gradient_part1 = d_gradient_part1 - S_ay_d_u_t/(1+C_1d(0, 0));
    }

    for(int j = 0; j < n-1; ++j){
      arma::mat ax_d_u = A*X_t.col(j)  - Omega2.t();
      arma::mat S_ax_d_u = Sigma2_inv*ax_d_u;
      arma::mat C_2u = X_t.col(j)*X.row(j)*estA*Sigma2_inv - X_t.col(j)*Omega2*Sigma2_inv;
      arma::mat C_2d = ax_d_u.t()*S_ax_d_u;
      A_gradient_part2 = A_gradient_part2 + C_2u/(1+C_2d(0, 0));
      arma::mat S_ax_d_u_t = S_ax_d_u.t();
      d_gradient_part2 = d_gradient_part2 - S_ax_d_u_t/(1+C_2d(0, 0));
    }


    A_gradient = A_gradient_part1/m + A_gradient_part2/n;
    d_gradient = d_gradient_part1/m + d_gradient_part2/n;

    Rcpp::Rcout << "d" << d_gradient << std::endl;
    Rcpp::Rcout << "A" << A_gradient << std::endl;

    if(kk==0){
      pre = arma::ceil(log10(abs(A_gradient)));
      pre2 = arma::ceil(log10(abs(d_gradient)));
      if(pre2.max() > pre.max()){
        start = pre2.max();
      } else {
        start = pre.max();
      }
      start = pow(10, -start);
      Rcpp::Rcout << "start" << start << std::endl;
    }
    // pre = arma::ceil(log10(abs(A_gradient)));

    //  pre2 = arma::ceil(log10(abs(d_gradient)));
    if(start >= 0.01){
      estA = estA - start*A_gradient/(kk+1);
      estd = estd - start*d_gradient/(kk+1);
    } else {
      estA = estA - start*A_gradient;
      estd = estd - start*d_gradient;
    }
  }
  //Rcpp::Rcout << "d" << d_gradient << std::endl;
  //Rcpp::Rcout << "A" << A_gradient << std::endl;

  return Rcpp::List::create(Rcpp::Named("A") = estA,
                            Rcpp::Named("d") = estd);

}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List run_simple_t_distribution_est_theta_cpp(arma::mat X, arma::mat Y, arma::mat x1, arma::mat y1,
                                             int steps = 1000, double nu = 5){


  double nu_factor = nu/(nu-2);
  uint32_t n = X.n_rows;
  uint32_t m = Y.n_rows;
  uint32_t p = X.n_cols;

  arma::rowvec mu_x = mean(x1, 0);
  arma::rowvec mu_y = mean(y1, 0);

  //  Rcpp::Rcout << "0  "    << std::endl;
  arma::mat Sigma1 = callmtCov(y1, nu)/nu_factor;
  arma::mat Sigma2 = callmtCov(y1, nu)/nu_factor;


  arma::rowvec estd = -(mean(X, 0) - mu_x + mean(Y, 0) - mu_y)/2;

  arma::mat X_t = X.t();
  arma::mat Y_t = Y.t();

  arma::mat Sigma1_inv = inv(Sigma1);
  arma::mat Sigma2_inv = inv(Sigma2);

  double estTheta = 0;
  arma::mat estA = rotation_matrix(estTheta);
  arma::mat A_T = estA.t();
  arma::mat A_dev = rotation_matrix_dev(estTheta);
  arma::mat A_dev_T = A_dev.t();
  double start = 0;

  arma::mat Theta_gradient = arma::zeros<arma::mat>(1, 1);
  arma::rowvec d_gradient = arma::zeros<arma::rowvec>(p);
  arma::mat pre = Theta_gradient;
  arma::rowvec pre2 = d_gradient;

  for(int kk = 0; kk < steps-1; ++kk){

    Rcpp::Rcout << "estd " << estd << std::endl;
    Rcpp::Rcout << "estTheta " << estTheta << std::endl;

    A_T = estA.t();
    A_dev = rotation_matrix_dev(estTheta);
    A_dev_T = A_dev.t();

    arma::mat Theta_gradient_part1 = arma::zeros<arma::mat>(1, 1);
    arma::mat Theta_gradient_part2 = arma::zeros<arma::mat>(1, 1);

    arma::rowvec d_gradient_part1 = arma::zeros<arma::rowvec>(p);
    arma::rowvec d_gradient_part2 = arma::zeros<arma::rowvec>(p);

    arma::mat omega1 = A_T*Sigma1_inv*A_dev;
    arma::mat omega2 = A_T*Sigma2_inv*A_dev;

    arma::rowvec mu_y_est_d = mu_y + estd;
    arma::rowvec mu_x_est_d = mu_x + estd;
    arma::mat omega3 = A_dev_T*mu_y_est_d.t();
    arma::mat omega4 = A_dev_T*mu_x_est_d.t();

    for(int j = 0; j < m-1; ++j){
      arma::mat ay_d_u = A_T*Y_t.col(j) - mu_y_est_d.t();
      arma::mat S_ay_d_u = Sigma1_inv*ay_d_u;
      arma::mat C_1u = Y.row(j)*omega1*Y_t.col(j) - Y.row(j)*omega3;
      arma::mat C_1d = ay_d_u.t()*S_ay_d_u;
      Theta_gradient_part1 = Theta_gradient_part1 + C_1u/(1+C_1d(0, 0));
      arma::mat S_ay_d_u_t = S_ay_d_u.t();
      d_gradient_part1 = d_gradient_part1 - S_ay_d_u_t/(1+C_1d(0, 0));
    }

    for(int j = 0; j < n-1; ++j){
      arma::mat ax_d_u = A_T*X_t.col(j)  - mu_x_est_d.t();
      arma::mat S_ax_d_u = Sigma2_inv*ax_d_u;
      arma::mat C_2u = X.row(j)*omega2*X_t.col(j) - X.row(j)*omega4;
      arma::mat C_2d = ax_d_u.t()*S_ax_d_u;
      Theta_gradient_part2 = Theta_gradient_part2 + C_2u/(1+C_2d(0, 0));
      arma::mat S_ax_d_u_t = S_ax_d_u.t();
      d_gradient_part2 = d_gradient_part2 - S_ax_d_u_t/(1+C_2d(0, 0));
    }

    Theta_gradient = Theta_gradient_part1/m + Theta_gradient_part2/n;
    d_gradient = d_gradient_part1/m + d_gradient_part2/n;

    //if(kk==0){
      pre = arma::ceil(log10(abs(Theta_gradient)));
      pre2 = arma::ceil(log10(abs(d_gradient)));
      if(pre2.max() > pre.max()){
        start = pre2.max();
      } else {
        start = pre.max();
      }
      start = pow(10, -start);
    //  Rcpp::Rcout << "start" << start << std::endl;
   // }
    Rcpp::Rcout << "d" << d_gradient << std::endl;
    Rcpp::Rcout << "Theta" << Theta_gradient << std::endl;
    // pre = arma::ceil(log10(abs(A_gradient)));

    //  pre2 = arma::ceil(log10(abs(d_gradient)));

//  if(start >= 0.01){
//      estTheta = estTheta - start*Theta_gradient(0, 0)/(kk+1);
//      estd = estd - start*d_gradient/(kk+1);
//    } else {

      estTheta = estTheta - start*Theta_gradient(0, 0)/(kk+1);
      estd = estd - start*d_gradient/(kk+1);
//    }

    estA = rotation_matrix(estTheta);
  }


  return Rcpp::List::create(Rcpp::Named("Theta") = estTheta,
                            Rcpp::Named("A") = estA,
                            Rcpp::Named("d") = estd);

}











// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List run_simple_angle_t_distribution_est_cpp(arma::mat X, arma::mat Y, arma::mat x1, arma::mat y1, arma::vec angle_list, int steps = 1000, double nu = 5){

  double nu_factor = nu/(nu-2);
  uint32_t n = X.n_rows;
  uint32_t m = Y.n_rows;
  uint32_t p = X.n_cols;
  uint32_t al = angle_list.n_elem;

  arma::rowvec mu_x = mean(x1, 0);
  arma::rowvec mu_y = mean(y1, 0);

  arma::mat Sigma1 = callmtCov(y1, nu)/nu_factor;
  arma::mat Sigma2 = callmtCov(y1, nu)/nu_factor;

  arma::rowvec estd = (mean(X, 0) - mu_x + mean(Y, 0) - mu_y)/2;

  arma::mat X_t = X.t();
  arma::mat Y_t = Y.t();
  arma::mat XX = X_t*X;
  arma::mat YY = Y.t()*Y;

  arma::mat Sigma1_inv = inv(Sigma1);
  arma::mat Sigma2_inv = inv(Sigma2);
  arma::mat Sigma12_inv_sum = Sigma1_inv + Sigma2_inv;

  arma::mat omega1 = arma::zeros<arma::mat>(m, p);

  for(int j = 0; j < m-1; ++j){
    omega1.row(j) = mu_y + estd;
  }

  arma::mat omega2 = arma::zeros<arma::mat>(n, p);
  for(int j = 0; j < n-1; ++j){
    omega2.row(j) = mu_x + estd;
  }

  arma::mat iniA = rotation_matrix(0);
  arma::mat estA = iniA;
  arma::mat A = estA.t();
  arma::vec changes = arma::zeros<arma::vec>(al);
  double selected = 0;

  for(int kk = 0; kk < steps-1; ++kk){

    A = estA.t();

    for(int j = 0; j < m-1; ++j){
      omega1.row(j) = mu_y + estd;
    }

    for(int j = 0; j < n-1; ++j){
      omega2.row(j) = mu_x + estd;
    }

    arma::mat A_gradient_part2 = arma::zeros<arma::mat>(2, 2);
    arma::mat A_gradient_part1 = arma::zeros<arma::mat>(2, 2);
    arma::rowvec d_gradient_part1 = arma::zeros<arma::vec>(2);
    arma::rowvec d_gradient_part2 = arma::zeros<arma::vec>(2);

    for(int j = 0; j < m-1; ++j){
      arma::mat ay_d_u = A*Y.row(j) - omega1;
      arma::mat S_ay_d_u = Sigma1_inv*ay_d_u;
      arma::mat Yrowj = Y.row(j);
      arma::mat C_1u = Y.row(j)*Yrowj.t()*estA*Sigma1_inv - Y.row(j)*omega1.t()*Sigma1_inv;
      arma::mat C_1d = ay_d_u.t()*S_ay_d_u;
      A_gradient_part1 = A_gradient_part1 + C_1u/(1+C_1d(0, 0));
      d_gradient_part1 = d_gradient_part1 - S_ay_d_u/(1+C_1d(0, 0));
    }

    for(int j = 0; j < n-1; ++j){
      arma::mat ax_d_u = A*X.row(j) - omega2;
      arma::mat S_ax_d_u = Sigma2_inv*ax_d_u;
      arma::mat Xrowj = X.row(j);
      arma::mat C_2u = X.row(j)*Xrowj.t()*estA*Sigma2_inv - X.row(j)*omega2.t()*Sigma2_inv;
      arma::mat C_2d = ax_d_u.t()*S_ax_d_u;
      A_gradient_part2 = A_gradient_part2 + C_2u/(1+C_2d(0, 0));
      d_gradient_part2 = d_gradient_part2 - S_ax_d_u/(1+C_2d(0, 0));
    }
    arma::mat A_gradient = A_gradient_part1/m + A_gradient_part2/n;
    arma::mat d_gradient = d_gradient_part1/m + d_gradient_part2/n;

    for(int i = 0; i < al-1; ++i){
      arma::mat A_prime_gradient_part2 = arma::zeros<arma::mat>(2, 2);
      arma::mat A_prime_gradient_part1 = arma::zeros<arma::mat>(2, 2);
      arma::mat estA_prime = rotation_matrix(angle_list(i));
      arma::mat A_prime = estA_prime.t();

      for(int j = 0; j < m-1; ++j){
        arma::mat ay_d_u = A_prime*Y.row(j) - omega1;
        arma::mat S_ay_d_u = Sigma1_inv*ay_d_u;
        arma::mat Yrowj = Y.row(j);
        arma::mat C_1u = Y.row(j)*Yrowj.t()*estA_prime*Sigma1_inv - Y.row(j)*omega1.t()*Sigma1_inv;
        arma::mat C_1d = ay_d_u.t()*S_ay_d_u;
        A_prime_gradient_part1 = A_prime_gradient_part1 + C_1u/(1+C_1d(0, 0));
      }
      for(int j = 0; j < n-1; ++j){
        arma::mat ax_d_u = A_prime*X.row(j) - omega2;
        arma::mat S_ax_d_u = Sigma2_inv*ax_d_u;
        arma::mat Xrowj = X.row(j);
        arma::mat C_2u = X.row(j)*Xrowj.t()*estA_prime*Sigma2_inv - X.row(j)*omega2.t()*Sigma2_inv;
        arma::mat C_2d = ax_d_u.t()*S_ax_d_u;
        A_prime_gradient_part2 = A_prime_gradient_part2 + C_2u/(1+C_2d(0, 0));
      }
      arma::mat A_prime_gradient = A_prime_gradient_part1/m + A_prime_gradient_part2/n;

      arma::mat med = median(abs(A_prime_gradient)/abs(A_gradient));
      changes(i) = med(0, 0);
    }

    arma::uvec aa_ind = arma::find(changes < 1);
    if(aa_ind.n_elem >= 1){
      selected = angle_list(aa_ind(0));
      estA = rotation_matrix(selected);
    } else {
      estA = iniA;
    }

    estd = estd - 0.01*d_gradient;

  }

  return Rcpp::List::create(Rcpp::Named("A") = estA,
                            Rcpp::Named("A") = estA,
                            Rcpp::Named("d") = estd);

}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List run_simple_angle_t_distribution_likelihood_est_cpp(arma::mat X, arma::mat Y, arma::mat x1, arma::mat y1, arma::vec angle_list, int steps = 1000, double nu = 5){


  double nu_factor = nu/(nu-2);
  uint32_t n = X.n_rows;
  uint32_t m = Y.n_rows;
  uint32_t p = X.n_cols;
  uint32_t al = angle_list.n_elem;

  arma::rowvec mu_x = mean(x1, 0);
  arma::rowvec mu_y = mean(y1, 0);

  arma::mat Sigma1 = callmtCov(y1, nu)/nu_factor;
  arma::mat Sigma2 = callmtCov(y1, nu)/nu_factor;

  arma::rowvec estd = (mean(X, 0) - mu_x + mean(Y, 0) - mu_y)/2;

  arma::mat X_t = X.t();
  arma::mat Y_t = Y.t();
  arma::mat XX = X_t*X;
  arma::mat YY = Y.t()*Y;

  arma::mat Sigma1_inv = inv(Sigma1);
  arma::mat Sigma2_inv = inv(Sigma2);
  arma::mat Sigma12_inv_sum = Sigma1_inv + Sigma2_inv;

  arma::mat omega1 = arma::zeros<arma::mat>(m, p);

  for(int j = 0; j < m-1; ++j){
    omega1.row(j) = mu_y + estd;
  }

  arma::mat omega2 = arma::zeros<arma::mat>(n, p);
  for(int j = 0; j < n-1; ++j){
    omega2.row(j) = mu_x + estd;
  }

  arma::mat iniA = rotation_matrix(0);
  arma::mat estA = iniA;
  arma::mat A = estA.t();
  arma::vec changes = arma::zeros<arma::vec>(al);
  double selected = 0;

  for(int kk = 0; kk < steps-1; ++kk){

    A = estA.t();

    for(int j = 0; j < m-1; ++j){
      omega1.row(j) = mu_y + estd;
    }

    for(int j = 0; j < n-1; ++j){
      omega2.row(j) = mu_x + estd;
    }

    arma::mat A_gradient_part2 = arma::zeros<arma::mat>(2, 2);
    arma::mat A_gradient_part1 = arma::zeros<arma::mat>(2, 2);
    arma::rowvec d_gradient_part1 = arma::zeros<arma::vec>(2);
    arma::rowvec d_gradient_part2 = arma::zeros<arma::vec>(2);

    for(int j = 0; j < m-1; ++j){
      arma::mat ay_d_u = A*Y.row(j) - omega1;
      arma::mat S_ay_d_u = Sigma1_inv*ay_d_u;
      arma::mat Yrowj = Y.row(j);
      arma::mat C_1u = Y.row(j)*Yrowj.t()*estA*Sigma1_inv - Y.row(j)*omega1.t()*Sigma1_inv;
      arma::mat C_1d = ay_d_u.t()*S_ay_d_u;
      A_gradient_part1 = A_gradient_part1 + C_1u/(1+C_1d(0, 0));
      d_gradient_part1 = d_gradient_part1 - S_ay_d_u/(1+C_1d(0, 0));
    }

    for(int j = 0; j < n-1; ++j){
      arma::mat ax_d_u = A*X.row(j) - omega2;
      arma::mat S_ax_d_u = Sigma2_inv*ax_d_u;
      arma::mat Xrowj = X.row(j);
      arma::mat C_2u = X.row(j)*Xrowj.t()*estA*Sigma2_inv - X.row(j)*omega2.t()*Sigma2_inv;
      arma::mat C_2d = ax_d_u.t()*S_ax_d_u;
      A_gradient_part2 = A_gradient_part2 + C_2u/(1+C_2d(0, 0));
      d_gradient_part2 = d_gradient_part2 - S_ax_d_u/(1+C_2d(0, 0));
    }
    arma::mat A_gradient = A_gradient_part1/m + A_gradient_part2/n;
    arma::mat d_gradient = d_gradient_part1/m + d_gradient_part2/n;

    for(int i = 0; i < al-1; ++i){
      arma::mat estA_prime = rotation_matrix(angle_list(i));
      changes(i) = likelihood_para_t_cpp(X, Y, mu_x, mu_y, estA_prime, estd, Sigma1, Sigma2, nu, true);
    }

    arma::uword aa_ind = index_max(changes);
    selected = angle_list(aa_ind);
    estA = rotation_matrix(selected);

    estd = estd - 0.01*d_gradient;

  }

  return Rcpp::List::create(Rcpp::Named("A") = estA,
                            Rcpp::Named("d") = estd);

}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat gradient_A(arma::mat XX, arma::mat YY, arma::mat X_t, arma::mat Y_t, arma::rowvec mu_x, arma::rowvec mu_y, arma::mat estA, arma::rowvec estd,
                     uint32_t n, uint32_t m, uint32_t p, arma::mat D1_inv, arma::mat D2_inv){

  arma::mat omega1 = arma::zeros<arma::mat>(m, p);
  arma::mat omega2 = arma::zeros<arma::mat>(n, p);

  for(int j = 0; j < m-1; ++j){
    omega1.row(j) = mu_y + estd;
  }

  for(int j = 0; j < n-1; ++j){
    omega2.row(j) = mu_x + estd;
  }

  arma::mat C1 = YY*estA*D1_inv/m + XX*estA*D2_inv/n;
  arma::mat C2 = Y_t*omega1*D1_inv/m + X_t*omega2*D2_inv/n;
  arma::mat A_gradient = C1 - C2;

  return A_gradient;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::rowvec gradient_d(arma::mat X, arma::mat Y, arma::rowvec mu_x, arma::rowvec mu_y, arma::mat estA, arma::rowvec estd,
                     uint32_t n, uint32_t m, uint32_t p, arma::mat D1_inv, arma::mat D2_inv, arma::mat D12_inv_sum){

  arma::mat XestA = X*estA;
  arma::mat YestA = Y*estA;
  arma::rowvec errorX = arma::zeros<arma::rowvec>(p);
  arma::rowvec errorY = arma::zeros<arma::rowvec>(p);
  for(int j = 0; j < n-1; ++j){
    errorX = errorX + XestA.row(j);
  }
  for(int j = 0; j < m-1; ++j){
    errorY = errorY + YestA.row(j);
  }
  arma::rowvec ave_errorX = errorX/n - mu_x;
  arma::rowvec ave_errorY = errorY/m - mu_y;
  arma::rowvec dmat_part = ave_errorX*D2_inv + ave_errorY*D1_inv;
  arma::rowvec d_gradient = estd*D12_inv_sum - dmat_part;

  return d_gradient;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat select_stepsize_for_A(arma::mat X, arma::mat Y, arma::rowvec mu_x, arma::rowvec mu_y, arma::mat estA, arma::rowvec estd, double ll,
                                   arma::mat A_gradient, uint32_t p, arma::mat D1, arma::mat D2, double gamma, double down){

  arma::mat A_prime = estA;
  arma::mat selected = estA;

  for(int i = 0; i < p; ++i){
    for(int j = 0; j < p; ++j){

      double gra_now = A_gradient(i, j);
      double gra_now2 = gra_now*gra_now*gamma;
      double start = sqrt(abs(estA(i, j)/gra_now))/2;
      double aa = start;

      while(aa > 0){
        double aa2 = aa*aa;
        A_prime(i, j) = estA(i, j) - aa2*gra_now;
        double lA_prime = likelihood_para_cpp(X, Y, mu_x, mu_y, A_prime, estd, D1, D2, true);
        if(lA_prime - ll + aa2*gra_now2 > 0){
          selected(i, j)= A_prime(i, j);
          Rcpp::Rcout << "  lA_prime" <<   lA_prime << std::endl;
          break;
        }
        aa = aa - start*down;
      }
   //   A_prime = estA;
    }
  }
  return selected;
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::rowvec select_stepsize_for_d(arma::mat X, arma::mat Y, arma::rowvec mu_x, arma::rowvec mu_y, arma::mat estA, arma::rowvec estd, double ll,
                             arma::rowvec d_gradient, uint32_t p, arma::mat D1, arma::mat D2, double gamma, double down){

  arma::rowvec d_prime = estd;
  arma::rowvec selected = estd;

  for(int j = 0; j < p; ++j){

    double gra_now = d_gradient(j);
    double gra_now2 = gra_now*gra_now*gamma;
    double start = sqrt(abs(estd(j)/gra_now))/2;
    double aa = start;

    while(aa > 0){
      double aa2 = aa*aa;
      d_prime(j) = estd(j) - aa2*gra_now;
      double ld_prime = likelihood_para_cpp(X, Y, mu_x, mu_y, estA, d_prime, D1, D2, true);
      if(ld_prime - ll + aa2*gra_now2 > 0){
        selected(j) = d_prime(j);
        Rcpp::Rcout << "  ld_prime" <<   ld_prime << std::endl;
        break;
      }
      aa = aa - start*down;
    }
 //   d_prime = estd;
  }

  return selected;
}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat select_stepsize_for_A_alt(arma::mat X, arma::mat Y, arma::rowvec mu_x, arma::rowvec mu_y, arma::mat estA, arma::rowvec estd, double ll,
                                arma::mat A_gradient, uint32_t p, arma::mat D1, arma::mat D2, double beta){

  arma::mat A_prime = estA;
  arma::mat selected = estA;

  for(int i = 0; i < p; ++i){
    for(int j = 0; j < p; ++j){

      double gra_now = A_gradient(i, j);
      double gra_now2 = gra_now*gra_now;
      double t = 1;
      A_prime(i, j) = estA(i, j) - t*gra_now;
      double ld_prime = likelihood_para_cpp(X, Y, mu_x, mu_y, A_prime, estd,  D1, D2, true);
      while(ld_prime - ll + t*gra_now2/2 > 0){
        t *= beta;
        A_prime(i, j) = estA(i, j) - t*gra_now;
        ld_prime = likelihood_para_cpp(X, Y, mu_x, mu_y, A_prime, estd,  D1, D2, true);
      }
      selected(i, j) = estA(i, j) - t*gra_now;
      A_prime = estA;
    }
  }
  return selected;
}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::rowvec select_stepsize_for_d_alt(arma::mat X, arma::mat Y, arma::rowvec mu_x, arma::rowvec mu_y, arma::mat estA, arma::rowvec estd, double ll,
                                   arma::rowvec d_gradient, uint32_t p, arma::mat D1, arma::mat D2, double beta){

  arma::rowvec d_prime = estd;
  arma::rowvec selected = estd;

  for(int j = 0; j < p; ++j){

    double gra_now = d_gradient(j);
    double gra_now2 = gra_now*gra_now;
    double t = 1;

    d_prime(j) = estd(j) - t*gra_now;
    double ld_prime = likelihood_para_cpp(X, Y, mu_x, mu_y, estA, d_prime, D1, D2, true);
    while(ld_prime - ll + t*gra_now2/2 < 0){
      t *= beta;
      d_prime(j) = estd(j) - t*gra_now;
      ld_prime = likelihood_para_cpp(X, Y, mu_x, mu_y, estA, d_prime, D1, D2, true);
    }
    selected(j) = estd(j) - t*gra_now;
    d_prime = estd;

  }

  return selected;
}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List run_est_with_gradient_descent_cpp(arma::mat X, arma::mat Y, arma::mat x1, arma::mat y1,
                                             int steps,  double gamma, double down){ // double beta) {//

  uint32_t n = X.n_rows;
  uint32_t m = Y.n_rows;
  uint32_t p = X.n_cols;

  arma::rowvec mu_x = mean(x1, 0);
  arma::rowvec mu_y = mean(y1, 0);
  arma::mat D1 = cov(y1);
  arma::mat D2 = cov(x1);

  arma::rowvec estd = (mean(X, 0) - mu_x + mean(Y, 0) - mu_y)/2;

  arma::mat X_t = X.t();
  arma::mat Y_t = Y.t();
  arma::mat XX = X_t*X;
  arma::mat YY = Y_t*Y;

  arma::mat D1_inv = inv(D1);
  arma::mat D2_inv = inv(D2);
  arma::mat D12_inv_sum = D1_inv + D2_inv;

  arma::mat estA = arma::zeros<arma::mat>(p, p);
  estA.eye();
  estA += 0.05;
  double likelihood_est = likelihood_para_cpp(X, Y, mu_x, mu_y, estA, estd, D1, D2, true);

  arma::mat A_gradient = gradient_A(XX, YY, X_t, Y_t, mu_x, mu_y, estA, estd, n, m, p, D1_inv, D2_inv);

  arma::rowvec d_gradient = gradient_d(X, Y, mu_x, mu_y, estA, estd, n, m, p, D1_inv, D2_inv, D12_inv_sum);

  Rcpp::Rcout << "  start " <<   likelihood_est  << std::endl;
  Rcpp::Rcout << "  start " <<   A_gradient  << std::endl;
  Rcpp::Rcout << "  start " <<   d_gradient  << std::endl;
  arma::mat A_prime = arma::zeros<arma::mat>(p, p);
  arma::rowvec d_prime = arma::zeros<arma::rowvec>(p);

  for(int kk = 0; kk < steps-1; ++kk){

    arma::mat absd = abs(d_gradient);
    double val = absd.min();

    if( val >= 0.0001){
      d_prime = select_stepsize_for_d(X, Y, mu_x, mu_y, estA, estd, likelihood_est, d_gradient, p, D1, D2, gamma, down);
   // d_prime = select_stepsize_for_d_alt(X, Y, mu_x, mu_y, estA, estd, likelihood_est, d_gradient, p, D1, D2, beta);
    } else {
      d_prime = estd;
    }

    arma::mat absA = abs(A_gradient);
    val = absA.min();

    if( val >= 0.0001){
      A_prime = select_stepsize_for_A(X, Y, mu_x, mu_y, estA, estd, likelihood_est, A_gradient, p, D1, D2, gamma, down);
   //  A_prime = select_stepsize_for_A_alt(X, Y, mu_x, mu_y, estA, estd, likelihood_est, A_gradient, p, D1, D2, beta);
    } else {
      A_prime = estA;
    }

    estA = A_prime; estd = d_prime;

    A_gradient = gradient_A(XX, YY, X_t, Y_t, mu_x, mu_y, estA, estd, n, m, p, D1_inv, D2_inv);
    Rcpp::Rcout << "  A_gradient " <<   A_gradient  << std::endl;

    d_gradient = gradient_d(X, Y, mu_x, mu_y, estA, estd, n, m, p, D1_inv, D2_inv, D12_inv_sum);
    Rcpp::Rcout << "  d_gradient " <<   d_gradient  << std::endl;

    likelihood_est = likelihood_para_cpp(X, Y, mu_x, mu_y, estA, estd, D1, D2, true);
    Rcpp::Rcout << "  likelihood_est " <<   likelihood_est  << std::endl;

  }

  return Rcpp::List::create(Rcpp::Named("A") = estA,
                            Rcpp::Named("d") = estd);

}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double gradient_theta(arma::mat X, arma::mat Y, arma::rowvec mu_x, arma::rowvec mu_y, double estTheta, arma::mat estA, arma::rowvec estd,
                      uint32_t n, uint32_t m, uint32_t p, arma::mat D1_inv, arma::mat D2_inv){

  arma::mat A_T = estA.t();
  arma::mat A_dev = rotation_matrix_dev(estTheta);
  arma::mat A_dev_T = A_dev.t();

  arma::mat omega1 = A_T*D1_inv*A_dev;
  arma::mat omega2 = A_T*D2_inv*A_dev;

  arma::mat omega3 = A_dev_T*(mu_y + estd);
  arma::mat omega4 = A_dev_T*(mu_x + estd);

  arma::mat Theta_gradient = arma::zeros<arma::mat>(1, 1);
  for(int j = 0; j < m-1; ++j){
    Theta_gradient +=  Y.row(0)*omega1*Y.col(0) - Y.row(0)*omega3;
  }

  arma::mat Theta_gradient2 = arma::zeros<arma::mat>(1, 1);
  for(int j = 0; j < n-1; ++j){
    Theta_gradient2 +=  X.row(0)*omega2*X.col(0) - X.row(0)*omega4;
  }

  double res = Theta_gradient(0, 0)/m + Theta_gradient2(0, 0)/n;
  return res;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double select_stepsize_for_theta(arma::mat X, arma::mat Y, arma::rowvec mu_x, arma::rowvec mu_y, double estTheta, arma::mat estA, arma::rowvec estd, double ll,
                                double Theta_gradient, uint32_t p, arma::mat D1, arma::mat D2, double gamma, double down){

  double Theta_prime = estTheta;
  double selected = estTheta;

  double gra_now = Theta_gradient;
  double gra_now2 = gra_now*gra_now*gamma;
  double start = sqrt(abs(estTheta/gra_now))/2;
  double aa = start;

  arma::mat A_prime = estA;
  while(aa > 0){
    double aa2 = aa*aa;
    Theta_prime = estTheta - aa2*gra_now;
    A_prime = rotation_matrix(Theta_prime);
    double lA_prime = likelihood_para_cpp(X, Y, mu_x, mu_y, A_prime, estd, D1, D2, true);
    if(lA_prime - ll - aa2*gra_now2 > 0){
        selected = Theta_prime;
        break;
      }
      aa = aa - start*down;
    }

  return selected;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List run_est_with_gradient_descent_angle_search_cpp(arma::mat X, arma::mat Y, arma::mat x1, arma::mat y1, int steps, double gamma, double down){

  uint32_t n = X.n_rows;
  uint32_t m = Y.n_rows;
  uint32_t p = X.n_cols;

  arma::rowvec mu_x = mean(x1, 0);
  arma::rowvec mu_y = mean(y1, 0);
  arma::mat D1 = cov(y1);
  arma::mat D2 = cov(x1);

  arma::rowvec estd = (mean(X, 0) - mu_x + mean(Y, 0) - mu_y)/2;

  arma::mat D1_inv = inv(D1);
  arma::mat D2_inv = inv(D2);
  arma::mat D12_inv_sum = D1_inv + D2_inv;

  double estTheta = 0;
  arma::mat estA = rotation_matrix(estTheta);

  double likelihood_est = likelihood_para_cpp(X, Y, mu_x, mu_y, estA, estd, D1, D2, true);

  double Theta_gradient = gradient_theta(X, Y, mu_x, mu_y, estTheta, estA, estd, n, m, p, D1_inv, D2_inv);

  arma::rowvec d_gradient = gradient_d(X, Y, mu_x, mu_y, estA, estd, n, m, p, D1_inv, D2_inv, D12_inv_sum);

  arma::mat A_prime = arma::zeros<arma::mat>(p, p);;
  arma::rowvec d_prime = arma::zeros<arma::rowvec>(p);
  double Theta_prime = 0;

  for(int kk = 0; kk < steps-1; ++kk){

    arma::mat absd = abs(d_gradient);
    double val = absd.min();

    if( val >= 0.0001){
      d_prime = select_stepsize_for_d(X, Y, mu_x, mu_y, estA, estd, likelihood_est, d_gradient, p, D1, D2, gamma, down);
    } else {
      d_prime = estd;
    }

    if( abs(Theta_gradient) >= 0.0001){
      Theta_prime = select_stepsize_for_theta(X, Y, mu_x, mu_y, estTheta, estA, estd, likelihood_est, Theta_gradient, p, D1, D2, gamma, down);
      A_prime = rotation_matrix(estTheta);
    } else {
      Theta_prime = estTheta;
      A_prime = estA;
    }

    estA = A_prime; estd = d_prime; estTheta = Theta_prime;

    Theta_gradient = gradient_theta(X, Y, mu_x, mu_y, estTheta, estA, estd, n, m, p, D1_inv, D2_inv);
    d_gradient = gradient_d(X, Y, mu_x, mu_y, estA, estd, n, m, p, D1_inv, D2_inv, D12_inv_sum);

    likelihood_est = likelihood_para_cpp(X, Y, mu_x, mu_y, estA, estd, D1, D2, true);

  }

  return Rcpp::List::create(Rcpp::Named("Theta") = estTheta,
                            Rcpp::Named("A") = estA,
                            Rcpp::Named("d") = estd);

}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat gradient_A_t_dis(arma::mat X, arma::mat Y, arma::mat Sigma1_inv, arma::mat Sigma2_inv, arma::rowvec mu_x, arma::rowvec mu_y, arma::mat estA, arma::rowvec estd,
                     uint32_t n, uint32_t m, uint32_t p){

  arma::mat omega1 = arma::zeros<arma::mat>(m, p);
  arma::mat omega2 = arma::zeros<arma::mat>(n, p);

  arma::mat A = estA.t();

  for(int j = 0; j < m-1; ++j){
    omega1.row(j) = mu_y + estd;
  }

  for(int j = 0; j < n-1; ++j){
    omega2.row(j) = mu_x + estd;
  }

  arma::mat A_gradient_part2 = arma::zeros<arma::mat>(2, 2);
  arma::mat A_gradient_part1 = arma::zeros<arma::mat>(2, 2);

  for(int j = 0; j < m-1; ++j){
    arma::mat ay_d_u = A*Y.row(j) - omega1;
    arma::mat S_ay_d_u = Sigma1_inv*ay_d_u;
    arma::mat Yrowj = Y.row(j);
    arma::mat C_1u = Y.row(j)*Yrowj.t()*estA*Sigma1_inv - Y.row(j)*omega1.t()*Sigma1_inv;
    arma::mat C_1d = ay_d_u.t()*S_ay_d_u;
    A_gradient_part1 = A_gradient_part1 + C_1u/(1+C_1d(0, 0));
  }

  for(int j = 0; j < n-1; ++j){
    arma::mat ax_d_u = A*X.row(j) - omega2;
    arma::mat S_ax_d_u = Sigma2_inv*ax_d_u;
    arma::mat Xrowj = X.row(j);
    arma::mat C_2u = X.row(j)*Xrowj.t()*estA*Sigma2_inv - X.row(j)*omega2.t()*Sigma2_inv;
    arma::mat C_2d = ax_d_u.t()*S_ax_d_u;
    A_gradient_part2 = A_gradient_part2 + C_2u/(1+C_2d(0, 0));
  }
  arma::mat A_gradient = A_gradient_part1/m + A_gradient_part2/n;

  return A_gradient;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::rowvec gradient_d_t_dis(arma::mat X, arma::mat Y, arma::mat Sigma1_inv, arma::mat Sigma2_inv, arma::rowvec mu_x, arma::rowvec mu_y, arma::mat estA, arma::rowvec estd,
                        uint32_t n, uint32_t m, uint32_t p){

  arma::mat omega1 = arma::zeros<arma::mat>(m, p);
  arma::mat omega2 = arma::zeros<arma::mat>(n, p);

  arma::mat A = estA.t();

  for(int j = 0; j < m-1; ++j){
    omega1.row(j) = mu_y + estd;
  }

  for(int j = 0; j < n-1; ++j){
    omega2.row(j) = mu_x + estd;
  }

  arma::rowvec d_gradient_part1 = arma::zeros<arma::vec>(2);
  arma::rowvec d_gradient_part2 = arma::zeros<arma::vec>(2);

  for(int j = 0; j < m-1; ++j){
    arma::mat ay_d_u = A*Y.row(j) - omega1;
    arma::mat S_ay_d_u = Sigma1_inv*ay_d_u;
    arma::mat Yrowj = Y.row(j);
    arma::mat C_1u = Y.row(j)*Yrowj.t()*estA*Sigma1_inv - Y.row(j)*omega1.t()*Sigma1_inv;
    arma::mat C_1d = ay_d_u.t()*S_ay_d_u;
    d_gradient_part1 = d_gradient_part1 - S_ay_d_u/(1+C_1d(0, 0));
  }

  for(int j = 0; j < n-1; ++j){
    arma::mat ax_d_u = A*X.row(j) - omega2;
    arma::mat S_ax_d_u = Sigma2_inv*ax_d_u;
    arma::mat Xrowj = X.row(j);
    arma::mat C_2u = X.row(j)*Xrowj.t()*estA*Sigma2_inv - X.row(j)*omega2.t()*Sigma2_inv;
    arma::mat C_2d = ax_d_u.t()*S_ax_d_u;
    d_gradient_part2 = d_gradient_part2 - S_ax_d_u/(1+C_2d(0, 0));
  }
  arma::mat d_gradient = d_gradient_part1/m + d_gradient_part2/n;

  return d_gradient;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat select_stepsize_for_A_t_dis(arma::mat X, arma::mat Y, arma::rowvec mu_x, arma::rowvec mu_y, arma::mat estA, arma::rowvec estd, double ll,
                                arma::mat A_gradient, uint32_t p, arma::mat Sigma1, arma::mat Sigma2, double nu, double gamma, double down){

  arma::mat A_prime = estA;
  arma::mat selected = estA;

  for(int i = 0; i < p; ++i){
    for(int j = 0; j < p; ++j){

      double gra_now = A_gradient(i, j);
      double gra_now2 = gra_now*gra_now*gamma;
      double start = sqrt(abs(estA(i, j)/gra_now))/2;
      double aa = start;

      while(aa > 0){
        double aa2 = aa*aa;
        A_prime(i, j) = estA(i, j) - aa2*gra_now;
        double lA_prime = likelihood_para_t_cpp(X, Y, mu_x, mu_y, A_prime, estd, Sigma1, Sigma2, nu, true);
        if(lA_prime - ll + aa2*gra_now2 > 0){
          selected(i, j)= A_prime(i, j);
          break;
        }
        aa = aa - start*down;
      }

    }
  }
  return selected;
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::rowvec select_stepsize_for_d_t_dis(arma::mat X, arma::mat Y, arma::rowvec mu_x, arma::rowvec mu_y, arma::mat estA, arma::rowvec estd, double ll,
                                   arma::rowvec d_gradient, uint32_t p, arma::mat Sigma1, arma::mat Sigma2, double nu, double gamma, double down){

  arma::rowvec d_prime = estd;
  arma::rowvec selected = estd;

  for(int j = 0; j < p; ++j){

    double gra_now = d_gradient(j);
    double gra_now2 = gra_now*gra_now*gamma;
    double start = sqrt(abs(estd(j)/gra_now))/2;

    double aa = start;

    while(aa > 0){
      double aa2 = aa*aa;
      d_prime(j) = estd(j) - aa2*gra_now;
      double ld_prime = likelihood_para_t_cpp(X, Y, mu_x, mu_y, estA, d_prime, Sigma1, Sigma2, nu, true);
      if(ld_prime - ll + aa2*gra_now2 > 0){
        selected(j) = d_prime(j);
        break;
      }
      aa = aa - start*down;
    }
  }

  return selected;
}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List run_est_with_gradient_descent_t_dis_cpp(arma::mat X, arma::mat Y, arma::mat x1, arma::mat y1,
                                                   double gamma, double down, int steps = 1000, double nu = 5){


  double nu_factor = nu/(nu-2);
  uint32_t n = X.n_rows;
  uint32_t m = Y.n_rows;
  uint32_t p = X.n_cols;

  arma::rowvec mu_x = mean(x1, 0);
  arma::rowvec mu_y = mean(y1, 0);

  arma::mat Sigma1 = callmtCov(y1, nu)/nu_factor;
  arma::mat Sigma2 = callmtCov(y1, nu)/nu_factor;

  arma::rowvec estd = (mean(X, 0) - mu_x + mean(Y, 0) - mu_y)/2;

  arma::mat X_t = X.t();
  arma::mat Y_t = Y.t();
  arma::mat XX = X_t*X;
  arma::mat YY = Y.t()*Y;

  arma::mat Sigma1_inv = inv(Sigma1);
  arma::mat Sigma2_inv = inv(Sigma2);
  arma::mat Sigma12_inv_sum = Sigma1_inv + Sigma2_inv;

  arma::mat omega1 = arma::zeros<arma::mat>(m, p);

  for(int j = 0; j < m-1; ++j){
    omega1.row(j) = mu_y + estd;
  }

  arma::mat omega2 = arma::zeros<arma::mat>(n, p);
  for(int j = 0; j < n-1; ++j){
    omega2.row(j) = mu_x + estd;
  }

  arma::mat beta1 = inv(YY)*Y_t*omega1;
  arma::mat beta2 = inv(XX)*X_t*omega2;

  arma::mat estA = (beta1 + beta2)/2;

  double likelihood_est = likelihood_para_t_cpp(X, Y, mu_x, mu_y, estA, estd, Sigma1, Sigma2, nu, true);

  arma::mat A_gradient = gradient_A_t_dis(X, Y, Sigma1_inv, Sigma2_inv, mu_x, mu_y, estA, estd, n,  m, p);

  arma::rowvec d_gradient = gradient_d_t_dis(X, Y, Sigma1_inv, Sigma2_inv, mu_x, mu_y, estA, estd, n, m, p);

  arma::mat A_prime = arma::zeros<arma::mat>(p, p);
  arma::rowvec d_prime = arma::zeros<arma::rowvec>(p);

  for(int kk = 0; kk < steps-1; ++kk){

    arma::mat absd = abs(d_gradient);
    double val = absd.min();

    if( val >= 0.0001){
      d_prime = select_stepsize_for_d_t_dis(X, Y,  mu_x, mu_y, estA, estd, likelihood_est,
                                            d_gradient, p, Sigma1, Sigma2, nu, gamma, down);
    } else {
      d_prime = estd;
    }

    arma::mat absA = abs(A_gradient);
    val = absA.min();

    if( val >= 0.0001){
      A_prime = select_stepsize_for_A_t_dis(X, Y, mu_x, mu_y, estA, estd, likelihood_est,
                                            A_gradient, p, Sigma1, Sigma2, nu, gamma, down);
    } else {
      A_prime = estA;
    }

    estA = A_prime; estd = d_prime;

    A_gradient = gradient_A_t_dis(X, Y, Sigma1_inv, Sigma2_inv, mu_x, mu_y, estA, estd, n,  m, p);
    d_gradient = gradient_d_t_dis(X, Y, Sigma1_inv, Sigma2_inv, mu_x, mu_y, estA, estd, n, m, p);

    likelihood_est = likelihood_para_t_cpp(X, Y, mu_x, mu_y, estA, estd, Sigma1, Sigma2, nu, true);

  }

  return Rcpp::List::create(Rcpp::Named("A") = estA,
                            Rcpp::Named("d") = estd);

}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double gradient_theta_t_dis(arma::mat X, arma::mat Y, arma::mat Sigma1_inv, arma::mat Sigma2_inv, arma::rowvec mu_x, arma::rowvec mu_y, double estTheta,
                               arma::mat estA, arma::rowvec estd, uint32_t n, uint32_t m, uint32_t p){

  arma::mat A_dev = rotation_matrix_dev(estTheta);
  arma::mat A_dev_T = A_dev.t();

  arma::mat A_T = estA.t();

  arma::mat omega1 = A_T*Sigma1_inv*A_dev;
  arma::mat omega2 = A_T*Sigma2_inv*A_dev;

  arma::mat omega3 = A_dev_T*(mu_y + estd);
  arma::mat omega4 = A_dev_T*(mu_x + estd);

  arma::mat omega5 = arma::zeros<arma::mat>(m, p);
  arma::mat omega6 = arma::zeros<arma::mat>(n, p);


  for(int j = 0; j < m-1; ++j){
    omega5.row(j) = mu_y + estd;
  }

  for(int j = 0; j < n-1; ++j){
    omega6.row(j) = mu_x + estd;
  }

  arma::mat theta_gradient_part2 = arma::zeros<arma::mat>(1, 1);
  arma::mat theta_gradient_part1 = arma::zeros<arma::mat>(1, 1);

  for(int j = 0; j < m-1; ++j){
    arma::mat ay_d_u = A_T*Y.row(j) - omega5;
    arma::mat S_ay_d_u = Sigma1_inv*ay_d_u;
    arma::mat C_1u =  Y.row(j)*omega1*Y.col(j) - Y.row(j)*omega3;
    arma::mat C_1d = ay_d_u.t()*S_ay_d_u;
    theta_gradient_part1 = theta_gradient_part1 + C_1u/(1+C_1d(0, 0));
  }

  for(int j = 0; j < n-1; ++j){
    arma::mat ax_d_u = A_T*X.row(j) - omega6;
    arma::mat S_ax_d_u = Sigma2_inv*ax_d_u;
    arma::mat C_2u = X.row(j)*omega2*X.col(j) - X.row(j)*omega4;
    arma::mat C_2d = ax_d_u.t()*S_ax_d_u;
    theta_gradient_part2 = theta_gradient_part2 + C_2u/(1+C_2d(0, 0));
  }
  double res = theta_gradient_part1(0, 0)/m + theta_gradient_part2(0, 0)/n;

  return res;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double select_stepsize_for_theta_t_dis(arma::mat X, arma::mat Y, arma::rowvec mu_x, arma::rowvec mu_y, double estTheta, arma::mat estA, arma::rowvec estd, double ll,
                                 double Theta_gradient, uint32_t p, arma::mat Sigma1, arma::mat Sigma2, double nu, double gamma, double down){

  double Theta_prime = estTheta;
  double selected = estTheta;

  double gra_now = Theta_gradient;
  double gra_now2 = gra_now*gra_now*gamma;
  double start = sqrt(abs(estTheta/gra_now))/2;
  double aa = start;

  arma::mat A_prime = estA;
  while(aa > 0){
    double aa2 = aa*aa;
    Theta_prime = estTheta + aa2*gra_now;
    A_prime = rotation_matrix(Theta_prime);
    double lA_prime = likelihood_para_t_cpp(X, Y, mu_x, mu_y, A_prime, estd, Sigma1, Sigma2, nu, true);
    if(lA_prime - ll - aa2*gra_now2 > 0){
      selected = Theta_prime;
      break;
    }
    aa = aa - start*down;
  }

  return selected;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List run_est_with_gradient_descent_angle_search_t_dis_cpp(arma::mat X, arma::mat Y, arma::mat x1, arma::mat y1,
                                                   double gamma, double down, int steps = 1000, double nu = 5){


  double nu_factor = nu/(nu-2);
  uint32_t n = X.n_rows;
  uint32_t m = Y.n_rows;
  uint32_t p = X.n_cols;

  arma::rowvec mu_x = mean(x1, 0);
  arma::rowvec mu_y = mean(y1, 0);

  arma::mat Sigma1 = callmtCov(y1, nu)/nu_factor;
  arma::mat Sigma2 = callmtCov(y1, nu)/nu_factor;

  arma::rowvec estd = (mean(X, 0) - mu_x + mean(Y, 0) - mu_y)/2;

  arma::mat Sigma1_inv = inv(Sigma1);
  arma::mat Sigma2_inv = inv(Sigma2);

  double estTheta = 0;
  arma::mat estA = rotation_matrix(estTheta);

  double likelihood_est = likelihood_para_t_cpp(X, Y, mu_x, mu_y, estA, estd, Sigma1, Sigma2, nu, true);

  double Theta_gradient = gradient_theta_t_dis(X, Y, Sigma1_inv, Sigma2_inv, mu_x, mu_y, estTheta, estA, estd, n, m, p);

  arma::rowvec d_gradient = gradient_d_t_dis(X, Y, Sigma1_inv, Sigma2_inv, mu_x, mu_y, estA, estd, n, m, p);

  arma::mat A_prime = arma::zeros<arma::mat>(p, p);;
  arma::rowvec d_prime = arma::zeros<arma::rowvec>(p);
  double Theta_prime = 0;

  for(int kk = 0; kk < steps-1; ++kk){

    arma::mat absd = abs(d_gradient);
    double val = absd.min();

    if( val >= 0.0001){
      d_prime = select_stepsize_for_d_t_dis(X, Y,  mu_x, mu_y, estA, estd, likelihood_est,
                                            d_gradient, p, Sigma1, Sigma2, nu, gamma, down);
    } else {
      d_prime = estd;
    }

    if( abs(Theta_gradient) >= 0.0001){
      Theta_prime = select_stepsize_for_theta_t_dis(X, Y, mu_x, mu_y, estTheta, estA, estd, likelihood_est,
                                                    Theta_gradient, p, Sigma1, Sigma2, nu, gamma, down);
      A_prime = rotation_matrix(estTheta);
    } else {
      Theta_prime = estTheta;
      A_prime = estA;
    }

    estA = A_prime; estd = d_prime; estTheta = Theta_prime;

    Theta_gradient = gradient_theta_t_dis(X, Y, Sigma1_inv, Sigma2_inv, mu_x, mu_y, estTheta, estA, estd, n, m, p);
    d_gradient = gradient_d_t_dis(X, Y, Sigma1_inv, Sigma2_inv, mu_x, mu_y, estA, estd, n, m, p);

    likelihood_est = likelihood_para_t_cpp(X, Y, mu_x, mu_y, estA, estd, Sigma1, Sigma2, nu, true);

  }

  return Rcpp::List::create(Rcpp::Named("Theta") = estTheta,
                            Rcpp::Named("A") = estA,
                            Rcpp::Named("d") = estd);

}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List run_simple_angle_est_cpp(arma::mat X, arma::mat Y, arma::mat x1, arma::mat y1, arma::vec angle_list, int steps = 1000){

  uint32_t n = X.n_rows;
  uint32_t m = Y.n_rows;
  uint32_t p = X.n_cols;
  uint32_t al = angle_list.n_elem;

  arma::rowvec mu_x = mean(x1, 0);
  arma::rowvec mu_y = mean(y1, 0);
  arma::mat D1 = cov(y1);
  arma::mat D2 = cov(x1);

  arma::rowvec estd = (mean(X, 0) + mean(Y, 0)- mu_x - mu_y)/2;

  arma::mat X_t = X.t();
  arma::mat Y_t = Y.t();
  arma::mat XX = X_t*X;
  arma::mat YY = Y.t()*Y;

  arma::mat D1_inv = inv(D1);
  arma::mat D2_inv = inv(D2);
  arma::mat D12_inv_sum = D1_inv + D2_inv;

  arma::mat omega1 = arma::zeros<arma::mat>(m, p);

  for(int j = 0; j < m-1; ++j){
    omega1.row(j) = mu_y + estd;
  }

  arma::mat omega2 = arma::zeros<arma::mat>(n, p);
  for(int j = 0; j < n-1; ++j){
    omega2.row(j) = mu_x + estd;
  }

  arma::mat iniA = rotation_matrix(0);
  arma::mat estA = iniA;
  arma::mat C1 = YY*estA*D1_inv/m + XX*estA*D2_inv/n;
  arma::mat C2 = Y_t*omega1*D1_inv/m + X_t*omega2*D2_inv/n;
  arma::mat A_gradient = C1 - C2;
  arma::rowvec d_gradient = estd;
  arma::vec changes = arma::zeros<arma::vec>(al);
  double selected = 0;

  for(int kk = 0; kk < steps-1; ++kk){

    for(int j = 0; j < m-1; ++j){
      omega1.row(j) = mu_y + estd;
    }

    for(int j = 0; j < n-1; ++j){
      omega2.row(j) = mu_x + estd;
    }

    C1 = YY*estA*D1_inv/m + XX*estA*D2_inv/n;
    C2 = Y_t*omega1*D1_inv/m + X_t*omega2*D2_inv/n;
    A_gradient = C1 - C2;

    for(int j = 0; j < al-1; ++j){
      arma::mat A_t_prime = rotation_matrix(angle_list(j));
      C1 = YY*A_t_prime*D1_inv/m + XX*A_t_prime*D2_inv/n;
      C2 = Y_t*omega1*D1_inv/m + X_t*omega2*D2_inv/n;
      arma::mat med = median(abs((C1 - C2)/A_gradient));
      changes(j) = med(0, 0);
    }

    arma::mat XestA = X*estA;
    arma::mat YestA = Y*estA;
    arma::rowvec errorX = arma::zeros<arma::rowvec>(p);
    arma::rowvec errorY = arma::zeros<arma::rowvec>(p);
    for(int j = 0; j < n-1; ++j){
      errorX = errorX + XestA.row(j) - mu_x;
    }
    for(int j = 0; j < m-1; ++j){
      errorY = errorY + YestA.row(j) - mu_y;
    }
    arma::rowvec ave_errorX = errorX/n;
    arma::rowvec ave_errorY = errorY/m;
    arma::rowvec dmat_part = ave_errorX*D2_inv + ave_errorY*D1_inv;
    d_gradient = estd*D12_inv_sum - dmat_part;
    arma::uvec aa_ind = arma::find(changes < 1);
    if(aa_ind.n_elem >= 1){
      selected = angle_list(aa_ind(0));
      estA = rotation_matrix(selected);
    } else {
      estA = iniA;
    }
    estd = estd - 0.01*d_gradient;

  }

  Rcpp::Rcout << "d" << d_gradient  << std::endl;
  Rcpp::Rcout << "A" << A_gradient  << std::endl;
  return Rcpp::List::create(Rcpp::Named("A") = estA,
                            Rcpp::Named("d") = estd);

}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List run_simple_angle_likelihood_est_cpp(arma::mat X, arma::mat Y, arma::mat x1, arma::mat y1, arma::vec angle_list, int steps = 1000){

  uint32_t n = X.n_rows;
  uint32_t m = Y.n_rows;
  uint32_t p = X.n_cols;
  uint32_t al = angle_list.n_elem;

  arma::rowvec mu_x = mean(x1, 0);
  arma::rowvec mu_y = mean(y1, 0);
  arma::mat D1 = cov(y1);
  arma::mat D2 = cov(x1);

  arma::rowvec estd = (mean(X, 0) + mean(Y, 0)- mu_x - mu_y)/2;

  arma::mat X_t = X.t();
  arma::mat Y_t = Y.t();
  arma::mat XX = X_t*X;
  arma::mat YY = Y.t()*Y;

  arma::mat D1_inv = inv(D1);
  arma::mat D2_inv = inv(D2);
  arma::mat D12_inv_sum = D1_inv + D2_inv;

  arma::mat omega1 = arma::zeros<arma::mat>(m, p);

  for(int j = 0; j < m-1; ++j){
    omega1.row(j) = mu_y + estd;
  }

  arma::mat omega2 = arma::zeros<arma::mat>(n, p);
  for(int j = 0; j < n-1; ++j){
    omega2.row(j) = mu_x + estd;
  }

  arma::mat iniA = rotation_matrix(0);
  arma::mat estA = iniA;
  arma::mat C1 = YY*estA*D1_inv/m + XX*estA*D2_inv/n;
  arma::mat C2 = Y_t*omega1*D1_inv/m + X_t*omega2*D2_inv/n;
  arma::mat A_gradient = C1 - C2;
  arma::rowvec d_gradient = estd;
  arma::vec changes = arma::zeros<arma::vec>(al);
  double selected = 0;

  for(int kk = 0; kk < steps-1; ++kk){

    for(int j = 0; j < m-1; ++j){
      omega1.row(j) = mu_y + estd;
    }

    for(int j = 0; j < n-1; ++j){
      omega2.row(j) = mu_x + estd;
    }

    C1 = YY*estA*D1_inv/m + XX*estA*D2_inv/n;
    C2 = Y_t*omega1*D1_inv/m + X_t*omega2*D2_inv/n;
    A_gradient = C1 - C2;

    for(int j = 0; j < al-1; ++j){
      arma::mat A_t_prime = rotation_matrix(angle_list(j));
      changes(j) =  likelihood_para_cpp(X, Y, mu_x, mu_y, A_t_prime, estd, D1, D2, true);
    }

    arma::mat XestA = X*estA;
    arma::mat YestA = Y*estA;
    arma::rowvec errorX = arma::zeros<arma::rowvec>(p);
    arma::rowvec errorY = arma::zeros<arma::rowvec>(p);
    for(int j = 0; j < n-1; ++j){
      errorX = errorX + XestA.row(j) - mu_x;
    }
    for(int j = 0; j < m-1; ++j){
      errorY = errorY + YestA.row(j) - mu_y;
    }
    arma::rowvec ave_errorX = errorX/n;
    arma::rowvec ave_errorY = errorY/m;
    arma::rowvec dmat_part = ave_errorX*D2_inv + ave_errorY*D1_inv;
    d_gradient = estd*D12_inv_sum - dmat_part;
    arma::uword aa_ind = index_max(changes);
    selected = angle_list(aa_ind);
    estA = rotation_matrix(selected);

    estd = estd - 0.01*d_gradient;

  }

  return Rcpp::List::create(Rcpp::Named("A") = estA,
                            Rcpp::Named("d") = estd);

}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List run_simple_est_theta_coordinate_descent_cpp(arma::mat X, arma::mat Y, arma::mat x1, arma::mat y1, int steps = 20, int gra_steps = 10){

  uint32_t n = X.n_rows;
  uint32_t m = Y.n_rows;
  uint32_t p = X.n_cols;

  arma::rowvec mu_x = mean(x1, 0);
  arma::rowvec mu_y = mean(y1, 0);
  arma::mat D1 = cov(y1);
  arma::mat D2 = cov(x1);

  arma::rowvec estd = -(mean(X, 0) - mu_x + mean(Y, 0) - mu_y)/2;

  arma::mat X_t = X.t();
  arma::mat Y_t = Y.t();

  arma::mat D1_inv = inv(D1);
  arma::mat D2_inv = inv(D2);
  arma::mat D12_inv_sum = D1_inv + D2_inv;
  arma::mat D12_sum_inv = inv(D12_inv_sum);



  arma::mat Theta_gradient = arma::zeros<arma::mat>(1, 1);
  double start = 0.01;
  double estTheta = 0;
  arma::mat estA = rotation_matrix(estTheta);

  for(int kk = 0; kk < steps-1; ++kk){

    arma::mat XestA = X*estA;
    arma::mat YestA = Y*estA;
    arma::rowvec errorX = arma::zeros<arma::rowvec>(p);
    arma::rowvec errorY = arma::zeros<arma::rowvec>(p);
    for(int j = 0; j < n-1; ++j){
      errorX = errorX + XestA.row(j);
    }
    for(int j = 0; j < m-1; ++j){
      errorY = errorY + YestA.row(j);
    }
    arma::rowvec ave_errorX = errorX/n - mu_x;
    arma::rowvec ave_errorY = errorY/m - mu_y;
    arma::rowvec dmat_part = ave_errorX*D2_inv + ave_errorY*D1_inv;
    estd = dmat_part*D12_sum_inv;

    for(int jj = 0; jj < gra_steps-1; ++jj){

      arma::mat A_T = estA.t();
      arma::mat A_dev = rotation_matrix_dev(estTheta);
      arma::mat A_dev_T = A_dev.t();

      arma::mat omega1 = A_T*D1_inv*A_dev;
      arma::mat omega2 = A_T*D2_inv*A_dev;

      arma::rowvec mu_y_est_d = mu_y + estd;
      arma::rowvec mu_x_est_d = mu_x + estd;
      arma::mat omega3 = A_dev_T*mu_y_est_d.t();
      arma::mat omega4 = A_dev_T*mu_x_est_d.t();

      arma::mat Theta_gradient1 = arma::zeros<arma::mat>(1, 1);
      for(int j = 0; j < m-1; ++j){
        Theta_gradient1 +=  Y.row(j)*omega1*Y_t.col(j) - Y.row(j)*omega3;
      }

      arma::mat Theta_gradient2 = arma::zeros<arma::mat>(1, 1);

      for(int j = 0; j < n-1; ++j){
        Theta_gradient2 +=  X.row(j)*omega2*X_t.col(j) - X.row(j)*omega4;
      }
      Theta_gradient = Theta_gradient1/m + Theta_gradient2/n;

  //    Rcpp::Rcout << "A_gradient " << Theta_gradient(0, 0) << std::endl;

      if(jj == 0){
      if(abs(Theta_gradient(0, 0)) > 0.1){
        arma::mat pre = arma::ceil(log10(abs(Theta_gradient)));
        start = pre.max();
        start = pow(10, -start-1);
      } else {
        start = 0.001;
      }
      }

     // Rcpp::Rcout << "start" << start << std::endl;
      estTheta = estTheta - start*Theta_gradient(0, 0)/(jj+1);
      estA = rotation_matrix(estTheta);

    }
  if(abs(Theta_gradient(0, 0)) < 0.01){
    break;
  }
  }



  return Rcpp::List::create(Rcpp::Named("Theta") = estTheta,
                            Rcpp::Named("A") = estA,
                            Rcpp::Named("d") = estd);

}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List run_simple_est_theta_coordinate_descent_one_cluster_cpp(arma::mat X, arma::mat x1,
                                                                   int steps = 20, int gra_steps = 20){

  uint32_t n = X.n_rows;
  uint32_t p = X.n_cols;

  arma::rowvec mu_x = mean(x1, 0);
  arma::mat D2 = cov(x1);

  arma::rowvec estd = -(mean(X, 0) - mu_x)/2;

  arma::mat X_t = X.t();

  arma::mat D2_inv = inv(D2);

  arma::rowvec pre2 = estd;
  arma::rowvec d_gradient = estd;
  arma::mat Theta_gradient = arma::zeros<arma::mat>(1, 1);
  double start = 0;
  double estTheta = 0;
  arma::mat estA = rotation_matrix(estTheta);

  for(int kk = 0; kk < steps-1; ++kk){

    arma::mat XestA = X*estA;
    arma::rowvec errorX = arma::zeros<arma::rowvec>(p);
    for(int j = 0; j < n-1; ++j){
      errorX = errorX + XestA.row(j);
    }

    arma::rowvec ave_errorX = errorX/n - mu_x;
    estd = ave_errorX;

    for(int jj = 0; jj < gra_steps-1; ++jj){
      arma::mat A_T = estA.t();
      arma::mat A_dev = rotation_matrix_dev(estTheta);
      arma::mat A_dev_T = A_dev.t();

      arma::mat omega2 = A_T*D2_inv*A_dev;

      arma::rowvec mu_x_est_d = mu_x + estd;
      arma::mat omega4 = A_dev_T*mu_x_est_d.t();

      arma::mat Theta_gradient2 = arma::zeros<arma::mat>(1, 1);
      for(int j = 0; j < n-1; ++j){
        Theta_gradient2 +=  X.row(j)*omega2*X_t.col(j) - X.row(j)*omega4;
      }
      Theta_gradient = Theta_gradient2/n;

      if(jj == 0){
        if(abs(Theta_gradient(0, 0)) > 0.1){
          arma::mat pre = arma::ceil(log10(abs(Theta_gradient)));
          start = pre.max();
          start = pow(10, -start-1);
        } else {
          start = 0.001;
        }
      }

  //    Rcpp::Rcout << "start" << start << std::endl;
      estTheta = estTheta - start*Theta_gradient(0, 0)/(jj+1);
      estA = rotation_matrix(estTheta);

    }

    if(abs(Theta_gradient(0, 0)) < 0.01){
      break;
    }

  }



  return Rcpp::List::create(Rcpp::Named("Theta") = estTheta,
                            Rcpp::Named("A") = estA,
                            Rcpp::Named("d") = estd);

}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List run_simple_est_theta_coordinate_descent_three_cluster_cpp(arma::mat X, arma::mat Y,
                            arma::mat Z, arma::mat x1, arma::mat y1, arma::mat z1, int steps = 20, int gra_steps = 10){

  uint32_t n = X.n_rows;
  uint32_t m = Y.n_rows;
  uint32_t l = Z.n_rows;
  uint32_t p = X.n_cols;

  arma::rowvec mu_x = mean(x1, 0);
  arma::rowvec mu_y = mean(y1, 0);
  arma::rowvec mu_z = mean(z1, 0);
  arma::mat D1 = cov(y1);
  arma::mat D2 = cov(x1);
  arma::mat D3 = cov(z1);

  arma::rowvec estd = -(mean(X, 0) - mu_x + mean(Y, 0) - mu_y + mean(Z, 0) - mu_z)/2;

  arma::mat X_t = X.t();
  arma::mat Y_t = Y.t();
  arma::mat Z_t = Z.t();

  arma::mat D1_inv = inv(D1);
  arma::mat D2_inv = inv(D2);
  arma::mat D3_inv = inv(D3);
  arma::mat D123_inv_sum = D1_inv + D2_inv + D3_inv;
  arma::mat D123_sum_inv = inv(D123_inv_sum);


  arma::mat Theta_gradient = arma::zeros<arma::mat>(1, 1);
  double start = 0.01;
  double estTheta = 0;
  arma::mat estA = rotation_matrix(estTheta);

  for(int kk = 0; kk < steps-1; ++kk){

    arma::mat XestA = X*estA;
    arma::mat YestA = Y*estA;
    arma::mat ZestA = Z*estA;
    arma::rowvec errorX = arma::zeros<arma::rowvec>(p);
    arma::rowvec errorY = arma::zeros<arma::rowvec>(p);
    arma::rowvec errorZ = arma::zeros<arma::rowvec>(p);
    for(int j = 0; j < n-1; ++j){
      errorX = errorX + XestA.row(j);
    }
    for(int j = 0; j < m-1; ++j){
      errorY = errorY + YestA.row(j);
    }
    for(int j = 0; j < l-1; ++j){
      errorZ = errorZ + ZestA.row(j);
    }
    arma::rowvec ave_errorX = errorX/n - mu_x;
    arma::rowvec ave_errorY = errorY/m - mu_y;
    arma::rowvec ave_errorZ = errorZ/l - mu_z;
    arma::rowvec dmat_part = ave_errorX*D2_inv + ave_errorY*D1_inv + ave_errorZ*D3_inv;
    estd = dmat_part*D123_sum_inv;

    for(int jj = 0; jj < gra_steps-1; ++jj){

      arma::mat A_T = estA.t();
      arma::mat A_dev = rotation_matrix_dev(estTheta);
      arma::mat A_dev_T = A_dev.t();

      arma::mat omega1 = A_T*D1_inv*A_dev;
      arma::mat omega2 = A_T*D2_inv*A_dev;
      arma::mat omega3 = A_T*D3_inv*A_dev;

      arma::rowvec mu_y_est_d = mu_y + estd;
      arma::rowvec mu_x_est_d = mu_x + estd;
      arma::rowvec mu_z_est_d = mu_z + estd;
      arma::mat omega4 = A_dev_T*mu_y_est_d.t();
      arma::mat omega5 = A_dev_T*mu_x_est_d.t();
      arma::mat omega6 = A_dev_T*mu_z_est_d.t();

      arma::mat Theta_gradient1 = arma::zeros<arma::mat>(1, 1);
      for(int j = 0; j < m-1; ++j){
        Theta_gradient1 +=  Y.row(j)*omega1*Y_t.col(j) - Y.row(j)*omega4;
      }

      arma::mat Theta_gradient2 = arma::zeros<arma::mat>(1, 1);
      for(int j = 0; j < n-1; ++j){
        Theta_gradient2 +=  X.row(j)*omega2*X_t.col(j) - X.row(j)*omega6;
      }
      arma::mat Theta_gradient3 = arma::zeros<arma::mat>(1, 1);
      for(int j = 0; j < l-1; ++j){
        Theta_gradient3 +=  Z.row(j)*omega3*Z_t.col(j) - Z.row(j)*omega6;
      }
      Theta_gradient = Theta_gradient1/m + Theta_gradient2/n + Theta_gradient3/l;


  // Rcpp::Rcout << "A_gradient " << Theta_gradient(0, 0) << std::endl;

   if(jj == 0){
     if(abs(Theta_gradient(0, 0)) > 1){
       arma::mat pre = arma::ceil(log10(abs(Theta_gradient)));
       start = pre.max();
       start = pow(10, -start-1);
     } else {
       start = 0.001;
     }
   }

   estTheta = estTheta - start*Theta_gradient(0, 0)/(jj+1);
   estA = rotation_matrix(estTheta);

    }
    if(abs(Theta_gradient(0, 0)) < 0.01 & kk > 5){
      break;
    }

  }

  return Rcpp::List::create(Rcpp::Named("Theta") = estTheta,
                            Rcpp::Named("A") = estA,
                            Rcpp::Named("d") = estd);

}
