#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
#include <iostream>
using namespace Rcpp;
using namespace Eigen;
using namespace std;
// [[Rcpp::depends(RcppEigen)]]

typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
typedef Eigen::Map<SparseMatrix<double> > MapMats;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Triplet<double> T;
#define max(a, b)  (((a) > (b)) ? (a) : (b))

void LassoGroup(VectorXd &x, double lambda, double rho){
  // soft thresholding rule.
  if(x.norm()<lambda/rho) x = MatrixXd::Zero(x.size(),1);
  else x = (1-(lambda/rho)/x.norm())*x;
}

void McpGroup(VectorXd &x, double lambda, double rho, double gamma){
  if(gamma > 1.0/rho){
    if(x.norm()<gamma*rho){
      LassoGroup(x, lambda, rho); x=x/(1-1/(gamma*rho));
    }
  } else {
    cout<< "error in McpGroup: gamma <1.0/rho"<<endl;
  }
}

void ScadGroup(VectorXd &x,double lambda, double rho, double gamma){
  if(gamma>1.0+1.0/rho){
    if( x.norm() < lambda + lambda/rho){
      LassoGroup(x, lambda, rho);
    }else if (x.norm() < gamma * lambda ) {
      LassoGroup(x, lambda*gamma, (gamma-1)*rho); x=x/(1-1/((gamma-1)*rho));
    }
  } else {
    cout<< "error in ScadGroup: gamma <1+1.0/rho"<<endl;
  }
}

void LassoGroup(SparseVector<double> &x, double lambda, double rho){
  // soft thresholding rule.
  if(x.norm()<lambda/rho) x.setZero();
  else x = (1-(lambda/rho)/x.norm())*x;
}

void McpGroup(SparseVector<double> &x, double lambda, double rho, double gamma){
  if(gamma > 1.0/rho){
    if(x.norm()<gamma*rho){
      LassoGroup(x, lambda, rho); x=x/(1-1/(gamma*rho));
    }
  } else {
    cout<< "error in McpGroup: gamma <1.0/rho"<<endl;
  }
}

void ScadGroup(SparseVector<double> &x,double lambda, double rho, double gamma){
  if(gamma>1.0+1.0/rho){
    if( x.norm() < lambda + lambda/rho){
      LassoGroup(x, lambda, rho);
    }else if (x.norm() < gamma * lambda ) {
      LassoGroup(x, lambda*gamma, (gamma-1)*rho); x=x/(1-1/((gamma-1)*rho));
    }
  } else {
    cout<< "error in ScadGroup: gamma <1+1.0/rho"<<endl;
  }
}

#define max(a, b)  (((a) > (b)) ? (a) : (b))

double softThresholding(double x, double lambda){
  //cout<<"LASSO"<<endl;
  if(abs(x)<lambda) return 0.0;
  else if(x>0) return x-lambda;
  else return x+lambda;
}

double MCPa(double x, double lambda, double rho, double gamma){
  double y;y=x;
  if(gamma>1.0/rho){
    cout<<"MCP"<<endl;
    double lambda1=lambda/rho; double tu1=rho*gamma;
    if(abs(x)<lambda*gamma) y=softThresholding(x,lambda1)/(1.0-1.0/tu1);
  } else {
    cout<<"error in MCPa: gamma < 1/rho"<<endl;
  }
  return y;
}

double SCADa(double x,double lambda, double rho, double gamma){
  // tuPara(0)-lambda, tuPara(1)-v, tuPara(2)-gamma.
  double y;y=x;
  if(gamma>1.0/rho+1){
    //cout<<"SCAD"<<endl;
    double lambda1=lambda/rho; double tu1=lambda*gamma;
    if(abs(x)<lambda+lambda1) y=softThresholding(x,lambda1);
    else if(abs(x)<tu1) y=softThresholding(x,tu1/((gamma-1.0)*rho))/(1.0-1.0/((gamma-1.0)*rho));
  }else {
    cout<<"error in MCPa: gamma < 1/rho"<<endl;
  }
  return y;
}

VectorXd update_eta(const VectorXd &eta1, int p_x, VectorXd &tuPara){
  VectorXd eta;
  eta = eta1;
  
  VectorXd tempeta(p_x);
  for(int k=0; k<eta.size()/p_x; k++){
    if(tuPara(3)==1.0){
      tempeta=eta.segment(k*p_x, p_x);
      LassoGroup(tempeta, tuPara(0), tuPara(1));
      eta.segment(k*p_x, p_x) = tempeta;
    }
    else if(tuPara(3)==2.0) {
      tempeta=eta.segment(k*p_x, p_x);
      ScadGroup(tempeta, tuPara(0), tuPara(1), tuPara(2));
      eta.segment(k*p_x,p_x) = tempeta;
    }
    else {
      tempeta=eta.segment(k*p_x, p_x);
      McpGroup(tempeta, tuPara(0), tuPara(1), tuPara(2));
      eta.segment(k*p_x, p_x)=tempeta;
    }
  }
  return eta;
}

int groupNoVarSelAdmm(VectorXd &Y, MatrixXd &Z, MatrixXd &X, VectorXd & struT, SpMat &A,
                      VectorXd &beta, VectorXd &eta, VectorXd &lagU, VectorXd &tuPara,
                      VectorXd &h_r_norm, VectorXd &h_s_norm, VectorXd &h_eps_pri, VectorXd &h_eps_dual,
                      VectorXd &precision, double &aa){
  // extract data from the dataset.
  int N=struT.size(); int p_z = Z.cols(); int p_x = X.cols(); // number of pcas

  // set some parameters
  int maxIter=(int)tuPara(4); int convergeStatus=-1;
  
  MatrixXd diagX(X.rows(), N*p_x); diagX.setZero();
  int  countr=0, countl=0;
  for (int i=0; i<N; i++){
    diagX.block(countr, countl, (int)struT(i), p_x) = X.block(countr, 0, (int)struT(i), p_x);
    countr = countr + struT(i);
    countl = countl + p_x;
  }
  
  MatrixXd Xt(diagX.rows(), Z.cols() + diagX.cols());
  Xt << Z, diagX;
  
  MatrixXd LS(Xt.cols(), Xt.cols()); LS.setZero();
  VectorXd ZY(Xt.cols());
  LS = Xt.transpose() * Xt;
  ZY = Xt.transpose() * Y;
  //cout << "ZY "<< ZY.transpose()<< endl;
  
  // some temp variables.
  VectorXd multA;
  VectorXd oldeta(eta); VectorXd tempeta;
  double primalR=0.0, dualS=0.0;
  
  // over-relaxition and convergence check by the program of Boyd
  double ABSTOL   = precision(0), RELTOL   = precision(1);
  double eps_pri, eps_dual;
  VectorXd Cright; //mubetaChoright.segment(0,p)=X.transpose()*Y;
  
  SpMat AA(A.cols(), A.cols());
  AA = tuPara(1)* A.transpose() * A;
  LS += AA;
  // Cholesky decomposition used in updating mu and beta.
  LDLT<MatrixXd> betaCho(LS);
  
  // The iteration of ADMM algorithm, the maximum number of iterations is maxIter,
  //  and we track the primal residual until it is less than absEps.
  
  for(int iter=0; iter<maxIter; iter++){
    // update of  beta.
    Cright = tuPara(1)* A.transpose()*(eta - lagU);
    Cright = Cright + ZY;
    //beta = cdBeta(LS, Cright, lambdabeta, beta, tuPara);
    beta = betaCho.solve(Cright);
    //cout<<"tuPara "<<tuPara.transpose()<<endl;
    //cout<<"beta "<< beta.transpose()<<endl;
    // updata of eta and theta with over-relaxition
    oldeta = eta;
    multA = A * beta;
    tempeta = multA + lagU;
    //cout<<"eta "<< eta.transpose()<<endl;
    eta = update_eta(tempeta, p_x, tuPara);
    
    //cout<<"eta "<< eta.transpose()<<endl;
    // update of lagU
    lagU=lagU+(multA-eta);
    
    // whether the algorithm has converged.
    primalR=(multA - eta).norm();
    //dualS=(tuPara(1)*A.transpose()*(eta-oldeta)).norm();
    dualS = (tuPara(1)*(eta-oldeta)).norm();
    h_r_norm(iter) = primalR;
    h_s_norm(iter) = dualS;
    
    //eps_pri = sqrt(eta.size())*ABSTOL + RELTOL*max(multA.norm(), eta.norm());
    //eps_dual= sqrt(p)*ABSTOL + RELTOL*(A.transpose()*lagU).norm();
    eps_pri = sqrt(p_x*N)*ABSTOL + RELTOL*max(beta.norm(), eta.norm());
    eps_dual= sqrt(p_x*N)*ABSTOL + RELTOL*(tuPara(1)*lagU).norm();
    
    h_eps_pri(iter)  = eps_pri;
    h_eps_dual(iter) = eps_dual;
    
    //cout<<"primal residual "<<primalR<<endl;
    //cout<<"dual residual "<<dualS<<endl;

    if (primalR < eps_pri && dualS < eps_dual)
    { convergeStatus=iter; break;}
  }
  // the RSS 
  aa = log((Y - Xt * beta).squaredNorm()/Y.size());
  
  return  convergeStatus;
}

//' @useDynLib fdagroup
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
SEXP fdagrouping(SEXP Y1, SEXP Z1, SEXP X1, SEXP struT1, SEXP A1, 
              SEXP beta1, SEXP eta1, SEXP lagU1, SEXP tuPara1, SEXP precision1){
  const MapMatd Ycopy(Rcpp::as < MapMatd>(Y1)); VectorXd Y(Ycopy);
  const MapMatd Zcopy(Rcpp::as < MapMatd>(Z1)); MatrixXd Z(Zcopy);
  const MapMatd Xcopy(Rcpp::as < MapMatd>(X1)); MatrixXd X(Xcopy);
  const MapMatd struTcopy(Rcpp::as < MapMatd>(struT1)); VectorXd struT(struTcopy);
  const MapMats Acopy(Rcpp::as < MapMats>(A1)); SpMat A(Acopy);
  
  const MapMatd betacopy(Rcpp::as < MapMatd>(beta1)); VectorXd beta(betacopy);
  const MapMatd etacopy(Rcpp::as < MapMatd>(eta1)); VectorXd eta(etacopy);
  const MapMatd lagUcopy(Rcpp::as < MapMatd>(lagU1)); VectorXd lagU(lagUcopy);
  const MapMatd tuParacopy(Rcpp::as < MapMatd>(tuPara1)); VectorXd tuPara(tuParacopy);
  const MapMatd precisioncopy(Rcpp::as < MapMatd>(precision1)); VectorXd precision(precisioncopy);
  // Z    - finite dim covariate, including mean, fixed effects, treatment and z
  // X    - functional PCA score
  // A    - augmented with 0 at the first Z.cols() columns
  // beta - all the paramters are concatenated in beta
  // eta  - augmented paramters
  // logU - Lagrangian multipliers

  // these variables are used to record the residuals 
  int maxIter=(int)tuPara(4);
  VectorXd h_r_norm(maxIter), h_s_norm(maxIter), h_eps_pri(maxIter), h_eps_dual(maxIter);
  h_r_norm.setZero(); h_s_norm.setZero(); h_eps_pri.setZero(); h_eps_dual.setZero();
  int convergeSta = -1; double aa=0.0;
  
  convergeSta=groupNoVarSelAdmm(Y, Z, X, struT, A,
                                beta, eta, lagU, tuPara,
                                h_r_norm, h_s_norm, h_eps_pri, h_eps_dual,
                                precision, aa);
  
  if(convergeSta >= 0){
    h_r_norm   = h_r_norm.head(convergeSta);
    h_s_norm   = h_s_norm.head(convergeSta);
    h_eps_pri  = h_eps_pri.head(convergeSta);
    h_eps_dual = h_eps_dual.head(convergeSta);
  }
  
  return wrap(List::create(Named("beta")=beta, Named("eta")=eta,Named("lagU")=lagU,
                                 Named("converge")=convergeSta, Named("aa")=aa,
                                 Named("h_r_norm")= h_r_norm, Named("h_s_norm")=h_s_norm,
                                 Named("h_eps_pri")=h_eps_pri, Named("h_eps_dual")=h_eps_dual ));
}





