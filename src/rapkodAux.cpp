#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//Auxilary function to extract a row from a Matrix while subsetting the columns 
NumericVector subr(NumericMatrix mat, const int& i, const Range& rg)
{
  NumericMatrix aux = mat( _, rg);
  return aux(i, _);
}

//Auxilary function to extract a column from a Matrix while subsetting the rows
NumericVector subc(NumericMatrix mat, const Range& rg, const int& i)
{
  NumericMatrix aux = mat(rg, _);
  return aux(_, i);
}



// [[Rcpp::export]]
List Rcpp_aux(const int& ref_n, const int& n, NumericMatrix X, NumericMatrix K, const bool& given_kern, const Function& CrossKern, //const Function& wrapf, 
              const int& p, NumericVector s2_i, const double& alpha, const bool& use_tested_inlier, const Function& rnormmat)
{
  
  NumericVector Sn_X(n - ref_n); //test statistic
  NumericVector pv(n - ref_n);  //p-values
  LogicalVector flag(n - ref_n);  //decisions to mark tested obs. as outliers or not
  
  NumericMatrix pV_X(p, n - ref_n);  //eventually will be a matrix whose columns are projections p_V(x) for each tested obs. x
  arma::mat H(p, ref_n); //auxilary matrix H (will be filled with standard gaussian entries)
  //arma::vec H_vec(H.memptr(), H.n_elem, false, false);  //vectorized version of the matrix H and that share the same memory
  
  NumericVector k_X(ref_n);
  NumericMatrix X_ref;
  X_ref = X(Range(0, ref_n -1), _);  //X_ref: corresponds to X[ref_I, ] in R writing 
         
  int oldest = 0; //row index of the oldest observation in matrix X_ref       
                                   
  for (int i = ref_n; i < n; ++i)
  {
  //i is the index of curently tested observation
  
  //Compute k_X = vector(k(x_i, y_1) ... k(x_i, y_ref.n)) where y_1,...,y_reg.n are current reference inliers
  if (!given_kern){k_X = CrossKern(X_ref, X(i, _)); }else{k_X = X_ref( _, i);}
  
    
      //H_vec = as<arma::vec>(rnorm(p*ref_n, 0, 1));
      H = as<arma::mat>(rnormmat(p, ref_n));
  
      pV_X( _, i - ref_n) = pow(ref_n, -0.5)*as<NumericVector>(wrap(H * as<arma::vec>(k_X)));
    
      Sn_X[(i - ref_n)] = sum(pow(pV_X( _, i - ref_n), 2));
       
      pv[(i - ref_n)] = mean(pchisq( Sn_X[i - ref_n] / s2_i , p) );
      flag[(i - ref_n)] = (pv[i - ref_n] < alpha);

      if (use_tested_inlier && !flag[i - ref_n]){



        if (oldest > 0) s2_i[Range(0, oldest-1)] = s2_i[Range(0, oldest-1)] - pow(subr(K, 1, Range(0, oldest-1)), 2)/(ref_n - 1) + pow(k_X[Range(0, oldest-1)], 2)/(ref_n - 1) ;
        if (oldest < ref_n-1) s2_i[Range(oldest+1, ref_n-1)] = s2_i[Range(oldest+1, ref_n-1)] - pow(subr(K, 1, Range(oldest+1, ref_n-1)), 2)/(ref_n - 1) + pow(k_X[Range(oldest+1, ref_n-1)], 2)/(ref_n - 1) ;
        s2_i[oldest] = mean(pow(k_X,2)) - pow(k_X[oldest], 2)/ref_n;



         if (oldest > 0) subc(K, Range(0, oldest-1) ,oldest) = k_X[Range(0, oldest-1)];
         if (oldest < ref_n-1) subc(K, Range(oldest+1, ref_n-1) ,oldest) = k_X[Range(oldest+1, ref_n-1)];
         if (oldest > 0) subr(K, oldest, Range(0, oldest-1)) = k_X[Range(0, oldest-1)];
         if (oldest < ref_n-1) subr(K, oldest, Range(oldest+1, ref_n-1)) = k_X[Range(oldest+1, ref_n-1)];




        X_ref(oldest, _) = X(i, _);


        oldest = (oldest + 1) % ref_n;


      }


      
  }
  
  return List::create(_["Sn.X"] = Sn_X, _["pv"] = pv, _["flag"] = flag, _["pvX"]= pV_X);
}






