#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar
//
//// [[Rcpp::export]]
//int timesTwo(int x) {
//   return x * 2;
//}
//
//// [[Rcpp::export]]
//int get_ai(int r,int i, int n, int *shift_ai) {
//  int ai = 0;
//	if (i>0) {
//	  ai = floor( (r  % shift_ai[i-1]) / (shift_ai[i]));
//  } else {
//	  ai = floor(r / (shift_ai[i]));
//  }
//  return(ai);
//}	
//
//
//// [[Rcpp::export]]
//void make_cheating_payoff (int *n_ptr, int *a_dim, int *shift_ai,double *g, double *c)
//{
//  int n = *n_ptr;
//	int R = 1;
//	int i = 0; int r = 0; int br = 0;
//
//
//  for(i = 0; i < n; i++) {
//	  R = R*a_dim[i];
//  }
//  for (i = 0; i<n; i++) {
//	  int i_add = R*i;
//	  for (r =0; r<R; r++) {
//		  int ai = get_ai(r,i,n,shift_ai);  // Index of the action of player i
//		  double cmax = g[i_add + r + (0-ai) * shift_ai[i]];
//		  for (br=1; br<a_dim[i]; br++) {
//		    double cact = g[i_add + r + (br-ai) * shift_ai[i]];
//		    if (cact > cmax) {
//		      cmax = cact;
//        }
//	    }
//	    c[i_add + r] = cmax;
//	    //c[i_add + r] = ai;
//    }
//  }
//  //c[1] = 5;
//}
//
//
//// [[Rcpp::export]]
//void make_cheating_payoff (int *n_ptr, int *a_dim, int *shift_ai,double *g, double *c)
//{
//	int n = *n_ptr;
//	int R = 1;
//	int i = 0; int r = 0; int br = 0;
//
//
//  for(i = 0; i < n; i++) {
//	  R = R*a_dim[i];
//  }
//  for (i = 0; i<n; i++) {
//	  int i_add = R*i;
//	  for (r =0; r<R; r++) {
//		  int ai = get_ai(r,i,n,shift_ai);  // Index of the action of player i
//		  double cmax = g[i_add + r + (0-ai) * shift_ai[i]];
//		  for (br=1; br<a_dim[i]; br++) {
//		    double cact = g[i_add + r + (br-ai) * shift_ai[i]];
//		    if (cact > cmax) {
//		      cmax = cact;
//        }
//	    }
//	    c[i_add + r] = cmax;
//	    //c[i_add + r] = ai;
//    }
//  }
//  //c[1] = 5;
//}
