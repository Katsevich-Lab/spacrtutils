#include <Rcpp.h>
using namespace Rcpp;

/*
 Assumptions in this file:
 1. pEmit is shaped [p, M, K].
 2. Q is shaped [p-1, K, K].
 3. Observations x[i] range from 0..(M-1) (0-based).
 4. Forward recursion uses pEmit(i+1, x[i+1], kNext).
 */

//-----------------------//
// 1) Helpers for indexing
//-----------------------//

// pEmit has dimension [p, M, K].
inline double get_pEmit(int i, int obs, int k,
                        const NumericVector &pEmit,
                        int p, int M, int K) {
  // 0-based: i in [0..p-1], obs in [0..M-1], k in [0..K-1].
  return pEmit[i + p * (obs + M * k)];
}

// Q has dimension [p-1, K, K].
inline double get_Q(int i, int k1, int k2,
                    const NumericVector &Q,
                    int pMinusOne,  // = p-1
                    int K_) {
  // 0-based: i in [0..(p-2)], k1 in [0..K-1], k2 in [0..K-1].
  // The total length is (p-1)*K*K.
  return Q[i + pMinusOne * (k1 + K_ * k2)];
}


//----------------------------//
// 2) Forward Probability
//----------------------------//
//
// Standard HMM forward recursion:
//   A(0, k) = pInit[k] * pEmit(0, x[0], k)
//   For i in [0..p-2]:
//     A(i+1, kNext) = sum_{kPrev} [ A(i, kPrev) * Q(i, kPrev, kNext ) * pEmit(i+1, x[i+1], kNext) ]
//
// [[Rcpp::export]]
NumericMatrix compute_forward_prob_cpp(
    const IntegerVector &x,     // Observations (size p)
    const NumericVector &pInit, // Initial prob vector (size K)
    const NumericVector &pEmit, // p x M x K array
    const NumericVector &Q,     // (p-1) x K x K array
    int p,                      // number of positions (SNPs)
    int M,                      // number of possible obs values
    int K                       // number of hidden states
) {
  // The dimension of Q is (p-1)*K*K
  int pTrans = p - 1;  // number of transitions

  // Forward matrix A: size p x K
  NumericMatrix A(p, K);

  // 1) Initialize the first row: multiply by emission at position 0
  for(int kState = 0; kState < K; kState++) {
    double init_val = pInit[kState];
    A(0, kState) = init_val;
  }

  // 2) Recurrence: iPos in [0..(p-2)]
  for(int iPos = 0; iPos < p - 1; iPos++) {
    for(int kNext = 0; kNext < K; kNext++) {
      double sumVal = 0.0;
      for(int kPrev = 0; kPrev < K; kPrev++) {
        double fVal = A(iPos, kPrev);
        double qVal = get_Q(iPos, kPrev, kNext, Q, pTrans, K);
        // Emission at the current position (iPos) for the next hidden state (kNext):
        double eVal = get_pEmit(iPos, x[iPos], kPrev, pEmit, p, M, K);
        sumVal += fVal * qVal * eVal;
      }
      A(iPos + 1, kNext) = sumVal;
    }
  }

  return A;
}


//----------------------------//
// 3) Backward Probability
//----------------------------//
//
// Standard HMM backward recursion:
//   B(p-1, k) = 1
//   For i in [p-2..0]:
//     B(i, kState) = sum_{kNext} [ Q(i, kState, kNext)* pEmit(i+1, x[i+1], kNext)* B(i+1, kNext) ]
//
// [[Rcpp::export]]
NumericMatrix compute_backward_prob_cpp(
    const IntegerVector &x,
    const NumericVector &pInit,  // same signature, though not strictly used here
    const NumericVector &pEmit,  // p x M x K
    const NumericVector &Q,      // (p-1) x K x K
    int p,
    int M,
    int K
) {
  int pTrans = p - 1;  // dimension for Q

  // Backward matrix B: size p x K
  NumericMatrix B(p, K);

  // 1) Initialize the last row
  for(int kState = 0; kState < K; kState++) {
    B(p - 1, kState) = 1.0;
  }

  // 2) Recurrence from (p-2) down to 0
  for(int iPos = p - 2; iPos >= 0; iPos--) {
    for(int kState = 0; kState < K; kState++) {
      double sumVal = 0.0;
      for(int kNext = 0; kNext < K; kNext++) {
        double bVal = B(iPos + 1, kNext);
        // Transition from kState -> kNext at position iPos:
        double qVal = get_Q(iPos, kState, kNext, Q, pTrans, K);
        // Emission at iPos+1, observation x[iPos+1], hidden state kNext
        double eVal = get_pEmit(iPos + 1, x[iPos + 1], kNext, pEmit, p, M, K);
        sumVal += bVal * qVal * eVal;
      }
      B(iPos, kState) = sumVal;
    }
  }

  return B;
}
