/*
 *
 * Copyright (c) 2014, Laurens van der Maaten (Delft University of Technology)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *    This product includes software developed by the Delft University of Technology.
 * 4. Neither the name of the Delft University of Technology nor the names of
 *    its contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY LAURENS VAN DER MAATEN ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 * EVENT SHALL LAURENS VAN DER MAATEN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 * IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 *
 */

// #include "winlibs/stdafx.h"
#ifdef _WIN32
#define _CRT_SECURE_NO_DEPRECATE
#endif

#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <sstream>
// #include "nbodyfft.h"
#include <math.h>
// #include "annoylib.h"
#include "kissrandom.h"
#include <thread>
#include <float.h>
#include <cstring>
// #include "vptree.h"
// #include "sptree.h"
#include "mytsne.h"
// #include "progress_bar/ProgressBar.hpp"
// #include "parallel_for.h"
#include "time_code.h"

#include <unistd.h>

using namespace std::chrono;

// #ifdef _WIN32
// #include "winlibs/unistd.h"
// #else
// #include <unistd.h>
// #endif

#include <functional>

#define _CRT_SECURE_NO_WARNINGS


int itTest = 0;

using namespace std;

//Helper function for printing Y at each iteration. Useful for debugging
void print_progress(int iter, double *Y, int N, int no_dims) {

    ofstream myfile;
    std::ostringstream stringStream;
    stringStream << "dat/intermediate" << iter << ".txt";
    std::string copyOfStr = stringStream.str();
    myfile.open(stringStream.str().c_str());
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < no_dims; i++) {
            myfile << Y[j * no_dims + i] << " ";
        }
        myfile << "\n";
    }
    myfile.close();
}


// Perform t-SNE
// int TSNE::run(double *X, int N, int D, double *Y, int no_dims, double perplexity, double theta, int rand_seed,
              // bool skip_random_init, int max_iter, int stop_lying_iter, int mom_switch_iter, 
              // double momentum, double final_momentum, double learning_rate, int K, double sigma,
              // int nbody_algorithm, int knn_algo, double early_exag_coeff, double *costs,
              // bool no_momentum_during_exag, int start_late_exag_iter, double late_exag_coeff, int n_trees, int search_k,
              // int nterms, double intervals_per_integer, int min_num_intervals, unsigned int nthreads, 
              // int load_affinities, int perplexity_list_length, double *perplexity_list, double df,
              // double max_step_norm) {
int runExactTSNE(double* P, int N, int D, double* Y, int no_dims, int rand_seed,
              bool skip_random_init, int max_iter, int mom_switch_iter, 
              double momentum, double final_momentum, 
              double* costs, double df, double max_step_norm, bool verbose) {
  
  
  // Determine whether we are using an exact algorithm
  // Allocate some memory
  auto *dY = (double *) malloc(N * no_dims * sizeof(double));
  auto *uY = (double *) malloc(N * no_dims * sizeof(double));
  auto *gains = (double *) malloc(N * no_dims * sizeof(double));
  if (dY == nullptr || uY == nullptr || gains == nullptr) throw std::bad_alloc();
  // Initialize gradient to zeros and gains to ones.
  for (int i = 0; i < N * no_dims; i++) uY[i] = .0;
  for (int i = 0; i < N * no_dims; i++) gains[i] = 1.0;

    // Set random seed
    if (skip_random_init != true) {
        if (rand_seed >= 0) {
            if(verbose){Rprintf("Using random seed: %d\n", rand_seed);}
            srand((unsigned int) rand_seed);
        } else {
            if(verbose){Rprintf("Using current time as random seed...\n");}
            srand(time(NULL));
        }
    }

    // Initialize solution (randomly)
    if (skip_random_init != true) {
		  if(verbose) {Rprintf("Randomly initializing the solution.\n");}
      for (int i = 0; i < N * no_dims; i++) Y[i] = randn() * .0001;
      if(verbose) {Rprintf("Y[0] = %lf\n", Y[0]);}
    } else {
		  if(verbose) {Rprintf("Using the given initialization.\n");}
    }

    if(verbose) {print_progress(0, Y, N, no_dims);}

    // Perform main training loop
    if(verbose) {Rprintf("Similarities loaded \nLearning embedding...\n");}

    std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

    if(verbose) { Rprintf("Running iterations: %d\n", max_iter); }
    for (int iter = 0; iter < max_iter; iter++) {
      itTest = iter;
      computeExactGradient(P, Y, N, no_dims, dY,df);
      // no_mementum_during_exag was = FALSE in .R
      for (int i = 0; i < N * no_dims; i++)
          gains[i] = (sign(dY[i]) != sign(uY[i])) ? (gains[i] + .2) : (gains[i] * .8);
      for (int i = 0; i < N * no_dims; i++) if (gains[i] < .01) gains[i] = .01;
      // for (int i = 0; i < N * no_dims; i++) uY[i] = momentum * uY[i] - learning_rate * gains[i] * dY[i];
      for (int i = 0; i < N * no_dims; i++) uY[i] = momentum * uY[i] - gains[i] * dY[i]; // try remove learning rate
	    // Clip the step sizes if max_step_norm is provided
      if (max_step_norm > 0) {
          for (int i=0; i<N; i++) {
              double step = 0;
              for (int j=0; j<no_dims; j++) {
                  step += uY[i*no_dims + j] * uY[i*no_dims + j];
              }
              step = sqrt(step);
              if (step > max_step_norm) {
                  for (int j=0; j<no_dims; j++) {
                      uY[i*no_dims + j] *= (max_step_norm/step);
                  }
		          }
	        }
      }
      for (int i = 0; i < N * no_dims; i++) Y[i] = Y[i] + uY[i];
        
      // Make solution zero-mean
      zeroMean(Y, N, no_dims);

      if (iter == mom_switch_iter) momentum = final_momentum;

        // Print out progress
      if ((iter+1) % 50 == 0 || iter == max_iter - 1) {
	      INITIALIZE_TIME;
        START_TIME;
        double C = .0;
        
        C = evaluateError(P, Y, N, no_dims,df, false);
        
        costs[iter] = C;       
        END_TIME("Computing Error");
        
        std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
        if(verbose) {
          Rprintf("Iteration %d (50 iterations in %.2f seconds), cost %f\n", iter+1, std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time).count()/(float)1000.0, C);
        }
        start_time = std::chrono::steady_clock::now();
      }
    }

    if(verbose) {Rprintf("All iterations done, cleaning now ...\n");}
    usleep(100000); // pause a little bit to print out information
    
    // Clean up memory
    free(dY);
    free(uY);
    free(gains);
    // free(P);
    if(verbose) {Rprintf("Cleanup done ...\n");}
    usleep(100000); // pause a little bit to print out information
    return 0;
}


void computeExactGradientTest(double* Y, int N, int D, double df ) {
  // Compute the squared Euclidean distance matrix
    double *DD = (double *) malloc(N * N * sizeof(double));
    if (DD == NULL) {
        Rprintf("Memory allocation failed!\n");
        exit(1);
    }
    computeSquaredEuclideanDistance(Y, N, D, DD);

    // Compute Q-matrix and normalization sum
    double *Q = (double *) malloc(N * N * sizeof(double));
    if (Q == NULL) {
        Rprintf("Memory allocation failed!\n");
        exit(1);
    }
    double sum_Q = .0;
    int nN = 0;
    for (int n = 0; n < N; n++) {
        for (int m = 0; m < N; m++) {
            if (n != m) {
                Q[nN + m] = 1.0 / pow(1.0 + DD[nN + m]/(double)df, (df));
                sum_Q += Q[nN + m];
            }
        }
        nN += N;
    }

    // Perform the computation of the gradient
    char buffer[500];
    sprintf(buffer, "temp/exact_gradient%d.txt", itTest);
    FILE *fp = fopen(buffer, "w"); // Open file for writing
    nN = 0;
    int nD = 0;
    for (int n = 0; n < N; n++) {
        double testQij = 0;
        double testPos = 0;
        double testNeg1 = 0;
        double testNeg2 = 0;
        double testdC = 0;
        int mD = 0;
        for (int m = 0; m < N; m++) {
            if (n != m) {
                testNeg1 += pow(Q[nN + m],(df +1.0)/df) * (Y[nD + 0] - Y[mD + 0]) / sum_Q;
                testNeg2 += pow(Q[nN + m],(df +1.0)/df) * (Y[nD + 1] - Y[mD + 1]) / sum_Q;
            }
            mD += D;
        }
        fprintf(fp, "%d, %.12e, %.12e\n", n, testNeg1,testNeg2);

        nN += N;
        nD += D;
    }
    fclose(fp);
    free(DD);
    free(Q);

}


// Compute the exact gradient of the t-SNE cost function
void computeExactGradient(double* P, double* Y, int N, int D, double* dC, double df) {
    // Make sure the current gradient contains zeros
    for (int i = 0; i < N * D; i++) dC[i] = 0.0;

    // Compute the squared Euclidean distance matrix
    auto *DD = (double *) malloc(N * N * sizeof(double));
    if (DD == nullptr) throw std::bad_alloc();
    computeSquaredEuclideanDistance(Y, N, D, DD);

    // Compute Q-matrix and normalization sum
    auto *Q = (double *) malloc(N * N * sizeof(double));
    if (Q == nullptr) throw std::bad_alloc();

    auto *Qpow = (double *) malloc(N * N * sizeof(double));
    if (Qpow == nullptr) throw std::bad_alloc();

    double sum_Q = .0;
    int nN = 0;
    for (int n = 0; n < N; n++) {
        for (int m = 0; m < N; m++) {
            if (n != m) {
                //Q[nN + m] = 1.0 / pow(1.0 + DD[nN + m]/(double)df, df);
                Q[nN + m] = 1.0 / (1.0 + DD[nN + m]/(double)df);
                Qpow[nN + m] = pow(Q[nN + m], df);
                sum_Q += Qpow[nN + m];
            }
        }
        nN += N;
    }

    // Perform the computation of the gradient
    nN = 0;
    int nD = 0;
    for (int n = 0; n < N; n++) {
        int mD = 0;
        for (int m = 0; m < N; m++) {
            if (n != m) {
                double mult = (P[nN + m] - (Qpow[nN + m] / sum_Q)) * (Q[nN + m]);
                for (int d = 0; d < D; d++) {
                    dC[nD + d] += (Y[nD + d] - Y[mD + d]) * mult;
                }
            }
            mD += D;
        }
        nN += N;
        nD += D;
    }
    free(Q);
    free(Qpow);
    free(DD);
}


// Evaluate t-SNE cost function (exactly)
double evaluateError(double* P, double* Y, int N, int D, double df, bool verbose) {
    // Compute the squared Euclidean distance matrix
    double *DD = (double *) malloc(N * N * sizeof(double));
    double *Q = (double *) malloc(N * N * sizeof(double));
    if (DD == NULL || Q == NULL) {
        Rprintf("Memory allocation failed!\n");
        exit(1);
    }
    if(verbose){ Rprintf("computeSquared\n"); }
    computeSquaredEuclideanDistance(Y, N, D, DD);

    // Compute Q-matrix and normalization sum
    if(verbose){ Rprintf("calculate Q-matrix\n"); }
    int nN = 0;
    double sum_Q = DBL_MIN;
    for (int n = 0; n < N; n++) {
        for (int m = 0; m < N; m++) {
            if (n != m) {
                //Q[nN + m] = 1.0 / pow(1.0 + DD[nN + m]/(double)df, df);
                Q[nN + m] = 1.0 / (1.0 + DD[nN + m]/(double)df);
                Q[nN +m ] = pow(Q[nN +m ], df);
                sum_Q += Q[nN + m];
            } else Q[nN + m] = DBL_MIN;
        }
        nN += N;
    }
    //Rprintf("sum_Q: %e", sum_Q);
    if(verbose){ Rprintf("normalize Q-matrix\n"); }
    for (int i = 0; i < N * N; i++) Q[i] /= sum_Q;
    //  for (int i = 0; i < N; i++) Rprintf("Q[%d]: %e\n", i, Q[i]);

//Rprintf("Q[N*N/2+1]: %e, Q[N*N-1]: %e\n", Q[N*N/2+1], Q[N*N/2+2]);

    // Sum t-SNE error
    if(verbose){ Rprintf("sum error, to %i\n", N*N); }
    double C = .0;
    for (int n = 0; n < N * N; n++) {
      C += P[n] * log((P[n] + FLT_MIN) / (Q[n] + FLT_MIN));
    }

    // Clean up memory
    free(DD);
    free(Q);
    return C;
}


// Compute squared Euclidean distance matrix
void computeSquaredEuclideanDistance(double* X, int N, int D, double* DD) {
    const double *XnD = X;
    for (int n = 0; n < N; ++n, XnD += D) {
        const double *XmD = XnD + D;
        double *curr_elem = &DD[n * N + n];
        *curr_elem = 0.0;
        double *curr_elem_sym = curr_elem + N;
        for (int m = n + 1; m < N; ++m, XmD += D, curr_elem_sym += N) {
            *(++curr_elem) = 0.0;
            for (int d = 0; d < D; ++d) {
                *curr_elem += (XnD[d] - XmD[d]) * (XnD[d] - XmD[d]);
            }
            *curr_elem_sym = *curr_elem;
        }
    }
}


// Makes data zero-mean
void zeroMean(double* X, int N, int D) {
    // Compute data mean
    double *mean = (double *) calloc(D, sizeof(double));
    if (mean == NULL) throw std::bad_alloc();

    int nD = 0;
    for (int n = 0; n < N; n++) {
        for (int d = 0; d < D; d++) {
            mean[d] += X[nD + d];
        }
        nD += D;
    }
    for (int d = 0; d < D; d++) {
        mean[d] /= (double) N;
    }

    // Subtract data mean
    nD = 0;
    for (int n = 0; n < N; n++) {
        for (int d = 0; d < D; d++) {
            X[nD + d] -= mean[d];
        }
        nD += D;
    }
    free(mean);
}


// Generates a Gaussian random number
double randn() {
    double x, y, radius;
    do {
        x = 2 * (rand() / ((double) RAND_MAX + 1)) - 1;
        y = 2 * (rand() / ((double) RAND_MAX + 1)) - 1;
        radius = (x * x) + (y * y);
    } while ((radius >= 1.0) || (radius == 0.0));
    radius = sqrt(-2 * log(radius) / radius);
    x *= radius;
    return x;
}

