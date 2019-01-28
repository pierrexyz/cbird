/* 
 * This file is part of the Copter library (http://mwhite.berkeley.edu/Copter/).
 * Copyright (C) 2007 Free Software Foundation, Inc. <http://fsf.org/>
 * 
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU General Public License as published by  
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <https://www.gnu.org/licenses/gpl-3.0.en.html>.
 */


#ifndef FFTLOG_H
#define FFTLOG_H


/* Compute the correlation function xi(r) from a power spectrum P(k), sampled
 * at logarithmically spaced points k[j]. */
void pk2xi(int N, const double k[], const double pk[], double r[], double xi[]);

/* Compute the power spectrum P(k) from a correlation function xi(r), sampled
 * at logarithmically spaced points r[i]. */
void xi2pk(int N, const double r[], const double xi[], double k[], double pk[]);

/* Compute the function
 *   \xi_l^m(r) = \int_0^\infty \frac{dk}{2\pi^2} k^m j_l(kr) P(k)
 * Note that the usual 2-point correlation function xi(r) is just xi_0^2(r)
 * in this notation.  The input k-values must be logarithmically spaced.  The
 * resulting xi_l^m(r) will be evaluated at the dual r-values
 *   r[0] = 1/k[N-1], ..., r[N-1] = 1/k[0]. */
void fftlog_ComputeXiLM(int l, int m, int N, const double k[], const double pk[], double r[], double xi[]);

#include <complex>
typedef std::complex<double> dcomplex;

/* Compute the discrete Hankel transform of the function a(r).  See the FFTLog
 * documentation (or the Fortran routine of the same name in the FFTLog
 * sources) for a description of exactly what this function computes.
 * If u is NULL, the transform coefficients will be computed anew and discarded
 * afterwards.  If you plan on performing many consecutive transforms, it is
 * more efficient to pre-compute the u coefficients. */
void fht(int N, const double r[], const dcomplex a[], double k[], dcomplex b[], double mu,
         double q = 0, double kcrc = 1, bool noring = true, dcomplex* u = NULL);

/* Pre-compute the coefficients that appear in the FFTLog implementation of
 * the discrete Hankel transform.  The parameters N, mu, and q here are the
 * same as for the function fht().  The parameter L is defined (for whatever
 * reason) to be N times the logarithmic spacing of the input array, i.e.
 *   L = N * log(r[N-1]/r[0])/(N-1) */
void compute_u_coefficients(int N, double mu, double q, double L, double kcrc, dcomplex u[]);

double goodkr(int N, double mu, double q, double L, double kr) ;

#endif // FFTLOG_H
