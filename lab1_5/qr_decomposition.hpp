#pragma once

#include <cmath>
#include <vector>
#include <complex>

#include "matrix.hpp"

template <class T>
void QR_Decomposition(const Matrix<T>& A, Matrix<T>& Q, Matrix<T>& R) {
	R = A;
	Q = Matrix<T>::Identity(A.Size());
	Vector<T> v(A.Size());
	Matrix<T> H(A.Size());
	for (size_t k = 0; k < A.Size() - 1; ++k) {
		for (size_t i = 0; i < k; ++i) {
			v[i] = 0;
		}
		v[k] = R[k][k] * R[k][k];
		for (size_t i = k+1; i < A.Size(); ++i) {
			v[i] = R[i][k];
			v[k] += R[i][k] * R[i][k];
		}
		v[k] = std::sqrt(v[k]) * (R[k][k] >= 0 ? 1 : -1) + R[k][k];

		H = Matrix<T>::Householder(A.Size(), v);
		R = H * R;
		Q = Q * H;
	}
}

template <class T>
size_t QR_Eigenvalues(Matrix<T> A, double eps, std::vector<std::complex<T>>& eigenvalues) {
	Matrix<T> Q(A.Size()), R(A.Size());
	size_t iter_count = 0;
	Vector<double> eps_1(A.Size());
	Vector<double> eps_2(A.Size());

	do {
		QR_Decomposition(A, Q, R);
		A = R * Q;
		
		for (size_t j = 0; j < A.Size(); ++j) {
			eps_1[j] = 0;
			for (size_t i = j+1; i < A.Size(); ++i) {
				eps_1[j] += A[i][j] * A[i][j]; 
			}
			eps_1[j] = std::sqrt(eps_1[j]);

			eps_2[j] = 0;
			for (size_t i = j+2; i < A.Size(); ++i) {
				eps_2[j] += A[i][j] * A[i][j]; 
			}
			eps_2[j] = std::sqrt(eps_2[j]);
		}
		++iter_count;
	} while (Norm(eps_1, 2) > eps && Norm(eps_2, 2) > eps);

	eigenvalues.assign(A.Size(), std::complex<T>());
	for (size_t i = 0; i < A.Size(); ++i) {
		if (i == A.Size()-1 || std::abs(A[i+1][i]) < eps) {
			eigenvalues[i] = std::complex(A[i][i]);
		} else {
			T d = (A[i][i] + A[i+1][i+1]) * (A[i][i] + A[i+1][i+1]) - 4 * (A[i][i+1] * A[i+1][i] - A[i][i] * A[i+1][i+1]);
			T re = (A[i][i] + A[i+1][i+1]) * 0.5;
			T im = std::sqrt(std::abs(d)) * 0.5;
			eigenvalues[i] = std::complex(re, im);
			eigenvalues[i+1] = std::complex(re, -im);
			++i;
		}
	}

	return iter_count;
}
