#pragma once

#include <vector>
#include <utility>

#include "matrix.hpp"

template <class T>
struct LUP {
	Matrix<T> L;
	Matrix<T> U;
	std::vector<std::pair<size_t, size_t>> P;

	LUP(const Matrix<T>& matrix): L(matrix.Size()), U(matrix.Size()), P() {
		matrix.Decompose(L, U, P);
	}

	Matrix<T> Compose() {
		return L * U;
	}

	T Det() {
		T res = 1;
		for (size_t i = 0; i < U.Size(); ++i) {
			res *= U[i][i];
		}
		if (P.size() & 1) {
			res = -res;
		}
		return res;
	}

	Vector<T> Solve(Vector<T> b) {
		if (b.size() != L.Size()) {
			return Vector<T>();
		}
		for (const std::pair<size_t, size_t>& p: P) {
			std::swap(b[p.first], b[p.second]);
		}

		Vector<T> x(L.Size()), z(L.Size());
		for (size_t i = 0; i < b.size(); ++i) {
			z[i] = b[i];
			for (size_t j = 0; j < i; ++j) {
				z[i] -= L[i][j] * z[j];
			}
		}

		for (int i = b.size()-1; i >= 0; --i) {
			x[i] = z[i];
			for (int j = b.size()-1; j > i; --j) {
				x[i] -= U[i][j] * x[j];
			}
			x[i] /= U[i][i];
		}

		return x;
	}

	Matrix<T> Invert() {
		Vector<T> e(L.Size());
		Matrix<T> inv(L.Size());

		for (size_t i = 0; i < L.Size(); ++i) {
			e[i] = 1;
			Vector<T> b = Solve(e);
			for (size_t j = 0; j < L.Size(); ++j) {
				inv[j][i] = b[j];
			}
			e[i] = 0;
		}
		
		return inv;
	}
};
