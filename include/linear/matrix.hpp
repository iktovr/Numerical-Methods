#pragma once

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "vector.hpp"

template <class T>
class Matrix {
private:
	std::vector<std::vector<T>> _data;
	size_t _size;

	Matrix() : _data(), _size(0) { }

public:
	Matrix(size_t n) : _data(n, std::vector<T>(n)), _size(n) { }

	std::vector<T>& operator[](size_t index) {
		return _data[index];
	}

	const std::vector<T>& operator[](size_t index) const {
		return _data[index];
	}

	size_t Size() const {
		return _size;
	}

	static Matrix<T> Identity(size_t n) {
		Matrix<T> result(n);
		for (size_t i = 0; i < n; ++i) {
			result[i][i] = 1;
		}
		return result;
	}

	static Matrix<T> Rotation(size_t n, size_t k, size_t l, double angle) {
		Matrix<T> result = Matrix<T>::Identity(n);
		if (k > l) {
			std::swap(k, l);
		}
		double c = std::cos(angle), s = std::sin(angle);
		result[k][k] = c;
		result[l][l] = c;
		result[k][l] = -s;
		result[l][k] = s;
		return result;
	}

	static Matrix<T> Householder(size_t n, const Vector<T>& v) {
		Matrix<T> vv(n);
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < n; ++j) {
				vv[i][j] = v[i] * v[j];
			}
		}

		return Matrix<T>::Identity(n) - vv * (2.0 / (v * v));
	}

	T Norm() const {
		T norm = 0;
		for (const std::vector<T>& row: _data) {
			T sum = 0;
			for (const T& a: row) {
				sum += std::abs(a);
			}
			if (sum > norm) {
				norm = sum;
			}
		}
		return norm;
	}

	Matrix<T> Transpose() {
		Matrix<T> result(*this);
		for (size_t i = 1; i < _size; ++i) {
			for (size_t j = 0; j < i; ++j) {
				std::swap(result[i][j], result[j][i]);
			}
		}
		return result;
	}

	void SwapRows(size_t i, size_t j) {
		if (i >= _size || j >= _size) {
			return;
		}

		std::swap(_data[i], _data[j]);
	}

	void SwapRows(std::vector<std::pair<size_t, size_t>>& P) {
		for (const std::pair<size_t, size_t>& p: P) {
			std::swap(_data[p.first], _data[p.second]);
		}
	}

	void SwapColumns(size_t i, size_t j) {
		if (i >= _size || j >= _size) {
			return;
		}

		for (size_t k = 0; k < _size; ++k) {
			std::swap(_data[k][i], _data[k][j]);
		}
	}

	bool Decompose(Matrix<T>& L, Matrix<T>& U, std::vector<std::pair<size_t, size_t>>& P) const {
		if (_size != L.Size() || _size != U.Size()) {
			throw "Dimension mismatch";
		}
		P.clear();
		for (size_t i = 0; i < _size; ++i) {
			for (size_t j = 0; j < _size; ++j) {
				L[i][j] = 0;
				U[i][j] = _data[i][j];
			}
			L[i][i] = 1;
		}

		for (size_t k = 0; k < _size-1; ++k) {
			T max = U[k][k];
			size_t max_ind = k;
			for (size_t i = k+1; i < _size; ++i) {
				if (U[i][k] > max) {
					max = U[i][k];
					max_ind = i;
				}
			}

			if (max_ind != k) {
				P.push_back(std::make_pair(k, max_ind));
				U.SwapRows(k, max_ind);
				L.SwapRows(k, max_ind);
				L.SwapColumns(k, max_ind);
			}

			for (size_t i = k+1; i < _size; ++i) {
				L[i][k] = U[i][k] / U[k][k];

				for (size_t j = 0; j < _size; ++j) {
					U[i][j] -= L[i][k] * U[k][j];
				}
			}
		}

		return true; 
	}

	template <class U>
	friend std::ostream& operator<<(std::ostream&, const Matrix<U>&);

	template <class U>
	friend std::istream& operator>>(std::istream&, Matrix<U>&);

	template <class U>
	friend Matrix<U> operator*(const Matrix<U>&, const Matrix<U>&);

	template <class U>
	friend Matrix<U> operator+(const Matrix<U>&, const Matrix<U>&);

	template <class U>
	friend Matrix<U> operator-(const Matrix<U>&, const Matrix<U>&);

	template <class U>
	friend Vector<U> operator*(const Matrix<U>&, const Vector<U>&);

	template <class U>
	friend Matrix<U> operator*(const Matrix<U>&, const U&);
};

template <class T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& matrix) {
	for (size_t i = 0; i < matrix.Size(); ++i) {
		for (size_t j = 0; j < matrix.Size(); ++j) {
			os << std::setw(8) << matrix[i][j] << ' ';
		}
		os << '\n';
	}
	return os;
}

template <class T>
std::istream& operator>>(std::istream& is, Matrix<T>& matrix) {
	for (size_t i = 0; i < matrix.Size(); ++i) {
		for (size_t j = 0; j < matrix.Size(); ++j) {
			is >> matrix[i][j];
		}
	}
	return is;
}

template <class T>
Matrix<T> operator*(const Matrix<T>& a, const Matrix<T>& b) {
	if (a.Size() != b.Size()) {
		throw "Dimension mismatch";
	}
	Matrix<T> res(a.Size());
	for (size_t i = 0; i < a.Size(); ++i) {
		for (size_t j = 0; j < a.Size(); ++j) {
			for (size_t k = 0; k < a.Size(); ++k) {
				res[i][j] += a[i][k] * b[k][j];
			}
		}
	}
	return res;
}

template <class T>
Matrix<T> operator+(const Matrix<T>& a, const Matrix<T>& b) {
	if (a.Size() != b.Size()) {
		throw "Dimension mismatch";
	}
	Matrix<T> res(a.Size());
	for (size_t i = 0; i < a.Size(); ++i) {
		for (size_t j = 0; j < a.Size(); ++j) {
			res[i][j] = a[i][j] + b[i][j];
		}
	}
	return res;
}

template <class T>
Matrix<T> operator-(const Matrix<T>& a, const Matrix<T>& b) {
	if (a.Size() != b.Size()) {
		throw "Dimension mismatch";
	}
	Matrix<T> res(a.Size());
	for (size_t i = 0; i < a.Size(); ++i) {
		for (size_t j = 0; j < a.Size(); ++j) {
			res[i][j] = a[i][j] - b[i][j];
		}
	}
	return res;
}

template <class T>
Vector<T> operator*(const Matrix<T>& a, const Vector<T>& b) {
	if (a.Size() != b.size()) {
		throw "Dimension mismatch";
	}
	Vector<T> res(a.Size());
	for (size_t i = 0; i < a.Size(); ++i) {
		for (size_t j = 0; j < a.Size(); ++j) {
			res[i] += a[i][j] * b[j];
		}
	}
	return res;
}

template <class T>
Matrix<T> operator*(const Matrix<T>& a, const T& b) {
	Matrix<T> res(a.Size());
	for (size_t i = 0; i < a.Size(); ++i) {
		for (size_t j = 0; j < a.Size(); ++j) {
			res[i][j] = a[i][j] * b;
		}
	}
	return res;
}
