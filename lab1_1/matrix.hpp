#pragma once

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

template <class T>
using Vector = std::vector<T>;

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
			return false;
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
	friend Vector<U> operator*(const Matrix<U>&, const Vector<U>&);
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
		return Matrix<T>();
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
Vector<T> operator*(const Matrix<T>& a, const Vector<T>& b) {
	if (a.Size() != b.size()) {
		return Vector<T>();
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
std::ostream& operator<<(std::ostream& os, const Vector<T>& vec) {
	for (size_t i = 0; i < vec.size(); ++i) {
		os << std::setw(8) << vec[i] << ' ';
	}
	return os;
}

template <class T>
std::istream& operator>>(std::istream& is, Vector<T>& vec) {
	for (size_t i = 0; i < vec.size(); ++i) {
		is >> vec[i];
	}
	return is;
}
