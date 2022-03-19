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

	Matrix<T> Transpose() {
		Matrix<T> result(*this);
		for (size_t i = 1; i < _size; ++i) {
			for (size_t j = 0; j < i; ++j) {
				std::swap(result[i][j], result[j][i]);
			}
		}
		return result;
	}

	size_t Size() const {
		return _size;
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
Matrix<T> operator+(const Matrix<T>& a, const Matrix<T>& b) {
	if (a.Size() != b.Size()) {
		return Matrix<T>();
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
		return Matrix<T>();
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
		return Vector<T>();
	}
	Matrix<T> res(a.Size());
	for (size_t i = 0; i < a.Size(); ++i) {
		for (size_t j = 0; j < a.Size(); ++j) {
			res[i] = a[i][j] * b[j];
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

template <class T>
T Norm(Vector<T> vec, int p = 0) {
	T norm = 0;
	if (p > 0) {
		for (const T& a: vec) {
			norm += std::pow(a, p);
		}
		norm = std::pow(norm, 1.0 / p);
	} else {
		for (const T& a: vec) {
			norm = std::max(norm, std::abs(a));
		}
	}
	return norm;
}

template<class T>
Vector<T> operator-(const Vector<T>& a, const Vector<T>& b) {
	if (a.size() != b.size()) {
		return Vector<T>();
	}
	Vector<T> res(a);
	for (size_t i = 0; i < b.size(); ++i) {
		res[i] -= b[i];
	}
	return res;
}

template<class T>
Vector<T> operator+(const Vector<T>& a, const Vector<T>& b) {
	if (a.size() != b.size()) {
		return Vector<T>();
	}
	Vector<T> res(a);
	for (size_t i = 0; i < b.size(); ++i) {
		res[i] += b[i];
	}
	return res;
}

template<class T>
Vector<T> operator-(const Vector<T>& a) {
	Vector<T> res(a);
	for (size_t i = 0; i < a.size(); ++i) {
		res[i] *= -1;
	}
	return res;
}

template<class T>
T operator*(const Vector<T>& a, const Vector<T>& b) {
	if (a.size() != b.size()) {
		return T();
	}
	T res = 0;
	for (size_t i = 0; i < b.size(); ++i) {
		res += a[i] * b[i];
	}
	return res;
}