#pragma once

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

template <class T>
using Vector = std::vector<T>;

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

template<class T>
Vector<T> operator-(const Vector<T>& a, const Vector<T>& b) {
	if (a.size() != b.size()) {
		throw "Dimension mismatch";
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
		throw "Dimension mismatch";
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
		throw "Dimension mismatch";
	}
	T res = 0;
	for (size_t i = 0; i < b.size(); ++i) {
		res += a[i] * b[i];
	}
	return res;
}