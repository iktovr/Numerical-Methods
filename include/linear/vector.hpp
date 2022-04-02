#pragma once

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

template <class T>
class Vector {
private:
	std::vector<T> _data;
	size_t _size;

	Vector() : _data(), _size(0) { }

public:
	Vector(size_t n) : _data(n), _size(n) { }

	T& operator[](size_t index) {
		return _data[index];
	}

	const T& operator[](size_t index) const {
		return _data[index];
	}

	size_t Size() const {
		return _size;
	}

	T Norm(int p = 0) {
		T norm = 0;
		if (p > 0) {
			for (const T& a: _data) {
				norm += std::pow(a, p);
			}
			norm = std::pow(norm, 1.0 / p);
		} else {
			for (const T& a: _data) {
				norm = std::max(norm, std::abs(a));
			}
		}
		return norm;
	}

	template <class U>
	friend std::ostream& operator<<(std::ostream&, const Vector<U>&);

	template <class U>
	friend std::istream& operator>>(std::istream&, Vector<U>&);

	template<class U>
	friend Vector<U> operator-(const Vector<U>&, const Vector<U>&);

	template<class U>
	friend Vector<U> operator+(const Vector<U>&, const Vector<U>&);

	template<class U>
	friend Vector<U> operator-(const Vector<U>&);

	template<class U>
	friend U operator*(const Vector<U>&, const Vector<U>&);
};


template <class T>
std::ostream& operator<<(std::ostream& os, const Vector<T>& vec) {
	for (size_t i = 0; i < vec.Size(); ++i) {
		os << std::setw(8) << vec[i] << ' ';
	}
	return os;
}

template <class T>
std::istream& operator>>(std::istream& is, Vector<T>& vec) {
	for (size_t i = 0; i < vec.Size(); ++i) {
		is >> vec[i];
	}
	return is;
}

template<class T>
Vector<T> operator-(const Vector<T>& a, const Vector<T>& b) {
	if (a.Size() != b.Size()) {
		throw "Dimension mismatch";
	}
	Vector<T> res(a);
	for (size_t i = 0; i < b.Size(); ++i) {
		res[i] -= b[i];
	}
	return res;
}

template<class T>
Vector<T> operator+(const Vector<T>& a, const Vector<T>& b) {
	if (a.Size() != b.Size()) {
		throw "Dimension mismatch";
	}
	Vector<T> res(a);
	for (size_t i = 0; i < b.Size(); ++i) {
		res[i] += b[i];
	}
	return res;
}

template<class T>
Vector<T> operator-(const Vector<T>& a) {
	Vector<T> res(a);
	for (size_t i = 0; i < a.Size(); ++i) {
		res[i] *= -1;
	}
	return res;
}

template<class T>
T operator*(const Vector<T>& a, const Vector<T>& b) {
	if (a.Size() != b.Size()) {
		throw "Dimension mismatch";
	}
	T res = 0;
	for (size_t i = 0; i < b.Size(); ++i) {
		res += a[i] * b[i];
	}
	return res;
}