#pragma once

#include <vector>
#include <iostream>
#include <initializer_list>

template <class T>
class Polynomial {
private:
	std::vector<T> _data;

public:
	Polynomial() : _data(1, T()) { }

	Polynomial(size_t degree, T value = T()) : _data(degree+1, value) { }

	Polynomial(const std::vector<T>& vec) : _data(vec) {
		if (vec.size() == 0) {
			_data.push_back(T());
		}
	}

	Polynomial(std::initializer_list<T> list) : _data(list) {
		if (list.size() == 0) {
			_data.push_back(T());
		}
	}

	size_t Degree() const {
		return _data.size() - 1;
	}

	size_t Size() const {
		return _data.size();
	}

	T& operator[](size_t index) {
		return _data[index];
	}

	const T& operator[](size_t index) const {
		return _data[index];
	}

	T operator()(T x) {
		T res = 0;
		for (size_t i = Degree(); i > 0; --i) {
			res = (res +  _data[i]) * x;
		}
		res += _data[0];
		return res;
	}

	void Assign(size_t degree, T value = T()) {
		_data.assign(degree + 1, value);
	}

	template <class U>
	Polynomial<T>& operator+=(const Polynomial<U>& other) {
		size_t len = std::min(_data.size(), other.Size());
		for (size_t i = 0; i < len; ++i) {
			_data[i] += other[i];
		}
		if (_data.size() < other.Size()) {
			for (size_t i = _data.size(); i < other.Size(); ++i) {
				_data.push_back(other[i]);
			}
		}
		return *this;
	}

	template <class U>
	Polynomial<T>& operator*=(const Polynomial<U>& other) {
		std::vector<T> self = _data;
		Assign(Degree() + other.Degree());
		for (size_t i = 0; i < other.Size(); ++i) {
			for (size_t j = 0; j < self.size(); ++j) {
				_data[j + i] += self[j] * other[i];
			}
		}
		return *this;
	}

	template <class U>
	Polynomial<T>& operator*=(const U& value) {
		for (size_t i = 0; i < _data.size(); ++i) {
			_data[i] *= value;
		}
		return *this;
	}
};

template <class T>
std::ostream& operator<<(std::ostream& os, const Polynomial<T>& p) {
	for (size_t i = p.Degree(); i > 0; --i) {
		if (p[i] < 0) {
			os << "- ";
		} else if (i != p.Degree()) {
			os << "+ ";
		}
		os << std::abs(p[i]) << " * x**" << i << ' ';
	}

	if (p[0] < 0) {
		os << "- ";
	} else if (p.Degree() != 0) {
		os << "+ ";
	}
	os << std::abs(p[0]);
	return os;
}

template<class T, class U>
Polynomial<T> operator*(const Polynomial<T>& a, const Polynomial<U>& b) {
	Polynomial<T> res(a.Degree() + b.Degree());
	for (size_t i = 0; i < a.Size(); ++i) {
		for (size_t j = 0; j < b.Size(); ++j) {
			res[j + i] += b[j] * a[i];
		}
	}
	return res;
}