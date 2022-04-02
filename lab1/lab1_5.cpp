#include <iostream>
#include <iomanip>
#include <complex>

#include "../include/linear/matrix.hpp"
#include "../include/linear/qr_decomposition.hpp"

template <class T>
std::ostream& operator<<(std::ostream& os, const std::vector<std::complex<T>>& vec) {
	for (size_t i = 0; i < vec.size(); ++i) {
		if (std::abs(vec[i].imag()) > 0) {
			os << vec[i].real() << std::showpos << vec[i].imag() << std::noshowpos << "*i ";
		} else {
			os << vec[i].real() << ' ';
		}
	}
	return os;
}

int main() {
	std::cout.precision(4);
	std::cout << std::fixed;

	size_t n;
	std::cin >> n;

	Matrix<double> A(n);
	double eps;
	std::cin >> A >> eps;

	std::vector<std::complex<double>> eigenvalues;
	size_t iter_count = QR_Eigenvalues(A, eps, eigenvalues);

	std::cout << "Собственные значения:\n" << eigenvalues << "\n";
	std::cout << "Количество итераций: " << iter_count << '\n';
}