#include <iostream>
#include <iomanip>
#include <complex>

#include "../include/linear/matrix.hpp"
#include "../include/linear/qr_decomposition.hpp"

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