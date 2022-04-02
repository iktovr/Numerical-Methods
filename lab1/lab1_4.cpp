#include <iostream>
#include <iomanip>

#include "../include/linear/matrix.hpp"
#include "../include/linear/rotation_method.hpp"

template <class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
	for (size_t i = 0; i < vec.size(); ++i) {
		os << std::setw(8) << vec[i] << ' ';
	}
	return os;
}

int main() {
	std::cout.precision(3);
	std::cout << std::fixed;

	size_t n;
	std::cin >> n;

	Matrix<double> A(n);
	double eps;
	std::cin >> A >> eps;
	
	std::vector<double> eigenvalues;
	Matrix<double> eigenvectors(n);
	size_t iter_count = RotationMethod(A, eps, eigenvalues, eigenvectors);

	std::cout << "Собственные значения:\n" << eigenvalues << "\n\n";
	std::cout << "Матрица собственных векторов:\n" << eigenvectors << '\n';
	std::cout << "Количество итераций: " << iter_count << '\n';
}