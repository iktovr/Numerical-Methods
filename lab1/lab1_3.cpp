#include <iostream>
#include <iomanip>

#include "../include/linear/matrix.hpp"
#include "../include/linear/vector.hpp"
#include "../include/linear/iteration_solvers.hpp"

int main() {
	std::cout.precision(6);
	std::cout << std::fixed;

	size_t n;
	std::cin >> n;

	Matrix<double> A(n);
	Vector<double> b(n);
	double eps;
	std::cin >> A >> b >> eps;

	Vector<double> x(n);
	size_t count = JacobiMethod(A, b, x, eps);
	std::cout << "Решение методом Якоби:\n" << x << '\n';
	std::cout << "Количество итераций: " << count << "\n\n";

	count = JacobiSeidelMethod(A, b, x, eps);
	std::cout << "Решение методом Зейделя:\n" << x << '\n';
	std::cout << "Количество итераций: " << count << '\n';
}