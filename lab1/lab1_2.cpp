#include <iostream>
#include <iomanip>

#include "../include/linear/tridiagonal_matrix.hpp"

int main() {
	std::cout.precision(3);
	std::cout << std::fixed;

	size_t n;
	std::cin >> n;

	TDMatrix<double> a(n);
	Vector<double> b(n);
	std::cin >> a >> b;
	std::cout << "Решение системы:\n" << a.Solve(b) << '\n';
}