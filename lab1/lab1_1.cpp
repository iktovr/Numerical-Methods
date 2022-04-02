#include <iostream>
#include <iomanip>

#include "../include/linear/matrix.hpp"
#include "../include/linear/vector.hpp"
#include "../include/linear/lup.hpp"

int main() {
	std::cout.precision(3);
	std::cout << std::fixed;

	size_t n;
	std::cin >> n;

	Matrix<double> a(n);
	Vector<double> b(n);
	std::cin >> a >> b;

	LUP<double> lup(a);

	std::cout << "Решение системы:\n" << lup.Solve(b) << '\n';
	std::cout << "Обратная матрица системы:\n" <<  lup.Invert() << '\n';
	std::cout << "Определитель матрицы системы:\n" << lup.Det() << '\n';
}