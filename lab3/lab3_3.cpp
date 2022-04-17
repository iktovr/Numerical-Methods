#include <functional>
#include <vector>
#include <iostream>

#include "../include/function/lsm.hpp"
#include "../include/function/polynomial.hpp"

int main() {
	std::cout.precision(4);
	std::cout << std::fixed;

	std::vector<double> x, y;
	double a, b;
	while (std::cin >> a >> b) {
		x.push_back(a);
		y.push_back(b);
	}

	std::vector<std::function<double(double)>> basis{
		[](double) { return 1; },
		[](double x) { return x; }
	};

	Polynomial<double> p1(LSM(x, y, basis));

	std::cout << "Приближающий многочлен 1-ой степени:\n" << p1 << '\n';
	std::cout << "Ошибка: " << SquareError(std::function<double(double)>{p1}, x, y) << "\n\n";

	basis.push_back(
		[](double x) { return x*x; }
	);

	Polynomial<double> p2(LSM(x, y, basis));

	std::cout << "Приближающий многочлен 2-ой степени:\n" << p2 << '\n';
	std::cout << "Ошибка: " << SquareError(std::function<double(double)>{p2}, x, y) << "\n\n";

	basis.push_back(
		[](double x) { return x*x*x; }
	);

	Polynomial<double> p3(LSM(x, y, basis));

	std::cout << "Приближающий многочлен 3-ой степени:\n" << p3 << '\n';
	std::cout << "Ошибка: " << SquareError(std::function<double(double)>{p3}, x, y) << '\n';
}