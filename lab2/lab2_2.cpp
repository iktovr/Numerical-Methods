#include <iostream>
#include <iomanip>
#include <functional>
#include <cmath>
#include <vector>
#include <utility>

#include "../include/nonlinear/system_solvers.hpp"

int main() {
	std::cout.precision(8);
	std::cout << std::fixed;

	Vector<function_t<double>> F{
		[](Vector<double> x){ return (x[0] * x[0] + 9) * x[1] - 27; },
		[](Vector<double> x){ return (x[0] - 1.5) * (x[0] - 1.5) + (x[1] - 1.5) * (x[1] - 1.5) - 9; }
	};

	Matrix<function_t<double>> J{
		{[](Vector<double> x){ return 2 * x[0] * x[1]; },  [](Vector<double> x){ return x[0] * x[0] + 9; }},
		{[](Vector<double> x){ return 2 * (x[0] - 1.5); }, [](Vector<double> x){ return 2 * (x[1] - 1.5); }}
	};

	std::vector<Vector<double>> start{{-2, 2}, {4, 1}};

	int k;
	double eps;
	std::cin >> eps >> k;

	size_t iter_count;

	std::cout << "Метод простой итерации:\n";
	for (Vector<double> x: start) {
		iter_count = IterationMethod(x, F, J, eps);

		std::cout << "Решение: " << x << "\nКоличество итераций: " << iter_count << '\n';
	}

	std::cout << '\n';
	std::cout << "Метод Ньютона:\n";
	for (Vector<double> x: start) {
		iter_count = NewtonMethod(x, F, J, eps, k);

		std::cout << "Решение: " << x << "\nКоличество итераций: " << iter_count << '\n';
	}
}