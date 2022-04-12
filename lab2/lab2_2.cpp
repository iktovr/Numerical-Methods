#include <iostream>
#include <iomanip>
#include <functional>
#include <cmath>
#include <vector>
#include <utility>

#include "../include/nonlinear/system_solvers.hpp"

template <class T>
using interval_t = std::pair<Vector<T>, Vector<T>>;

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

	std::vector<interval_t<double>> intervals{{{-2, 2}, {-1, 3}}, {{4, 0.5}, {4.8, 1.2}}};

	int k;
	double eps;
	std::cin >> eps >> k;

	size_t iter_count;
	Vector<double> x(2);

	std::cout << "Метод простой итерации:\n";
	for (const auto& [a, b]: intervals) {
		x = a;
		iter_count = IterationMethod(x, a, b, F, J, eps);

		std::cout << "Решение: " << x << "\nКоличество итераций: " << iter_count << '\n';
	}

	std::cout << '\n';
	std::cout << "Метод Ньютона:\n";
	for (auto [x, _]: intervals) {
		iter_count = NewtonMethod(x, F, J, eps, k);

		std::cout << "Решение: " << x << "\nКоличество итераций: " << iter_count << '\n';
	}
}