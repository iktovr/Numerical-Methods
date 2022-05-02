#include <iostream>
#include <functional>
#include <cmath>

#include "../include/function/integration.hpp"

int main() {
	std::cout.precision(5);
	std::cout << std::fixed;

	std::function<double(double)> f = [](double x){ return x / std::pow(3 * x + 4, 2); };
	double a, b, h;

	std::cin >> a >> b >> h;

	std::cout << "Метод прямоугольников: " << RectangleMethod(f, a, b, h) << '\n' ;
	std::cout << "Погрешность: " << RungeRombergError<double>(RectangleMethod<double>, f, a, b, h, 2, 1) << "\n\n";
	std::cout << "Метод трапеций: " << TrapeziumMethod(f, a, b, h) << '\n';
	std::cout << "Погрешность: " << RungeRombergError<double>(TrapeziumMethod<double>, f, a, b, h, 2, 2) << "\n\n";
	std::cout << "Метод Симпсона: " << SimpsonMethod(f, a, b, h) << '\n';
	std::cout << "Погрешность: " << RungeRombergError<double>(SimpsonMethod<double>, f, a, b, h, 2, 4) << '\n';
}