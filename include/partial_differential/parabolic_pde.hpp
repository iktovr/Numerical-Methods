#pragma once

#include <functional>
#include <vector>
#include <tuple>

template <class T>
using grid_t = std::vector<std::vector<T>>;

template <class T>
std::tuple<std::vector<T>, std::vector<T>, grid_t<T>>
ParabolicPDESolver(T a, T b, T c, const std::function<T(T, T)>& f, T start, T end, T t_end, const std::function<T(T)>& psi, 
                   const std::function<T(T)>& phi_start, const std::function<T(T)>& phi_end, int h_count, double sigma) {
	// u_t = a * u_xx + b * u_x + c * u + f
	// x in [start, end]
	// t in [0, t_end]
	// u(x, 0) = psi(x)
	// u(start, t) = phi_start(t)
	// u(end, t) = phi_end(t)
	// h = (end - start) / h_count
	// tau_count = t_end * a * h_count^2 / (end^2 * sigma)

	int tau_count = t_end * a * h_count * h_count / (end * end * sigma);
	grid_t<T> u(tau_count + 1, std::vector<T>(h_count + 1));
	std::vector<T> x(h_count + 1), t(tau_count + 1);
	T h = (end - start) / h_count;
	T tau = t_end / tau_count;

	for (int i = 0; i <= h_count; ++i) {
		x[i] = start + h * i;
	}
	for (int k = 0; k <= tau_count; ++k) {
		t[k] = tau * k;
	}

	for (int i = 0; i <= h_count; ++i) {
		u[0][i] = psi(x[i]);
	}
	for (int k = 1; k <= tau_count; ++k) {
		u[k][0] = phi_start(t[k]);
		u[k][h_count] = phi_end(t[k]);
	}

	for (int k = 0; k < tau_count; ++k) {
		for (int i = 1; i < h_count; ++i) {
			double ddu = (u[k][i-1] - 2 * u[k][i] + u[k][i+1]) / (h * h);
			double du = (u[k][i+1] - u[k][i-1]) / (2 * h);
			u[k+1][i] = (a * ddu + b * du + c * u[k][i] + f(x[i], t[k])) * tau + u[k][i];
		}
	}
	return {t, x, u};
}
