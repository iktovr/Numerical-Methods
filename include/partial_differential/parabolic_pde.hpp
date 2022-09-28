#pragma once

#include <functional>
#include <vector>
#include <tuple>

#include "../linear/tridiagonal_matrix.hpp"
#include "../linear/vector.hpp"

namespace ParabolicPDE {
	template <class T>
	using grid_t = std::vector<std::vector<T>>;

	template <class T>
	std::tuple<std::vector<T>, std::vector<T>, grid_t<T>>
	ExplicitSolver(T a, T b, T c, const std::function<T(T, T)>& f, T start, T end, T t_end, const std::function<T(T)>& psi, 
	               const std::function<T(T)>& phi_start, const std::function<T(T)>& phi_end, int h_count, double sigma) {

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

	template <class T>
	std::tuple<std::vector<T>, std::vector<T>, grid_t<T>>
	ImplicitSolver(T a, T b, T c, const std::function<T(T, T)>& f, T start, T end, T t_end, const std::function<T(T)>& psi, 
	               const std::function<T(T)>& phi_start, const std::function<T(T)>& phi_end, int h_count, double sigma) {

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

		T alpha = (a / h - b / 2) * tau / h,
		  beta = -1 - 2 * a * tau / (h * h) + c * tau,
		  gamma = (a / h + b / 2) * tau / h;

		TDMatrix<T> matrix(h_count+1);
		for (int i = 1; i < h_count; ++i) {
			matrix.a[i] = alpha;
			matrix.b[i] = beta;
			matrix.c[i] = gamma;
		}
		matrix.b[0] = 1;
		matrix.b[h_count] = 1;

		Vector<T> v(h_count+1);
		for (int k = 0; k < tau_count; ++k) {
			for (int i = 1; i < h_count; ++i) {
				v[i] = -u[k][i] - tau * f(x[i], t[k]);
			}
			v[0] = phi_start(t[k+1]);
			v[h_count] = phi_end(t[k+1]);

			u[k+1] = matrix.Solve(v);
		}
		return {t, x, u};
	}

	template <class T>
	std::tuple<std::vector<T>, std::vector<T>, grid_t<T>>
	CombinedSolver(T a, T b, T c, const std::function<T(T, T)>& f, T start, T end, T t_end, const std::function<T(T)>& psi, 
	               const std::function<T(T)>& phi_start, const std::function<T(T)>& phi_end, int h_count, double sigma, double theta) {

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

		T alpha = theta * (a / h - b / 2) * tau / h,
		  beta = -1 + theta * (-2 * a * tau / (h * h) + c * tau),
		  gamma = theta * (a / h + b / 2) * tau / h;

		TDMatrix<T> matrix(h_count+1);
		for (int i = 1; i < h_count; ++i) {
			matrix.a[i] = alpha;
			matrix.b[i] = beta;
			matrix.c[i] = gamma;
		}
		matrix.b[0] = 1;
		matrix.b[h_count] = 1;

		Vector<T> v(h_count+1);
		for (int k = 0; k < tau_count; ++k) {
			for (int i = 1; i < h_count; ++i) {
				double ddu = (u[k][i-1] - 2 * u[k][i] + u[k][i+1]) / (h * h);
				double du = (u[k][i+1] - u[k][i-1]) / (2 * h);
				double uik = a * ddu + b * du + c * u[k][i] + f(x[i], t[k]);
				v[i] = -u[k][i] - tau * (theta * f(x[i], t[k]) + (1 - theta) * uik);
			}
			v[0] = phi_start(t[k+1]);
			v[h_count] = phi_end(t[k+1]);

			u[k+1] = matrix.Solve(v);
		}
		return {t, x, u};
	}
}

