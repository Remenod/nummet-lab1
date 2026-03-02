#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include "tinyexpr.h"

constexpr double eps = 1e-12;

typedef struct
{
    double begin;
    double end;
} range_t;

std::vector<range_t> bracket_roots(
    std::function<double(double)> func,
    range_t range,
    int separation_parts = 10,
    double precision = 0.01)
{
    std::vector<range_t> result = {};
    const double step = std::abs(range.end - range.begin) / separation_parts;

    auto is_derivative_change_sign = [func](range_t range, int separation_parts = 5)
    {
        auto find_derivative = [func](double x, double dx = 0.0005)
        {
            return (func(x + dx) - func(x)) / dx;
        };

        const double step = std::abs(range.end - range.begin) / separation_parts;

        const bool prev_sign = find_derivative(range.begin) >= 0;
        for (int i = 0; i <= separation_parts; i++)
        {
            double d = find_derivative(range.begin + i * step);

            if (std::abs(d) < eps)
                return true;

            bool sign = d >= 0;
            if (sign != prev_sign)
                return true;

            if (std::abs(d) < eps ||
                std::abs(find_derivative(range.begin + (i + 1) * step)) < eps)
                return true;
        }
        return false;
    };

    for (int i = 0; i < separation_parts; i++)
    {
        range_t current_range =
            {
                .begin = range.begin + step * i,
                .end = range.begin + step * (i + 1),
            };

        if (is_derivative_change_sign(current_range) && step > precision)
        {
            auto recursion_result = bracket_roots(func, current_range, separation_parts, precision);
            result.insert(result.end(), recursion_result.begin(), recursion_result.end());
        }
        else
        {
            const double prod = func(current_range.begin) * func(current_range.end);

            const bool sign_change = prod < 0;
            const bool touches_zero = std::abs(prod) < eps;
            const bool no_duplicate = result.empty() || result.back().end != current_range.begin;

            if (sign_change || (touches_zero && no_duplicate))
            {
                result.emplace_back(current_range);
            }
        }
    }

    return result;
}

double refine_roots(
    std::function<double(double)> func,
    range_t range,
    double precision = 0.01)
{
    double fa = func(range.begin);
    double fb = func(range.end);

    if (fa * fb > 0)
        return NAN; // root exesting is not garateed

    while ((range.end - range.begin) > precision)
    {
        double mid = (range.begin + range.end) / 2.0;
        double fm = func(mid);

        if (std::abs(fm) < eps)
            return mid;

        if (fa * fm > 0)
        {
            range.begin = mid;
            fa = fm;
        }
        else
        {
            range.end = mid;
            fb = fm;
        }
    }

    return (range.begin + range.end) / 2.0;
}

int main()
{
    const double precision = 0.01;

    std::cout << "Enter an expression (e.g., sin(x) - 0.5*cos(x^2)):" << std::endl;
    const char *expr;

    std::string line;
    std::getline(std::cin, line);
    expr = line.size() < 1 ? "sin(x) - 0.5*cos(x^2)" : line.c_str();

    double te_x = 0.5;
    te_variable vars[] = {{"x", &te_x, TE_VARIABLE, NULL}};
    int err;
    te_expr *e = te_compile(expr, vars, 1, &err);

    if (err)
    {
        std::string spaces(err - 1, ' ');
        std::cerr << "\033[31m"
                  << "\nAn error occured while parsing expression\n"
                  << line << "\n"
                  << spaces << "^"
                  << "\033[0m"
                  << std::endl;
        return 1;
    }

    auto func = [&e, &te_x](double x)
    {
        te_x = x;
        return te_eval(e);
    };

    std::cout << "Enter operating range" << std::endl;
    range_t range = {};
    std::cin >> range.begin >> range.end;

    std::vector<range_t> bracketed_roots = bracket_roots(func, range, 10, precision);

    auto bracketed_roots_size = bracketed_roots.size();
    std::cout << "Found " << bracketed_roots_size
              << ((bracketed_roots_size == 1) ? " root" : " roots")
              << "\n";
    for (auto el : bracketed_roots)
    {
        auto root = refine_roots(func, el, precision);
        std::cout << root << " on range (" << el.begin << ", " << el.end << ")\n";
    }
    te_free(e);
    return 0;
}
