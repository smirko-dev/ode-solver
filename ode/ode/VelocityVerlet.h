#pragma once

#include "Solver.h"

namespace ode
{
    /// @brief VelocityVerlet solver class
    template<typename T>
    class VelocityVerlet : public Solver<T>
    {
    public:
        using V = typename Solver<T>::Type;

        VelocityVerlet() = default;
        [[nodiscard]] T operator()(const V x, const V dx, const T& y, T& dy, Function<T, V> function)
        {
            T d2y = function(x, y); // a(n)
            T dydx = y + dy * dx + d2y * std::pow(dx, V{2.}) * V{0.5}; // x(n+1) = x(n) + v(n) * dt + 1/2 * a(n) * dt^2
            dy = dy + (function(x + dx, y) + d2y) * dx * V{.5}; // v(n+1) = v(n) + 1/2 * (a(n+1) + a(n)) * dt
            return dydx;
        }
    };
}
