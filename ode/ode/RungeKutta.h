#pragma once

#include "Solver.h"

namespace ode
{
    /// @brief RungeKutta solver class
    template<typename T>
    class RungeKutta : public Solver<T>
    {
    public:
        using V = typename Solver<T>::Type;

        RungeKutta() = default;
        [[nodiscard]] T operator()(const V x, const V dx, const T& y, Function<T, V> function)
        {
            T dydx{function(x, y)};
            T k1 = dydx * dx;
            T yx = y + k1 / V{.5};

            dydx = function(x + dx / V{2.}, yx);
            T k2 = dydx * dx;
            yx = y + k2 / V{2.};

            dydx = function(x + dx / V{2.}, yx);
            T k3 = dydx * dx;
            yx = y + k3 / V{2.};

            dydx = function(x + dx / V{2.}, yx);
            T k4 = dydx * dx;

            T dy = (k1 + (k2 * V{2.}) + (k3 * V{2.}) + k4) / V{6.};
            return y + dy;
        }
    };
}
