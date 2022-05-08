#pragma once

#include "Solver.h"

namespace ode
{
    /// @brief MidPoint solver class
    template<typename T>
    class MidPoint : public Solver<T>
    {
    public:
        using V = typename Solver<T>::Type;

        MidPoint() = default;
        [[nodiscard]] T operator()(const V x, const V dx, const T& y, Function<T, V> function)
        {
            T dydx{function(x, y)};

            T k1 = dydx * dx;
            T yt = y + k1 / V{2.};
            dydx = function(x + dx / V{2.}, yt);

            return y + dydx * dx;
        }
    };
}