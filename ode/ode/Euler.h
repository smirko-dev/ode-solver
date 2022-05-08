#pragma once

#include "Solver.h"

namespace ode
{
    /// @brief Euler solver class
    template<typename T>
    class Euler : public Solver<T>
    {
    public:
        using V = typename Solver<T>::Type;

        Euler() = default;
        [[nodiscard]] T operator()(const V x, const V dx, const T& y, Function<T, V> function)
        {
            return y + function(x, y) * dx;
        }
    };
}
