#pragma once

#include "OdeSolver.h"

/**
 * @brief Euler class
 */
template<typename T>
class Euler : public OdeSolver<T>
{
public:
    Euler() = default;

    Vector<T> calc(T x, T dx, OdeFunction<T>& function) final
    {
        Vector<T> y{function.getParams()};
        Vector<T> dy(y.size());
        Vector<T> dydx{function.derive(x, y)};
        for (size_t i{0U}; i < y.size(); ++i)
        {
            dy[i] = dydx[i] * dx;
        }
        function.setParams(dy);
        return dy;
    }
};
