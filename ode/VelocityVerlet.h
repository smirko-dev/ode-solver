#pragma once

#include "OdeSolver.h"

/**
 * @brief VelocityVerlet class
 */
template<typename T>
class VelocityVerlet : public OdeSolver<T>
{
public:
    VelocityVerlet() = default;

    Vector<T> calc(T x, T dx, OdeFunction<T>& function) final
    {
        Vector<T> y{function.getParams()};
        Vector<T> dydx{function.derive(x + dx, y)};
        Vector<T> dyd2x{function.derive2(x + dx, y, dydx)};
        function.setParams(dyd2x);
        return dyd2x;
    }
};
