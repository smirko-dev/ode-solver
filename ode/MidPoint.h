#pragma once

#include "OdeSolver.h"

/**
 * @brief MidPoint class
 */
template<typename T>
class MidPoint : public OdeSolver<T>
{
public:
    MidPoint() = default;

    Vector<T> calc(T x, T dx, OdeFunction<T>& function) final
    {
        Vector<T> y{function.getParams()};
        Vector<T> dy(y.size());
        Vector<T> dydx{function.derive(x, y)};

        Vector<T> k1(y.size());
        for (size_t i{0U}; i < y.size(); ++i)
        {
            k1[i] = dx * dydx[i];
        }
        Vector<T> yt(y.size());
        for (size_t i{0U}; i < y.size(); ++i)
        {
            yt[i] = y[i] + k1[i] / 2.F;
        }

        dydx = function.derive(x + dx / 2.F, yt);
        for (size_t i{0U}; i < y.size(); ++i)
        {
            dy[i] = dx * dydx[i];
        }
        function.setParams(dy);
        return dy;
    }
};
