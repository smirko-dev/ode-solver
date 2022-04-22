#pragma once

#include "OdeSolver.h"

/**
 * @brief RungeKutta class
 */
template<typename T>
class RungeKutta : public OdeSolver<T>
{
public:
    RungeKutta() = default;

    Vector<T> calc(T x, T dx, OdeFunction<T>& function) final
    {
        Vector<T> y{function.getParams()};
        Vector<T> dy(y.size());
        Vector<T> yx(y.size());
        Vector<T> k1(y.size());
        Vector<T> k2(y.size());
        Vector<T> k3(y.size());
        Vector<T> k4(y.size());

        Vector<T> dydx{function.derive(x, y)};
        for (size_t i{0U}; i < y.size(); ++i)
        {
            k1[i] = dx * dydx[i];
        }
        for (size_t i{0U}; i < y.size(); ++i)
        {
            yx[i] = y[i] + k1[i] / 2.F;
        }

        dydx = function.derive(x + dx / 2.F, yx);
        for (size_t i{0U}; i < y.size(); ++i)
        {
            k2[i] = dx * dydx[i];
        }
        for (size_t i{0U}; i < y.size(); ++i)
        {
            yx[i] = y[i] + k1[i] / 2.F;
        }

        dydx = function.derive(x + dx / 2.F, yx);
        for (size_t i{0U}; i < y.size(); ++i)
        {
            k3[i] = dx * dydx[i];
        }
        for (size_t i{0U}; i < y.size(); ++i)
        {
            yx[i] = y[i] + k3[i] / 2.F;
        }

        dydx = function.derive(x + dx, yx);
        for (size_t i{0U}; i < y.size(); ++i)
        {
            k4[i] = dx * dydx[i];
        }

        for (size_t i{0U}; i < y.size(); ++i)
        {
            dy[i] = (k1[i] + (k2[i] * 2.F) + (k3[i] * 2.F) + k4[i]) / 6.F;
        }
        function.setParams(dy);
        return dy;
    }
};
