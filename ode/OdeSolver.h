#pragma once

#include "OdeFunction.h"

namespace ode
{
/**
 * @brief OdeSolver class
 */
template<typename T, typename Enable = void>
class Solver;

template<typename T>
class Solver<T, typename std::enable_if<std::is_floating_point<T>::value>::type>
{
public:
    /**
     * Calculate integration step
     * @param x          Variable
     * @param dx         Variable step
     * @param function   Ode function
     * @return calculated parameters
     */
    virtual Vector<T> calc(T x, T dx, Function<T>& function) = 0;

    /**
     * Calculate integration step
     * @param x0         Start variable
     * @param y0         Start parameters
     * @param x          Variable
     * @param dx         Variable step
     * @param function   Ode function
     * @return calculated parameters
     */
    virtual Vector<T> calcRange(T x0, const Vector<T>& y0, T x, T dx, Function<T>& function)
    {
        Vector<T> y{y0};
        for (T t{x0}; t <= x; t += dx)
        {
            y += calc(t, dx, function);
        }
        return y;
    }
};
}
