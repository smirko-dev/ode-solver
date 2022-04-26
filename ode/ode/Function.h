#pragma once

#include "Vector.h"

namespace ode
{
/**
 * @brief OdeFunction class
 */
template<typename T, typename Enable = void>
class Function;

template<typename T>
class Function<T, typename std::enable_if<std::is_floating_point<T>::value>::type>
{
public:
    /**
     * Calculate derivative
     * @param x      Step variable
     * @param y      List of parameters
     * @return List of calculated parameters
     */
    virtual Vector<T> derive(T x, Vector<T>& y) = 0;

    /**
     * Calculate derivative
     * @param x      Step variable
     * @param y      List of parameters
     * @param dy     List of derived parameters
     * @return List of calculated parameters
     */
    virtual Vector<T> derive2([[maybe_unused]] T x, [[maybe_unused]] Vector<T>& y, [[maybe_unused]] Vector<T>& dy)
    {
        return Vector<T>{};
    }

    /**
     * Return a vector with the parameters to the solver
     */
    virtual Vector<T> getParams() const = 0;

    /**
     * Callback of OdeSolver result parameters
     * @param y      Result parameters
     */
    virtual void setParams(const Vector<T>& y) = 0;
};
}
