# ODE

## ode::Function

The function for an ODE solver needs to provide the derivative (1st and optional 2nd oder) of an equation with the methods `derive` and `derive2`.

The parameters are given in a single vector. The `setParams` and `getParams` methods implement the mapping of the parameters since the solver just iterates of the given vector.

## ode::Solver

The ODE solver provides an interface for certain implementation´s.

## Example

Calculate Lorenz attractor by using [Runke Kutta](ode/RungeKutta.h).

```cpp
#include "ode/RungeKutta.h"
#include <cmath>
#include <iostream>

using Vector = ode::Vector<float_t>;
using Function = ode::Function<float_t>;
using RungeKutta = ode::RungeKutta<float_t>;

// Lorenz ODE
class Lorenz : public Function
{
public:
    explicit Lorenz(const float_t dt)
        : m_data{1.F, 0.F, 0.F}
        , m_dt{dt}
    {
    }

    Vector derive(float_t x, [[maybe_unused]] Vector& y) final
    {
        static constexpr float_t sigma = 10.F;
        static constexpr float_t R = 28.F;
        static constexpr float_t b = 8.F / 3.F;

        Vector dydx(3u);
        dydx[0u] = sigma * (y[1u] - y[0u]);
        dydx[1u] = R * y[0u] - y[1u] - y[0u] * y[2u];
        dydx[2u] = y[0u] * y[1u] - b * y[2u];
        return dydx * m_dt;
    }

    Vector getParams() const final
    {
        return m_data;
    }

    void setParams(const Vector& y) final
    {
        m_data += y;
    }

private:
    Vector m_data;
    float_t m_dt;
};

// Main function
int main(int argc, char** argv)
{
    static constexpr float_t dt{0.01F};

    RungeKutta rk{};
    Lorenz y(dt);

    for (float_t t{0.0F}; t < 10'000.F; t += dt)
    {
        rk.calc(t, dt, y);
        // Print out x and z to get the Lorenz butterfly
        std::cout << y.getParams()[0u]<< "," << y.getParams()[2u] << std::endl;
    }
}
```
