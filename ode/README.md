# ODE

## ode::Function

The function for an ODE solver needs to provide the derivative (1st and optional 2nd oder) of an equation with the methods `derive` and `derive2`.

The parameters are given in a single vector. The `setParams` and `getParams` methods implement the mapping of the parameters since the solver just iterates of the given vector.

## ode::Solver

The ODE solver provides an interface for certain implementation´s.

## Example

Solve `f(x) = sin(x)` with `f'(x) = cos(x)` by using the Euler algorithm.

```cpp
#include "ode/Euler.h"
#include <cmath>
#include <iostream>

using Vector = ode::Vector<float_t>;
using Function = ode::Function<float_t>;
using Euler = ode::Euler<float_t>;

class Function : public Function
{
public:
    Function()
        : m_data(1u)
    {
    }

    Vector derive(float_t x, [[maybe_unused]] Vector& y) final
    {
        return Vector{std::cos(x)};
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
};

int main(int argc, char** argv)
{
    static constexpr float_t dt{0.001F};

    Euler euler{};
    Function y{};

    for (float_t t{0.0F}; t < 1.F; t += dt)
    {
        euler.calc(t, dt, y);
        std::cout << t << "," << y.getParams()[0u] << "," << std::sin(x) << std::endl;
    }
}
```
