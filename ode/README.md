# ODE

## ode::Function

The function for an ODE solver needs to provide the derivative (1st and optional 2nd oder) of an equation.

The parameters are given in a single vector. The setParams and getParams methods implement the mapping of the parameters. 

## ode::Solver

The ODE solver provides an interface for certain implementation´s. 

## Example

Calculate `y=sin(x)` by `dy(x)=cos(x)` using the Euler method. 

```cpp
#include "Euler.h"
#include <cmath>

using Vector = ode::Vector<float_t>;
using Function = ode::Function<float_t>;
using Euler = ode::Euler<float_t>;

class Derivative : public Function
{
public:
    Derivative()
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
    Derivative derivative{};
    Euler solver;
    
    static constexpr float_t dt{0.001F};
    for (float_t t{0.0F}; t < 1.F; t += dt)
    {
        solver.calc(t, dt, derivative);
    }
}
```
