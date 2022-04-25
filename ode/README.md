# ODE Solver

## OdeFunction

The function for an ODE solver needs to provide the derivative (1st and optional 2nd oder) of an equation.

The parameters are given in a single vector. The setParams and getParams methods implement the mapping of the parameters. 

## OdeSolver

The ODE solver provides an interface for certain implementation´s. 

## Example

Calculate `y(x)=sin(x)` by `dy(x)=cos(x)` using the Euler method. 

```cpp
#include "Euler.h"
#include <cmath>

class Derivative : public OdeFunction<float_t>
{
public:
    Derivative()
        : m_data(1u)
    {
    }

    Vector<float_t> derive(float_t x, [[maybe_unused]] Vector<float_t>& y) final
    {
        return Vector<float_t>{std::cos(x)};
    }
    
    Vector<float_t> getParams() const final
    {
        return m_data;
    }
    
    void setParams(const Vector<float_t>& y) final
    {
        m_data += y;
    }

private:
    Vector<float_t> m_data;
};

int main(int argc, char** argv)
{
    Derivative derivative{};
    Euler<float_t> solver;
    
    static constexpr float_t dt{0.001F};
    for (float_t t{0.0F}; t < 1.F; t += dt)
    {
        solver.calc(t, dt, derivative);
    }
}
```
