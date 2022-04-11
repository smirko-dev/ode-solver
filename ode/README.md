# ODE Solver

## OdeFunction

The function for an ODE solver needs to provide the derivative (1st and optional 2nd oder) of an equation.

The parameters are given in a single vector. The setParams and getParams methods implement the mapping of the parameters. 

## OdeSolver

The ODE solver provides an interface for certain implementation´s. 

## Example

```cpp
#include "Euler.h"
#include <iostream>

class Derivative : public OdeFunction<float_t>
{
public:
    Derivative()
        : m_data(1u)
        , m_solver{}
    {
    }

    void step(const float_t t, const float_t dt)
    {
        m_solver.calc(t, dt, *this);
        std::cout << t << "\t" << m_data[0]<< "\t" << std::sin(t) << std::endl;
    }

protected:
    Vector<float_t> derive(float_t x, Vector<float_t>& y) final
    {
        Vector<float_t> dy(1u);
        dy[0] = std::cos(x);
        return dy;
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
    Euler<float_t> m_solver;
};

int main(int argc, char** argv)
{
    Derivative derivative{};
    static constexpr float_t dt{0.001F};
    for (float_t t{0.0F}; t < 1.F; t += dt)
    {
        derivative.step(t, dt);
    }
}
```
