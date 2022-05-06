# ODE

## Example

Solve `f(x) = sin(x)` with `f'(x) = cos(x)` by using the Euler algorithm.

```cpp
#include "ode/Euler.h"
#include <cmath>
#include <iostream>

int main(int argc, char** argv)
{
    static constexpr float_t dt{0.001F};

    float_t y{0.F};
    ode::Euler<float_t> euler{};

    for (float_t t{0.0F}; t < 1.F; t += dt)
    {
        y = euler.calc(t, dt, [](const float_t x, [[maybe_unused]] const float_t& y)
        {
            return std::cos(x);
        });
        std::cout << t << "," << y << "," << std::sin(x) << std::endl;
    }
}
```
