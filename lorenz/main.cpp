#include "ode/RungeKutta.h"
#include <cmath>
#include <iostream>

int main(int argc, char** argv)
{
    static constexpr float_t dt{0.05F};

    using Type = ode::Vector<float_t>;

    Type y{1.F, 0.F, 0.F};
    ode::RungeKutta<Type> rk{};

    for (float_t t{0.0F}; t < 2'000.F; t += dt)
    {
        y = rk(t, dt, y, [](float_t x, const Type& y) -> Type
        {
            static constexpr float_t a = 10.F;
            static constexpr float_t b = 28.F;
            static constexpr float_t c = 8.F / 3.F;

            Type dydx(3u);
            dydx[0u] = a * (y[1u] - y[0u]);
            dydx[1u] = b * y[0u] - y[1u] - y[0u] * y[2u];
            dydx[2u] = y[0u] * y[1u] - c * y[2u];
            return dydx * dt;
        });
        std::cout << y[0u] << "," << y[2u] << std::endl;
    }

    return 0;
}
