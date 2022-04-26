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

    Vector derive(float_t x, Vector& y) final
    {
        static constexpr float_t a = 10.F;
        static constexpr float_t b = 28.F;
        static constexpr float_t c = 8.F / 3.F;

        Vector dydx(3u);
        dydx[0u] = a * (y[1u] - y[0u]);
        dydx[1u] = b * y[0u] - y[1u] - y[0u] * y[2u];
        dydx[2u] = y[0u] * y[1u] - c * y[2u];
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
    static constexpr float_t dt{0.05F};

    RungeKutta rk{};
    Lorenz y(dt);

    for (float_t t{0.0F}; t < 2'000.F; t += dt)
    {
        rk.calc(t, dt, y);
        std::cout << y.getParams()[0u]<< "," << y.getParams()[2u] << std::endl;
    }

    return 0;
}
