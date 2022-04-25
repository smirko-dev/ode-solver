#include "Euler.h"
#include "MidPoint.h"
#include "RungeKutta.h"
#include <cmath>
#include <iostream>
#include <string>

using Vector = ode::Vector<float_t>;
using Function = ode::Function<float_t>;
using Euler = ode::Euler<float_t>;
using MidPoint = ode::MidPoint<float_t>;
using RungeKutta = ode::RungeKutta<float_t>;

// Derivative of a function
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

// Main funtion
int main(int argc, char** argv)
{
    bool silent{false};
    if (argc > 1 && std::string("--silent") == argv[1])
    {
        silent = true;
    }

    Euler euler{};
    Derivative y1{};

    MidPoint mp{};
    Derivative y2{};

    RungeKutta rk{};
    Derivative y3{};

    // Calculate and print results in range [0..1]
    bool errors{false};
    static constexpr float_t dt{0.001F};
    static constexpr float_t e{0.001F};
    for (float_t t{0.0F}; t < 1.F; t += dt)
    {
        euler.calc(t, dt, y1);
        mp.calc(t, dt, y2);
        rk.calc(t, dt, y3);
        auto y = sinf(t);
        
        if (!ode::equal(y1.getParams()[0u], y, e))
        {
            errors = true;
            std::cerr << "Mismatch Euler(" << t << ")=" << y1.getParams()[0u] << " != " << y << std::endl;
        }
        if (!ode::equal(y2.getParams()[0u], y, e))
        {
            errors = true;
            std::cerr << "Mismatch MidPoint(" << t << ")=" << y1.getParams()[0u] << " != " << y << std::endl;
        }
        if (!ode::equal(y3.getParams()[0u], y, e))
        {
            errors = true;
            std::cerr << "Mismatch RungeKutta(" << t << ")=" << y1.getParams()[0u] << " != " << y << std::endl;
        }

        if (!silent)
        {
            std::cout << t << "\t" << y1.getParams()[0u]<< "\t" << y2.getParams()[0u]<< "\t" << y3.getParams()[0u] << "\t" << std::sin(t) << std::endl;
        }
    }

    return (errors ? -1 : 0);
}
