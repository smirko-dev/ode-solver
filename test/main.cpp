#include "Euler.h"
#include "MidPoint.h"
#include "RungeKutta.h"
#include <cmath>
#include <iostream>
#include <string>

// Derivative of a function
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

// Main funtion
int main(int argc, char** argv)
{
    bool silent{false};
    if (argc > 1 && std::string("--silent") == argv[1])
    {
        silent = true;
    }

    Euler<float_t> euler{};
    Derivative y1{};

    MidPoint<float_t> mp{};
    Derivative y2{};

    RungeKutta<float_t> rk{};
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
        auto y = std::sinf(t);
        
        if (!equal(y1.getParams()[0u], y, e))
        {
            errors = true;
            std::cerr << "Mismatch Euler(" << t << ")=" << y1.getParams()[0u] << " != " << y << std::endl;
        }
        if (!equal(y2.getParams()[0u], y, e))
        {
            errors = true;
            std::cerr << "Mismatch MidPoint(" << t << ")=" << y1.getParams()[0u] << " != " << y << std::endl;
        }
        if (!equal(y3.getParams()[0u], y, e))
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
