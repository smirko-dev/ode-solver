#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <string>

#include "ode/Euler.h"
#include "ode/MidPoint.h"
#include "ode/RungeKutta.h"
#include "ode/VelocityVerlet.h"

// Main funtion
int main(int argc, char** argv)
{
    ode::Euler<float_t> euler1{};
    ode::Euler<ode::Vector<float_t>> euler2{};

    ode::MidPoint<float_t> midpoint1{};
    ode::MidPoint<ode::Vector<float_t>> midpoint2{};

    ode::RungeKutta<float_t> rungekutta1{};
    ode::RungeKutta<ode::Vector<float_t>> rungekutta2{};

    ode::VelocityVerlet<float_t> velocityverlet1{};
    ode::VelocityVerlet<ode::Vector<float_t>> velocityverlet2{};

    float_t y1{0.F}, y2{.0F}, y3{0.}, y4{0.}, dy{1.};
    ode::Vector<float_t> yv{0.F, 1.F}, dyv{0.F, 1.F};

    auto func1 = [](const float_t x, [[maybe_unused]] const float_t& y) -> float_t
    {
            return std::cos(x);
    };
    auto func2 = [](const float_t x, [[maybe_unused]] const ode::Vector<float_t>& y) -> ode::Vector<float_t>
    {
        ode::Vector<float_t> res(y.size());
        for (size_t i{0u}; i < y.size(); ++i)
        {
            res[i] = std::cos(x);
        }
        return res;
    };
    auto func3 = [](const float_t x, [[maybe_unused]] const float_t& y) -> float_t
    {
            return -std::sin(x);
    };
    auto func4 = [](const float_t x, [[maybe_unused]] const ode::Vector<float_t>& y) -> ode::Vector<float_t>
    {
        ode::Vector<float_t> res(y.size());
        for (size_t i{0u}; i < y.size(); ++i)
        {
            res[i] = -std::sin(x);
        }
        return res;
    };

    yv = euler2(.5F, .1F, yv, func2);
    yv = midpoint2(.5F, .1F, yv, func2);
    yv = rungekutta2(.5F, .1F, yv, func2);
    yv = velocityverlet2(.5F, .1F, yv, dyv, func4);

    static constexpr float_t dx{.001F};
    for (float_t x{0.F}; x < 1.F; x+= dx)
    {
        y1 = euler1(x, dx, y1, func1);
        y2 = midpoint1(x, dx, y2, func1);
        y3 = rungekutta1(x, dx, y3, func1);
        y4 = velocityverlet1(x, dx, y4, dy, func3);
        std::cout << std::sin(x) << "\t" << y1 << "\t" << y2 << "\t" << y3 << "\t" << y4 << std::endl;
    }

    return 0;
}
