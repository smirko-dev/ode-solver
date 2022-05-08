
#include "ode/VelocityVerlet.h"
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <functional>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>

class PlanetSystem
{
public:
    using Type = double_t;
    using Vector = ode::Vector<Type>;
    using Solver = ode::VelocityVerlet<Vector>;

    static constexpr Type Day = 60. * 60. * 24.; // Day in sec
    static constexpr Type Au = 1.496e11; // Astronomical unit in m
    static constexpr Type Sol = 1.989e30; // Mass of the sun in kg
    static constexpr Type G = 6.67408e-11 / (Au * Au * Au) * Sol * (Day * Day); // Gravtitational constant
    static constexpr Type Eps = Sol * (Au * Au) / (Day * Day); // Typical energy scale

    class Body
    {
    public:
        Body() = default;
        std::string name{};
        Type position{}; // Position x value
        Type velocity{}; // Velocity y value
        Type mass{};
    };

public:
    PlanetSystem() = default;

    void step(const Type t, const Type dt)
    {
        using namespace std::placeholders;
        m_pos += m_solver(t, dt, m_pos, m_vel, std::bind(&PlanetSystem::calculate, this, _1, _2));

        Type Ekin{0.0}, Epot{0.0};
        for (size_t i{0u}; i < m_pos.size(); ++i)
        {
            Ekin += m_ekin[i];
            Epot += m_epot[i];
            m_output << m_pos[i * 3u] << ",";
        }
        m_output << Ekin << "," << Epot << std::endl;

        potentialEnergy();
        kineticEnergy();
    }

    void init(const std::string filename, const std::vector<Body>& bodies)
    {
        m_pos.resize(bodies.size() * 3u);
        m_vel.resize(bodies.size() * 3u);
        m_acc.resize(bodies.size() * 3u);
        for (size_t i{0u}; i < bodies.size(); ++i)
        {
            m_pos[i * 3u] = bodies[i].position;
            m_vel[i * 3u + 1u] = bodies[i].velocity;
        }
        m_ekin.resize(bodies.size());
        m_epot.resize(bodies.size());
        m_mass.resize(bodies.size());
        
        potentialEnergy();
        kineticEnergy();
        
        m_output.open(filename, std::ios::out | std::ios::trunc);
    }

    void finish()
    {
        m_output.close();
    }

protected:
    Vector calculate(const Type x, const Vector& y)
    {
        Vector acc(y.size() / 3u);
        for (size_t i{0u}; i < y.size(); ++i)
        {
            for (size_t j{0u}; j < y.size(); ++j)
            {
                if (i != j)
                {
                    Type r3{0.0};
                    Vector rx(y.size() / 3u);
                    for (size_t k{0u}; k < 3; ++k)
                    {
                        rx[k] = y[i * 3u + k] - y[j * 3u + k];
                        r3 += std::pow(rx[k] * rx[k], 1.5);
                    }

                    acc += rx * (-PlanetSystem::G * m_mass[j]) / r3;
                }
            }
        }
        return acc;
    }

    void kineticEnergy()
    {
        for (size_t i{0u}; i < m_ekin.size(); ++i)
        {
            Vector vel{m_vel[i*3u], m_vel[i*3u+1u], m_vel[i*3u+2u]};
            m_ekin[i] = m_mass[i] * vel.norm(2u);
        }
    }

    void potentialEnergy()
    {
        for (size_t i{0u}; i < m_mass.size(); ++i)
        {
            m_epot[i] = 0.0;
            for (size_t j{0u}; j < m_mass.size(); ++i)
            {
                if (i != j)
                {
                    Type r{0.0};
                    Vector rx(3u);
                    for (size_t k{0u}; k < 3; ++k)
                    {
                        rx[k] = m_pos[i * 3u + k] - m_pos[j * 3u + k];
                        r += std::pow(rx[k] * rx[k], 0.5);
                    }
                    m_epot[i] += -PlanetSystem::G * m_mass[i] * m_mass[j] / r / 2.0;
                }
            }
        }
    }

private:
    Vector m_pos{};
    Vector m_vel{};
    Vector m_acc{};
    Vector m_ekin{};
    Vector m_epot{};
    Vector m_mass{};
    Solver m_solver{};
    std::ofstream m_output{};
};

/**
 * Console class
 */
class Console
{
public:
    Console() = default;
    
    void operator()(std::atomic<bool>& run)
    {
        std::string input{};
        while (run.load())
        {
            std::cin >> input;
            if (input == "exit")
            {
                run.store(false);
            }
        }
    }
};

/**
 * main function
 */
int main(int argc, char** argv)
{
    auto m_earth = 5.9742e24 / PlanetSystem::Sol; // Mass of earth in units of Sol
    auto v_earth = 29783.0 / PlanetSystem::Au * PlanetSystem::Day; // Average orbital velocity of the earth in Au/Day

    PlanetSystem planets{};
    planets.init("result.csv", {
        {"Sun", 0.0, 0.0, 1.0},
        {"Earth", -1.0, v_earth, m_earth}
    });

    Console console{};
    std::atomic<bool> run(true);
    std::thread console_t(console, std::ref(run));

    PlanetSystem::Type t{0.};
    static constexpr PlanetSystem::Type dt{0.001};
    while (run.load())
    {
        t += dt;
        planets.step(t, dt);
    }
    planets.finish();
    run.store(false);
    console_t.join();
    return 0;
}