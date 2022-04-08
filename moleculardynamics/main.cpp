
#include "VelocityVerlet.h"
#include <atomic>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>

/**
 * Body class
 */
class Body
{
public:
    Body()
         : position(3U)
         , velocity(3U)
         , force(3U)
    {
    }
    Vector<float_t> position{}; //!< Position vector
    Vector<float_t> velocity{}; //!< Velocity vector
    Vector<float_t> force{}; //!< Force vector
    float_t mass{0.F}; //!< Mass
};

/**
 * Energy class
 */
class Energy
{
public:
    Energy() = default;
    float_t kin{0.F}; //!< Kinetic energy
    float_t pot{0.F}; //!< Potential energy
    float_t total{0.F}; //!< Total energy
    float_t temp{0.F}; //!< Temperature
};

/**
 * World class
 */
class World : public OdeFunction<float_t>
{
public:
    World() = default;

    void step(const float_t t, const float_t dt)
    {
        // Calculate new values
        m_solver.calc(t, dt, *this);

        // Calculate energy
        kineticEnergy();

        // Print results to files
        print();

        // Print energy and temperature
        std::cout << m_energy.total << "\t" << m_energy.temp << std::endl;
    }

    void print()
    {
        for (auto& body : m_bodies)
        {
            m_plotfile << body.position[0] << "\t";
            m_plotfile << body.position[1] << "\t";
            m_plotfile << body.position[2] << "\t";

            if (m_rangeX[0] > body.position[0])
            {
                m_rangeX[0] = body.position[0];
            }
            if (m_rangeX[1] < body.position[0])
            {
                m_rangeX[1] = body.position[0];
            }
            if (m_rangeY[0] > body.position[1])
            {
                m_rangeY[0]= body.position[1];
            }
            if (m_rangeY[1] < body.position[1])
            {
                m_rangeY[1] = body.position[1];
            }
        }
        m_frames++;
        m_plotfile << std::endl;
    }

    bool initialize(const std::string& filename)
    {
        std::ifstream file(filename, std::ios::in);
        if (file.good())
        {
            uint32_t count{0U};
            file >> count;
            std::cout << "Number of bodies = " << count << std::endl;
            m_bodies.resize(count);
            for (uint32_t i{0U}; i < count; ++i)
            {
                file >> m_bodies[i].position[0];
                file >> m_bodies[i].position[1];
                file >> m_bodies[i].position[2];
                file >> m_bodies[i].mass;
            }
            m_plotfile.open("Moleculesystem.dat", std::ios::out | std::ios::trunc);
            return true;
        }
        std::cout << "Invalid file " << filename.c_str() << std::endl;
        return false;
    }

    void finish()
    {
        m_plotfile.close();
        std::cout << "Range = [" << m_rangeX[0] << ":" << m_rangeX[1] << ", " << m_rangeY[0] << ":" << m_rangeY[1] << "]" << std::endl;
        std::cout << "Frames = " << m_frames << std::endl;
    }

protected:
    Vector<float_t> derive(float_t x, Vector<float_t>& y) final
    {
        const size_t size = m_bodies.size() * 9U;
        Vector<float_t> dydx = y;

        // Calculate Lennard Jones Potential
        lennardJones(y);

        // Calculate position / velocity
        std::vector<Body>::iterator it = m_bodies.begin();
        for (uint32_t a{0U}; a < size; a += 9)
        {
            float_t mass{it->mass};
            it++;

            for (uint32_t k{0U}; k < 3; ++k)
            {
                // R += dR * dt + (.5 / m) * F * dt^2
                y[a + k] += y[a + k + 3] * x + .5F / mass * y[a + k + 6] * std::pow(x, 2.F);

                // dR = dR + (.5 / m) * F * dt
                dydx[a + k + 3] = y[a + k + 3] + .5F / mass * y[a + k + 6] * x;
            }
        }

        return dydx;
    }

    Vector<float_t> derive2(float_t x, Vector<float_t>& y, Vector<float_t>& dy) final
    {
        const size_t size = m_bodies.size() * 9U;
        Vector<float_t> dydx = y;

        // Calculate Lennard Jones Potential
        lennardJones(dydx);

        // Calculate velocity
        std::vector<Body>::iterator it = m_bodies.begin();
        for (uint32_t a{0U}; a < size; a += 9)
        {
			float_t mass{it->mass};
            it++;

            for (uint32_t k{0U}; k < 3; ++k)
            {
                // dR = dR + (.5 / m) * F * dt
                dydx[a + k + 3] = dy[a + k + 3] + .5F / mass * dydx[a + k + 6] * x;
            }
        }

        return dydx;
    }

    Vector<float_t> getParams() const final
    {
        Vector<float_t> y(m_bodies.size() * 9);

        uint32_t i{0U};
        for (auto& body : m_bodies)
        {
            // Position
            y[i++] = body.position[0];
            y[i++] = body.position[1];
            y[i++] = body.position[2];

            // Velocity
            y[i++] = body.velocity[0];
            y[i++] = body.velocity[1];
            y[i++] = body.velocity[2];

            // Force
            y[i++] = body.force[0];
            y[i++] = body.force[1];
            y[i++] = body.force[2];
        }
        return y;
    }

    void setParams(const Vector<float_t>& y) final
    {
        if (y.size() == m_bodies.size() * 9)
        {
            uint32_t i{0U};
            for (auto& body : m_bodies)
            {
                // Position
                body.position[0] = y[i++];
                body.position[1] = y[i++];
                body.position[2] = y[i++];

                // Velocity
                body.velocity[0] = y[i++];
                body.velocity[1] = y[i++];
                body.velocity[2] = y[i++];

                // Force
                body.force[0] = y[i++];
                body.force[1] = y[i++];
                body.force[2] = y[i++];
            }
        }
    }
    
    void lennardJones(Vector<float_t>& y)
    {
        const size_t size = m_bodies.size() * 9U;
        Vector<float_t> dr(size);
        float_t pot = 0.F;
        float_t rho = 0.F;

        // Initalize Force
        for (size_t a{0U}; a < size; a += 9)
        {
            for (size_t k{0U}; k < 3; ++k)
            {
                y[a + k + 6] = 0.F;
            }
        }

        // Calculate Lennard Jones Potential
        static constexpr float_t SIGMA_2{40.F * 40.F};
        static constexpr float_t EPSILON{20.F};
        for (size_t a{0U}; a < size; a += 9)
        {
            for (size_t b{a + 9U}; b < size; b += 9)
            {
                for (size_t k{0U}; k < 3; ++k)
                {
                    dr[k] = y[a + k] - y[b + k];
                }
                rho = SIGMA_2 / (std::pow(dr[0], 2.F) + std::pow(dr[1], 2.F) + std::pow(dr[2], 2.F));
                pot = 1.F * (2.F * std::pow(rho, 7.F) - std::pow(rho, 4.F));
                for (size_t k{0U}; k < 3; ++k)
                {
                    y[a + k + 6] += 24.F * EPSILON / SIGMA_2 * pot * dr[k];
                    y[b + k + 6] -= 24.F * EPSILON / SIGMA_2 * pot * dr[k];
                }
                m_energy.pot += 4.F * EPSILON * (std::pow(rho, 6.F) - std::pow(rho, 3.F));
            }
        }
    }

    void kineticEnergy()
    {
        static constexpr float_t KB{1.F};
        m_energy.kin = 0.F;

        for (auto& body : m_bodies)
        {
            // Ekin = 1/2 * m * v^2
            for (uint32_t k{0U}; k < 3; ++k)
            {
                m_energy.kin += .5F * body.mass * (std::pow(body.velocity[k], 2.F));
            }
        }

        const float_t size = static_cast<float_t>(m_bodies.size()) * 9.F;
        m_energy.total = m_energy.kin + m_energy.pot;
        m_energy.temp = m_energy.kin * 2.F / (3.F * size * KB);
    }

private:
    Energy m_energy{};
    std::vector<Body> m_bodies{};
    VelocityVerlet<float_t> m_solver{};
    std::ofstream m_plotfile{};
    size_t m_frames{0U};
    float_t m_rangeX[2];
    float_t m_rangeY[2];
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
    std::string filename{};
    if (argc >= 2)
    {
        World world{};
        if (world.initialize(argv[1]))        
        {
            Console console{};
            std::atomic<bool> run(true);
            std::thread console_t(console, std::ref(run));
            while (run.load())
            {
                world.step(0.0001F, 0.F);
            }
            world.finish();
            run.store(false);
            console_t.join();
        }
    }
    return 0;
}