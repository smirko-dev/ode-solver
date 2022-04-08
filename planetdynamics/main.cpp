
#include "RungeKutta.h"
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
    {
    }
    std::string name{}; //!< Planet name
    Vector<float_t> position{}; //!< Position vector
    Vector<float_t> velocity{}; //!< Velocity vector
    float_t radius{0.F}; //!< Radius
    float_t mass{0.F}; //!< Mass
    std::ofstream file{}; //!< Output stream
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

        // Print results to files
        print();
    }

    void print()
    {
        for (auto& body : m_bodies)
        {
            body.file << body.position[0] << "\t";
            body.file << body.position[1] << "\t";
            body.file << body.position[2] << "\n";

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
                file >> m_bodies[i].name;
                file >> m_bodies[i].position[0];
                file >> m_bodies[i].velocity[1];
                file >> m_bodies[i].mass;
                file >> m_bodies[i].radius;
                std::cout << m_bodies[i].name << " d=" << m_bodies[i].position[0] << " v=" << m_bodies[i].velocity[1] << " m=" << m_bodies[i].mass << " r=" << m_bodies[i].radius << std::endl;
                m_bodies[i].file.open(m_bodies[i].name + ".dat", std::ios::out | std::ios::trunc);
            }
            m_plotfile.open("Solarsystem.dat", std::ios::out | std::ios::trunc);
            return true;
        }
        std::cout << "Invalid file " << filename.c_str() << std::endl;
        return false;
    }

    void finish()
    {
        for (auto& body : m_bodies)
        {
            body.file.close();
        }
        m_plotfile.close();
        std::cout << "Range = [" << m_rangeX[0] << ":" << m_rangeX[1] << ", " << m_rangeY[0] << ":" << m_rangeY[1] << "]" << std::endl;
        std::cout << "Frames = " << m_frames << std::endl;
    }

protected:
    Vector<float_t> derive(float_t x, Vector<float_t>& y) final
    {
        const size_t size = m_bodies.size() * 6U;
        Vector<float_t> dydx(size);

        for (uint32_t a{0U}; a < size; a += 6U)
        {
            for (uint32_t k{0U}; k < 3; ++k)
            {
                // Position
                dydx[a + k] = y[a + k + 3];
                // Velocity
                dydx[a + k + 3] = 0.F;
            }

            for (uint32_t b{0U}; b < size; b += 6)
            {
                if (a != b)
                {
                    // Mass
                    std::vector<Body>::iterator it = m_bodies.begin();
                    it += (b / 6);
                    float_t m = it->mass;
                    // Distance
                    float_t d{0.F};
                    for (uint32_t k{0U}; k < 3; ++k)
                    {
                        d += std::pow(y[a + k] - y[b + k], 2.F);
                    }
                    d = std::sqrt(d);
                    // Velocity
                    for (uint32_t k{0U}; k < 3; ++k)
                    {
                        dydx[a + k + 3] += (y[b + k] - y[a + k]) * m / std::pow(d, 3.F);
                    }
                }
            }
        }
        return dydx;
    }

    Vector<float_t> getParams() const final
    {
        Vector<float_t> y(m_bodies.size() * 6U);

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
        }
        return y;
    }

    void setParams(const Vector<float_t>& y) final
    {
        if (y.size() == m_bodies.size() * 6U)
        {
            uint32_t i{0U};
            for (auto& body : m_bodies)
            {
                // Position
                body.position[0] += y[i++];
                body.position[1] += y[i++];
                body.position[2] += y[i++];

                // Velocity
                body.velocity[0] += y[i++];
                body.velocity[1] += y[i++];
                body.velocity[2] += y[i++];
            }
        }
    }
    
private:
    std::vector<Body> m_bodies{};
    RungeKutta<float_t> m_solver{};
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

            float_t t{0.F};
            static constexpr float_t dt{0.001F};
            while (run.load())
            {
                t += dt;
                world.step(t, dt);
            }
            world.finish();
            run.store(false);
            console_t.join();
        }
    }
    return 0;
}