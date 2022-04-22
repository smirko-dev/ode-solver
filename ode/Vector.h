#pragma once

#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <vector>

/**
 * @brief Floating point comparison
 */
template<typename T, typename std::enable_if<std::is_floating_point<T>::value>* = nullptr>
bool equal(const T a, const T b, const T e = std::numeric_limits<T>::epsilon())
{
    return std::abs(a - b) <= e;
}

/**
 * @brief Vector class
 */
template<typename T, typename Enable = void>
class Vector;

template<typename T>
class Vector<T, typename std::enable_if<std::is_floating_point<T>::value>::type> : public std::vector<T>
{
    using Base = std::vector<T>;

public:
    Vector() = default;
    Vector(const size_t size)
        : Base(size)
    {
    }

    Vector(std::initializer_list<T> rhs)
        : Base{rhs}
    {
    }

    Vector operator+() const
    {
        return Vector(*this);
    }

    Vector operator-() const
    {
        Vector result(Base::size());
        for (size_t i{0U}; i < Base::size(); ++i)
        {
            result[i] = (*this)[i] * -1.F;
        }
        return result;
    }

    Vector operator+(T x) const
    {
        Vector result(Base::size());
        for (size_t i{0U}; i < Base::size(); ++i)
        {
            result[i] = (*this)[i] + x;
        }
        return result;
    }

    Vector operator-(T x) const
    {
        Vector result(Base::size());
        for (size_t i{0U}; i < Base::size(); ++i)
        {
            result[i] = (*this)[i] - x;
        }
        return result;
    }

    Vector operator*(T x) const
    {
        Vector result(Base::size());
        for (size_t i{0U}; i < Base::size(); ++i)
        {
            result[i] = (*this)[i] * x;
        }
        return result;
    }

    Vector operator/(T x) const
    {
        assert(!equal(x, T{0}));
        Vector result(Base::size());
        for (size_t i{0U}; i < Base::size(); ++i)
        {
            result[i] = (*this)[i] / x;
        }
        return result;
    }

    Vector operator+(const Vector& vec) const
    {
        assert(Base::size() == vec.size());
        Vector result(Base::size());
        for (size_t i{0U}; i < Base::size(); ++i)
        {
            result[i] = (*this)[i] + vec[i];
        }
        return result;
    }

    Vector operator-(const Vector& vec) const
    {
        assert(Base::size() == vec.size());
        Vector result(Base::size());
        for (size_t i{0U}; i < Base::size(); ++i)
        {
            result[i] = (*this)[i] - vec[i];
        }
        return result;
    }

    Vector& operator+=(T x)
    {
        *this = *this + x;
        return *this;
    }

    Vector& operator-=(T x)
    {
        *this = *this - x;
        return *this;
    }

    Vector& operator*=(T x)
    {
        *this = *this * x;
        return *this;
    }

    Vector& operator/=(T x)
    {
        *this = *this / x;
        return *this;
    }

    Vector& operator+=(const Vector& vec)
    {
        *this = *this + vec;
        return *this;
    }

    Vector& operator-=(const Vector& vec)
    {
        *this = *this - vec;
        return *this;
    }

    [[nodiscard]] bool isUnity() const
    {
        for (const auto& value : *this)
        {
            if (!equal(value, T{1}))
            {
                return false;
            }
        }
        return true;
    }

    [[nodiscard]] bool isZero() const
    {
        for (const auto& value : *this)
        {
            if (!equal(value, T{0}))
            {
                return false;
            }
        }
        return true;
    }

    Vector& makeZero()
    {
        for (auto& value : *this)
        {
            value = T{0};
        }
        return *this;
    }

    [[nodiscard]] T norm(const uint32_t p) const
    {
        if (0U == p)
        {
            return T{0};
        }
        else if (1U == p)
        {
            T res = T{0};
            for (const auto& value : *this)
            {
                res += std::abs(value);
            }
            return res;
        }
        else if (2U == p)
        {
            T res = T{0};
            for (const auto& value : *this)
            {
                res += value * value;
            }
            return std::sqrt(res);
        }
        else
        {
            T res = T{0};
            for (const auto& value : *this)
            {
                res += std::pow(std::abs(value), static_cast<T>(p));
            }
            return std::pow(res, T{T{1} / p});
        }
        return T{0};
    }

    [[nodiscard]] T length() const
    {
        return norm(2);
    }

    void normalize()
    {
        operator/=(length());
    }

    void normalize(const Vector& vec)
    {
        assert(vec.isZero());
        operator/=(T{1} / vec.length());
    }

    [[nodiscard]] T dot(const Vector& vec) const
    {
        assert(Base::size() == vec.size());
        T result{};
        for (size_t i{0U}; i < Base::size(); ++i)
        {
            result += (*this)[i] * vec[i];
        }
        return result;
    }

};
