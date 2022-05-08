#pragma once

#include "Vector.h"
#include <functional>

namespace ode
{
    namespace internal
    {
        template<typename T>
        using is_float = std::is_floating_point<T>;

        template<typename T>
        inline constexpr bool is_float_v = is_float<T>::value;

        template<typename T>
        using is_float_t = typename is_float<T>::type;

        template<typename T> 
        struct is_vector
        {
            constexpr static bool value = false;
        };

        template<typename T> 
        struct is_vector<Vector<T>>
        {
            constexpr static bool value = true;
        };

        template<typename T>
        inline constexpr bool is_vector_v = is_vector<T>::value;

        template<typename T>
        using is_vector_t = typename is_vector<T>::type;
    }

    /// @brief Function for solver; can be derivative of 1st or 2nd order
    template<typename T, typename V>
    using Function = std::function<T(const V x, const T& y)>;
    
    template<typename T, typename Enable = void>
    class Solver;

    /// @brief Solver class specialized for floating point types
    template<typename T>
    class Solver<T, typename std::enable_if<internal::is_float<T>::value>::type>
    {
    public:
        using Type = T;
    };

    /// @brief Solver class specialized for type ode::Vector
    template<typename T>
    class Solver<T, typename std::enable_if<internal::is_vector<T>::value>::type>
    {
    public:
        using Type = typename T::value_type;
    };
}
