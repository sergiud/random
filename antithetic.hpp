// Copyright (c) Sergiu Dotenco 2010
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @brief Antithetic variates.
 * @file antithetic.hpp
 */

#ifndef ANTITHETIC_HPP
#define ANTITHETIC_HPP

#include <boost/call_traits.hpp>
#include <boost/operators.hpp>
#include <boost/type_traits/remove_reference.hpp>

/**
 * @brief Antithetic variate used for variance reduction.
 *
 * @tparam UniformRandomNumberGenerator The engine that generates uniformly
 *         distributed random numbers.
 */
template<class UniformRandomNumberGenerator>
class Antithetic
    : boost::equality_comparable
      <
          Antithetic<UniformRandomNumberGenerator>,
          Antithetic<UniformRandomNumberGenerator>
      >
{
    bool toggle_;
    UniformRandomNumberGenerator generator_;

    typedef typename
        boost::remove_reference<UniformRandomNumberGenerator>::type Engine;

public:
    typedef Engine::result_type result_type;

private:
    result_type value_;

public:
    Antithetic(typename
               boost::call_traits<UniformRandomNumberGenerator>::param_type u)
        : toggle_(false)
        , generator_(u)
    {
    }

    result_type operator()()
    {
        if (!toggle_) {
            value_ = generator_();
            toggle_ = true;
        }
        else {
            value_ = max() - value_ + min();
            toggle_ = false;
        }

        return value_;
    }

    result_type min() const
    {
        return generator_.min();
    }

    result_type max() const
    {
        return generator_.max();
    }

    bool operator==(const Antithetic& other) const
    {
        return generator_ == generator_ && toggle_ == toggle_ &&
            value_ == value_;
    }
};

/**
 * Returns an antithetic variate for the specified uniform random number
 * generator @a u.
 *
 * @tparam UniformRandomNumberGenerator The engine that generates uniformly
 *         distributed random numbers.
 */
template<class UniformRandomNumberGenerator>
inline Antithetic<UniformRandomNumberGenerator>
    make_antithetic(UniformRandomNumberGenerator& u)
{
    return Antithetic<UniformRandomNumberGenerator>(u);
}

#endif // ANTITHETIC_HPP
