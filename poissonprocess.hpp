// Copyright (c) Sergiu Dotenco 2010
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @brief Poisson-process implementation.
 * @file poissonprocess.hpp
 */

#ifndef POISSONPROCESS_HPP
#define POISSONPROCESS_HPP

#include <istream>
#include <ostream>

#include <boost/operators.hpp>
#include <boost/random/exponential_distribution.hpp>

/**
 * @brief Homogeneous Poisson process.
 */
template<class RealType = double>
class HomogeneousPoissonProcess
    : boost::equality_comparable
      <
          HomogeneousPoissonProcess<RealType>,
          HomogeneousPoissonProcess<RealType>
      >
{
    typedef boost::exponential_distribution<RealType> exponential_distribution;
    exponential_distribution exp_;

public:
    typedef typename exponential_distribution::result_type result_type;
    typedef result_type param_type;

    explicit HomogeneousPoissonProcess(result_type lambda)
        : exp_(lambda)
    {
    }

    template<class UniformRandomNumberGenerator>
    result_type operator()(UniformRandomNumberGenerator& u, result_type s)
    {          
        return s - exp_(u);
    }

    result_type min() const
    {
        return exp_.min();
    }

    result_type max() const
    {
        return exp_.max();
    }

    result_type lambda() const
    {
        return exp_.lambda();
    }

    void reset()
    {
        exp_.reset();
    }

    bool operator==(const HomogeneousPoissonProcess& other) const
    {
        return exp_ == other.exp_;
    }

    template<class E, class T>
    friend std::basic_ostream<E, T>& 
        operator<<(std::basic_ostream<E, T>& out, 
                   const HomogeneousPoissonProcess& hpp)
    {
        out << hpp.exp_;
        return out;
    }

    template<class E, class T>
    friend std::basic_istream<E, T>& 
        operator>>(std::basic_istream<E, T>& out, 
        HomogeneousPoissonProcess& hpp)
    {
        in >> hpp.exp_;
        return in;
    }
};

#endif // POISSONPROCESS_HPP
