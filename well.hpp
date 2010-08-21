// Copyright (c) Sergiu Dotenco 2010
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @brief Implementation of the Well-Equidistributed Long-period Linear (WELL)
 *        pseudo-random number generator.
 * @file well.hpp
 */

#ifndef WELL_HPP
#define WELL_HPP

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iomanip>
#include <istream>
#include <limits>
#include <ostream>
#include <stdexcept>

#include <boost/cstdint.hpp>
#include <boost/mpl/apply.hpp>
#include <boost/mpl/bitand.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/minus.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/ref.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

//! @cond hide_private
namespace Detail {
    template<class UIntType, unsigned N>
    struct Left
    {
        static UIntType shift(UIntType a)
        {
            return a << N;
        }
    };

    template<class UIntType, unsigned N>
    struct Right
    {
        static UIntType shift(UIntType a)
        {
            return a >> N;
        }
    };

    template<int N, class UIntType>
    inline UIntType shift(UIntType a)
    {
        return boost::mpl::if_c<(N < 0), 
                    Left<UIntType, -N>, 
                    Right<UIntType, N> 
               >::type::shift(a);
    }

    // Transformation matrices M0,...,M6 from Table I

    struct M0
    {
        template<class T>
        static T transform(T)
        {
            return T(0);
        }
    };

    struct M1
    {
        template<class T>
        static T transform(T x)
        {
            return x;
        }
    };

    template<int N>
    struct M2
    {
        template<class T>
        static T transform(T x)
        {
            return shift<N>(x);
        }
    };

    template<int N>
    struct M3
    {
        template<class T>
        static T transform(T x)
        {
            return x ^ shift<N>(x);
        }
    };

    template<unsigned a>
    struct M4
    {
        template<class T>
        static T transform(T x)
        {
            T result = x >> 1;

            if ((x & 1) == 1)
                result ^= a;

            return result;
        }
    };

    template<int N, unsigned b>
    struct M5
    {
        template<class T>
        static T transform(T x)
        {
            return x ^ (shift<N>(x) & b);
        }
    };

    template<std::size_t w, unsigned q, unsigned a, unsigned ds, unsigned dt>
    struct M6
    {
        template<class T>
        static T transform(T x)
        {
            T result = ((x << q) ^ (x >> (w - q))) & ds;

            if ((x & dt) != 0)
                result ^= a;

            return result;
        }
    };

    /**
     * Conditional expression of type (r & (r - 1)) == 0 which allows to check
     * whether a number r is of type 2^n.
     */
    typedef boost::mpl::equal_to<
                    boost::mpl::bitand_<
                        boost::mpl::_, 
                        boost::mpl::minus<boost::mpl::_, boost::mpl::int_<1> 
                    >
                >,
                boost::mpl::int_<0>
            > IsPowerOfTwo;  

    template<class UIntType, UIntType r>
    struct Power2Modulo
    {
        typedef typename 
            boost::mpl::apply<
                IsPowerOfTwo, 
                boost::mpl::integral_c<UIntType, r> 
            >::type type;

        BOOST_STATIC_ASSERT(type::value);

        template<class T>
        static T calc(T value)
        {
            return value & (r - 1);
        }
    };

    template<class UIntType, UIntType r>
    struct GenericModulo
    {
        /**
         * Determines @a value modulo @a r.
         * 
         * @pre value >= 0 and value < 2 * r
         * @post value >= 0 and value < r
         */
        template<class T>
        static T calc(T value)
        {
            BOOST_STATIC_ASSERT(boost::is_unsigned<T>::value);
            assert(value < 2 * r);

            if (value >= r)
                value -= r;

            return value;
        }
    };

    template<unsigned b, unsigned c>
    struct MatsumotoKuritaTempering
    {
        template<class UIntType>
        static UIntType apply(UIntType x)
        {
            x ^= (x << 7) & b;
            x ^= (x << 15) & c;

            return x;
        }
    };

    struct NoTempering
    {
        template<class UIntType>
        static UIntType apply(UIntType x)
        {
            return x;
        }
    };
} // namespace Detail
//! @endcond hide_private

/**
 * @brief Well-Equidistributed Long-period Linear (WELL) pseudo-random number
 * generator.
 *
 * The WELL pseudo-random number generator has been characterized in "Improved
 * Long-Period Generators Based on Linear Recurrences Modulo 2", by Francois
 * Panneton, Pierre L'Ecuyer and Makoto Matsumoto from ACM Transactions on
 * Mathematical Software, 32 (1, March) 2006, 1-16.
 *
 * @tparam UIntType The unsigned integer type.
 * @tparam w Word size.
 * @tparam r State size.
 */
template
<
    class UIntType,
    std::size_t w,
    std::size_t r,
    std::size_t p,
    std::size_t m1,  
    std::size_t m2,
    std::size_t m3,
    class T0,
    class T1,
    class T2,
    class T3,
    class T4,
    class T5,
    class T6,
    class T7,
    class Tempering = Detail::NoTempering
>
class Well
{
    BOOST_STATIC_ASSERT(std::numeric_limits<UIntType>::digits == w);

    UIntType state_[r];
    std::size_t index_;
    
    typedef typename 
            boost::mpl::apply<
                Detail::IsPowerOfTwo, 
                boost::mpl::integral_c<UIntType, r> 
            >::type is_power_of_two;

    template<class T>
    static T mod(T value)
    {
        // Use the bitwise and for power 2 modulo arithmetic or subtraction
        // otherwise. Subtraction is about two times faster than direct modulo
        // calculation.
        return boost::mpl::if_<
                    is_power_of_two,
                        Detail::Power2Modulo<UIntType, r>,
                        Detail::GenericModulo<UIntType, r>
                >::type::calc(value);
    }

    UIntType state(std::size_t index) const
    {
        return state_[mod(index)];
    }

    UIntType compute(std::size_t index) const
    {
        return state_[(index_ + index + r) % r];
    }

public:
    typedef UIntType result_type;

    static const std::size_t word_size = w;
    static const std::size_t state_size = r;
    static const std::size_t mask_bits = p;
    static const result_type default_seed = 5489;

    explicit Well(result_type value = default_seed)
    {
        seed(value);
    }

    template<class InputIterator>
    Well(InputIterator& first, InputIterator last)
    {
        seed(first, last);
    }
    
    template<class Generator>
    explicit Well(Generator& g)
    {
        seed(g);
    }

    template<class Generator>
    void seed(Generator& g)
    {
        // Ensure std::generate_n doesn't copy the generator g by using 
        // boost::reference_wrapper
        std::generate_n(state_, state_size, boost::ref(g));
    }

    void seed(UIntType value)
    {
        if (value == 0)
            value = default_seed;

        state_[0] = value;

        std::size_t i = 1;
        UIntType *const s = state_;

        // Same generator used to seed Mersenne twister
        for ( ; i < state_size; ++i)
            s[i] = (1812433253U * (s[i - 1] ^ (s[i - 1] >> (w - 2))) + i);

        index_ = i;
    }

    void seed()
    {
        seed(default_seed);
    }

    template<class InputIterator>
    void seed(InputIterator& first, InputIterator last)
    {
        index_ = 0;
        std::size_t i = 0;

        for ( ; i < r && first != last; ++i, ++first)
            state_[i] = *first;

        if (first == last && i < state_size)
            throw std::invalid_argument("Seed sequence too short");
    }

    result_type operator()()
    {
        const UIntType upper_mask = ~0U << p;
        const UIntType lower_mask = ~upper_mask;

        // v[i,j] = state[(r-i+j) mod r]
        std::size_t i = index_;
        // Equivalent to r-i but allows to avoid negative values in the
        // following two expressions
        std::size_t j = i + r;
        std::size_t k = mod(j - 1); // [i,r-1]
        std::size_t l = mod(j - 2); // [i,r-2]

        UIntType z0, z1, z2, z3, z4;

        z0 = (state_[k] & upper_mask) | (state_[l] & lower_mask);
        z1 = T0::transform(state_[i]) ^ 
             T1::transform(state(i + m1));
        z2 = T2::transform(state(i + m2)) ^ 
             T3::transform(state(i + m3));
        z3 = z1 ^ z2;
        z4 = T4::transform(z0) ^ T5::transform(z1) ^ 
             T6::transform(z2) ^ T7::transform(z3);

        state_[i] = z3; // v[i+1,1]
        state_[k] = z4; // v[i+1,0]

        index_ = k;

        return Tempering::apply(z4);
    }

    result_type min() const
    {
        return 0;
    }

    result_type max() const
    {
        return ~0U >> (std::numeric_limits<UIntType>::digits - w);
    }

    friend bool operator==(const Well& lhs, const Well& rhs)
    {
        for (std::size_t i = 0; i < state_size; ++i)
            if (lhs.compute(i) != rhs.compute(i))
                return false;

        return true;
    }

    friend bool operator!=(const Well& lhs, const Well& rhs)
    {
        return !(lhs == rhs);
    }

    template<class E, class T>
    friend std::basic_ostream<E, T>& operator<<(std::basic_ostream<E, T>& out,
                                                const Well& well)
    {
        E space = out.widen(' ');

        for (std::size_t i = 0; i < state_size; ++i)
            out << well.compute(i) << space;

        return out;
    }

    template<class E, class T>
    friend std::basic_istream<E, T>& operator>>(std::basic_istream<E, T>& in,
                                                Well& well)
    {
        for (std::size_t i = 0; i < state_size; ++i)
            in >> well.state_[i] >> std::ws;

        well.index_ = state_size;

        return in;
    }
};

typedef Well<boost::uint32_t, 32, 16, 0, 13, 9, 5, 
    Detail::M3<-16>, Detail::M3<-15>, Detail::M3<11>, Detail::M0, 
    Detail::M3<-2>, Detail::M3<-18>, Detail::M2<-28>, // FIXME: M2 not M3?
    Detail::M5<-5, 0xda442d24> > Well512a;

typedef Well<boost::uint32_t, 32, 17, 23, 13, 11, 10, 
    Detail::M3<-13>, Detail::M3<-15>, Detail::M1, Detail::M2<-21>, 
    Detail::M3<-13>, Detail::M2<1>, Detail::M0, Detail::M3<11> > Well521a;

typedef Well<boost::uint32_t, 32, 17, 23, 11, 10, 7, 
    Detail::M3<-21>, Detail::M3<6>, Detail::M0, Detail::M3<-13>, 
    Detail::M3<13>, Detail::M2<-10>, Detail::M2<-5>, Detail::M3<13> > Well521b;

typedef Well<boost::uint32_t, 32, 19, 1, 16, 15, 14, 
    Detail::M3<19>, Detail::M3<11>, Detail::M3<-14>, Detail::M1, 
    Detail::M3<18>, Detail::M1, Detail::M0, Detail::M3<-5> > Well607a;

typedef Well<boost::uint32_t, 32, 19, 1, 16, 18, 13, 
    Detail::M3<-18>, Detail::M3<-14>, Detail::M0, Detail::M3<18>, 
    Detail::M3<-24>, Detail::M3<5>, Detail::M3<-1>, Detail::M0> Well607b;

typedef Well<boost::uint32_t, 32, 25, 0, 14, 18, 17, 
    Detail::M1, Detail::M3<-15>, Detail::M3<10>, Detail::M3<-11>, 
    Detail::M3<16>, Detail::M2<20>, Detail::M1, Detail::M3<-28> > Well800a;

typedef Well<boost::uint32_t, 32, 25, 0, 9, 4, 22, 
    Detail::M3<-29>, Detail::M2<-14>, Detail::M1, Detail::M2<19>, 
    Detail::M1, Detail::M3<10>, Detail::M4<0xd3e43ffd>, Detail::M3<-25> > 
    Well800b;

typedef Well<boost::uint32_t, 32, 32, 0, 3, 24, 10, 
    Detail::M1, Detail::M3<8>, Detail::M3<-19>, Detail::M3<-14>, 
    Detail::M3<-11>, Detail::M3<-7>, Detail::M3<-13>, Detail::M0> Well1024a;

typedef Well<boost::uint32_t, 32, 32, 0, 22, 25, 26, 
    Detail::M3<-21>, Detail::M3<17>, Detail::M4<0x8bdcb91e>, Detail::M3<15>, 
    Detail::M3<-14>, Detail::M3<-21>, Detail::M1, Detail::M0> Well1024b;

typedef Well<boost::uint32_t, 32, 624, 31, 70, 179, 449, 
    Detail::M3<-25>, Detail::M3<27>, Detail::M2<9>, Detail::M3<1>, 
    Detail::M1, Detail::M3<-9>, Detail::M3<-21>, Detail::M3<21> > Well19937a;

typedef Well<boost::uint32_t, 32, 624, 31, 203, 613, 123, 
    Detail::M3<7>, Detail::M1, Detail::M3<12>, Detail::M3<-10>, 
    Detail::M3<-19>, Detail::M2<-11>, Detail::M3<4>, Detail::M3<-10> > 
    Well19937b;

typedef Well<boost::uint32_t, 32, 624, 31, 70, 179, 449, 
    Detail::M3<-25>, Detail::M3<27>, Detail::M2<9>, Detail::M3<1>, 
    Detail::M1, Detail::M3<-9>, Detail::M3<-21>, Detail::M3<21>,
    Detail::MatsumotoKuritaTempering<0xe46e1700, 0x9b868000> > Well19937c;

typedef Well<boost::uint32_t, 32, 679, 27, 151, 327, 84, 
    Detail::M1, Detail::M3<-26>, Detail::M3<19>, Detail::M0, 
    Detail::M3<27>, Detail::M3<-11>, 
    Detail::M6<32, 15, 0x86a9d87e, 0xffffffef, 0x00200000>, Detail::M3<-16> > 
    Well21701a;

typedef Well<boost::uint32_t, 32, 726, 23, 667, 43, 462, 
    Detail::M3<28>, Detail::M1, Detail::M3<18>, Detail::M3<3>, 
    Detail::M3<21>, Detail::M3<-17>, Detail::M3<-28>, Detail::M3<-1> > 
    Well23209a;

typedef Well<boost::uint32_t, 32, 726, 23, 610, 175, 662, 
    Detail::M4<0xa8c296d1>, Detail::M1, 
    Detail::M6<32, 15, 0x5d6b45cc, 0xfffeffff, 0x00000002>, 
    Detail::M3<-24>, Detail::M3<-26>, Detail::M1, Detail::M0, Detail::M3<16> > 
    Well23209b;

typedef Well<boost::uint32_t, 32, 1391, 15, 23, 481, 229, 
    Detail::M3<-24>, Detail::M3<30>, Detail::M3<-10>, Detail::M2<-26>, 
    Detail::M1, Detail::M3<20>, 
    Detail::M6<32, 9, 0xb729fcec, 0xfbffffff, 0x00020000>, Detail::M1> 
    Well44497a;

typedef Well<boost::uint32_t, 32, 1391, 15, 23, 481, 229, 
    Detail::M3<-24>, Detail::M3<30>, Detail::M3<-10>, Detail::M2<-26>, 
    Detail::M1, Detail::M3<20>, 
    Detail::M6<32, 9, 0xb729fcec, 0xfbffffff, 0x00020000>, Detail::M1,
    Detail::MatsumotoKuritaTempering<0x93dd1400, 0xfa118000> > Well44497b;

#endif // WELL_HPP
