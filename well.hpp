// boost random/well.hpp header file
//
// Copyright (c) Sergiu Dotenco 2010-2012, 2014
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @brief Implementation of the Well Equidistributed Long-period Linear (WELL)
 * pseudo-random number generator.
 */

#ifndef BOOST_RANDOM_WELL_HPP
#define BOOST_RANDOM_WELL_HPP

#include <cstddef>
#include <iomanip>
#include <istream>
#include <limits>
#include <ostream>

#include <boost/assert.hpp>
#include <boost/config.hpp>
#include <boost/cstdint.hpp>
#include <boost/mpl/apply.hpp>
#include <boost/mpl/bitand.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/minus.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/random/detail/config.hpp>
#include <boost/random/detail/generator_seed_seq.hpp>
#include <boost/random/detail/seed.hpp>
#include <boost/random/detail/seed_impl.hpp>
#include <boost/static_assert.hpp>

namespace boost {
namespace random {

//! @cond hide_private

namespace detail {

template<class UIntType, uint_least32_t N>
struct bitshift_left
{
    static UIntType compute(UIntType a)
    {
        return a << N;
    }
};

template<class UIntType, uint_least32_t N>
struct bitshift_right
{
    static UIntType compute(UIntType a)
    {
        return a >> N;
    }
};

template<class UIntType, int_least32_t N>
struct bitshift
{
    static UIntType compute(UIntType a)
    {
        return mpl::if_c< (N < 0),
            bitshift_left<UIntType, static_cast<uint_least32_t>(-N)>,
            bitshift_right<UIntType, static_cast<uint_least32_t>(N)>
        >::type::compute(a);
    }
};

/**
 * @name Transformation matrices @f$M0,\dotsc,M6@f$ from Table I
 * @{
 */

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
        return bitshift<T, N>::compute(x);
    }
};

template<int N>
struct M3
{
    template<class T>
    static T transform(T x)
    {
        return x ^ bitshift<T, N>::compute(x);
    }
};

template<uint_least32_t a>
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

template<int N, uint_least32_t b>
struct M5
{
    template<class T>
    static T transform(T x)
    {
        return x ^ (bitshift<T, N>::compute(x) & b);
    }
};

template
<
    std::size_t w,
    uint_least32_t q,
    uint_least32_t s,
    uint_least32_t t,
    uint_least32_t a
>
struct M6
{
    template<class T>
    static T transform(T x)
    {
        // Set all bits to 1 except the (s+1)th bit. Ensure only w bits are used.
        BOOST_CONSTEXPR_OR_CONST T ds = (~0U >> (std::numeric_limits<T>::digits - w)) & ~(1 << (w - (s + 1)));

        T result = ((x << q) ^ (x >> (w - q))) & ds;

        if (((x >> (w - (t + 1))) & 1) == 1)
            result ^= a;

        return result;
    }
};

//! @}

/**
 * Conditional expression of type (r & (r - 1)) == 0 which allows to check
 * whether a number @f$r@f$ is of type @f$2^n@f$.
 */
typedef mpl::equal_to<
            mpl::bitand_<
                mpl::_,
                mpl::minus<mpl::_, mpl::int_<1>
            >
        >,
        mpl::int_<0>
    > is_power_of_two;

template<class UIntType, UIntType r>
struct power_of_two_modulo
{
    typedef typename mpl::apply<
            is_power_of_two,
            mpl::integral_c<UIntType, r>
        >::type type;

    BOOST_STATIC_ASSERT_MSG(type::value, "value is not power of two");

    static UIntType compute(UIntType value)
    {
        return value & (r - 1);
    }
};

template<class UIntType, UIntType r>
struct generic_modulo
{
    BOOST_STATIC_ASSERT_MSG(!std::numeric_limits<UIntType>::is_signed,
        "UIntType must be unsigned");

    /**
     * @brief Determines @a value modulo @a r.
     *
     * @pre value >= 0 and value < 2 * r
     * @post value >= 0 and value < r
     */
    static UIntType compute(UIntType value)
    {
        BOOST_ASSERT_MSG(value < 2 * r, "value out of range");

        if (value >= r)
            value -= r;

        return value;
    }
};

template<class UIntType, UIntType r>
struct selective_modulo
{
    typedef typename mpl::apply<
            is_power_of_two,
            mpl::integral_c<UIntType, r>
        >::type r_is_power_of_two;

    static UIntType compute(UIntType value)
    {
        // Use the bitwise AND for power 2 modulo arithmetic, or subtraction
        // otherwise. Subtraction is about two times faster than direct modulo
        // computation.
        return mpl::if_<
                    r_is_power_of_two,
                        power_of_two_modulo<UIntType, r>,
                        generic_modulo<UIntType, r>
                >::type::compute(value);
    }
};

template<uint_least32_t b, uint_least32_t c>
struct matsumoto_kurita_tempering
{
    template<class UIntType, std::size_t r, std::size_t N>
    static UIntType apply(UIntType x, UIntType (&)[N], std::size_t)
    {
        x ^= (x << 7) & b;
        x ^= (x << 15) & c;

        return x;
    }
};

template<uint_least32_t mask>
struct harase_tempering
{
    template<class UIntType, std::size_t r, std::size_t N>
    static UIntType apply(UIntType x, UIntType (&s)[N], std::size_t m2)
    {
        return x ^ (s[selective_modulo<std::size_t, r>::compute(m2 + 1)] & mask);
    }
};

struct no_tempering
{
    template<class UIntType, std::size_t r, std::size_t N>
    static UIntType apply(UIntType x, UIntType (&)[N], std::size_t)
    {
        return x;
    }
};

} // namespace detail

//! @endcond

/**
 * @brief well_engine Equidistributed Long-period Linear (WELL) pseudo-random
 * number generator.
 *
 * The implementation is based on the "Improved Long-Period Generators Based on
 * Linear Recurrences selective_modulo 2" paper by Francois Panneton, Pierre
 * L'Ecuyer and Makoto Matsumoto from ACM Transactions on Mathematical Software,
 * 32 (1, March) 2006, pp. 1-16.
 *
 * @tparam UIntType The unsigned integer type.
 * @tparam w Word size.
 * @tparam r State size.
 * @tparam p Mask bits.
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
    class Tempering // mpl pluggable
>
class well_engine
{
    BOOST_STATIC_ASSERT_MSG(!std::numeric_limits<UIntType>::is_signed,
        "UIntType must be unsigned");
    BOOST_STATIC_ASSERT_MSG(w <=
        static_cast<std::size_t>(std::numeric_limits<UIntType>::digits),
            "word size cannot be represented using UIntType");
    BOOST_STATIC_ASSERT_MSG(r > 0, "state size must be non-zero");
    BOOST_STATIC_ASSERT_MSG(p < w,
        "number of mask bits cannot be greater than word size");
    BOOST_STATIC_ASSERT(m1 > 0 && m1 < r);
    BOOST_STATIC_ASSERT(m2 > 0 && m2 < r);
    BOOST_STATIC_ASSERT(m3 > 0 && m3 < r);

public:
    //! The unsigned integer type.
    typedef UIntType result_type;

    //! Word size.
    BOOST_STATIC_CONSTANT(std::size_t, word_size = w);
    //! State size.
    BOOST_STATIC_CONSTANT(std::size_t, state_size = r);
    //! Number of mask bits.
    BOOST_STATIC_CONSTANT(std::size_t, mask_bits = p);
    //! Default seed value.
    BOOST_STATIC_CONSTANT(UIntType, default_seed = 5489U);
    BOOST_STATIC_CONSTANT(bool, has_fixed_range = false);

    /**
     * Constructs a @c mersenne_twister_engine and calls @c seed().
     */
    well_engine()
    {
        seed();
    }

    /**
     * Constructs a @c mersenne_twister_engine and calls @c seed(value).
     */
    BOOST_RANDOM_DETAIL_ARITHMETIC_CONSTRUCTOR(well_engine, UIntType, value)
    {
        seed(value);
    }

    template<class InputIterator>
    well_engine(InputIterator& first, InputIterator last)
    {
        seed(first, last);
    }

    /**
    * Constructs a well_engine and calls @c seed(gen).
    *
    * @xmlnote
    * The copy constructor will always be preferred over
    * the templated constructor.
    * @endxmlnote
    */
    BOOST_RANDOM_DETAIL_SEED_SEQ_CONSTRUCTOR(well_engine, SeedSeq, seq)
    {
        seed(seq);
    }

    /**
    * Seeds a well_engine using values produced by seq.generate().
    */
    BOOST_RANDOM_DETAIL_SEED_SEQ_SEED(well_engine, SeedSeq, seq)
    {
        detail::seed_array_int<w>(seq, state_);
        index_ = state_size;

        // fix up the state if it's all zeros.
        if ((state_[0] & (~static_cast<UIntType>(0) << mask_bits)) == 0) {
            for (std::size_t j = 1; j != state_size; ++j) {
                if (state_[j] != 0)
                    return;
            }

            state_[0] = static_cast<UIntType>(1) << (w - 1);
        }
    }

    /**
    * Sets the state x(0) to v mod 2w. Then, iteratively,
    * sets x(i) to
    * (i + f * (x(i-1) xor (x(i-1) rshift w-2))) mod 2<sup>w</sup>
    * for i = 1 .. n-1. x(n) is the first value to be returned by operator().
    */
    BOOST_RANDOM_DETAIL_ARITHMETIC_SEED(well_engine, UIntType, value)
    {
        // New seeding algorithm from 
        // http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html
        // In the previous versions, MSBs of the seed affected only MSBs of the
        // state x[].
        const UIntType mask = (max)();

        state_[0] = value & mask;
        
        std::size_t i = 1;
        UIntType (&s)[state_size] = state_;

        for ( ; i != state_size; ++i) {
            // See Knuth "The Art of Computer Programming"
            // Vol. 2, 3rd ed., page 106
            s[i] = static_cast<UIntType>((1812433253U * (s[i - 1] ^ (s[i - 1] >> (w - 2))) + i) & mask);
        }

        index_ = i;
    }

    /**
     * Calls @c seed(default_seed).
     */
    void seed()
    {
        seed(default_seed);
    }

    template<class InputIterator>
    void seed(InputIterator& first, InputIterator last)
    {
        detail::fill_array_int<w>(first, last, state_);
        index_ = state_size;

        // fix up the state if it's all zeros.
        if ((state_[0] & (~static_cast<UIntType>(0) << mask_bits)) == 0) {
            for (std::size_t j = 1; j < state_size; ++j) {
                if (state_[j] != 0)
                    return;
            }

            state_[0] = static_cast<UIntType>(1) << (w - 1);
        }
    }

    /**
     * @brief Generates a random number.
     */
    result_type operator()()
    {
        BOOST_CONSTEXPR_OR_CONST UIntType upper_mask = ~0U << p;
        BOOST_CONSTEXPR_OR_CONST UIntType lower_mask = ~upper_mask;

        // v[i,j] = state[(r-i+j) mod r]
        const std::size_t i = mod(index_);

        // Equivalent to r-i but allows to avoid negative values in the
        // following two expressions
        const std::size_t j = i + r;
        const std::size_t k = mod(j - 1); // [i,r-1]
        const std::size_t l = mod(j - 2); // [i,r-2]

        const std::size_t im1 = mod(i + m1);
        const std::size_t im2 = mod(i + m2);
        const std::size_t im3 = mod(i + m3);

        UIntType z0;
        UIntType z1;
        UIntType z2;
        UIntType z3;
        UIntType z4;

        z0 = (state_[k] & upper_mask) | (state_[l] & lower_mask);
        z1 = T0::transform(state_[i]) ^
             T1::transform(state_[im1]);
        z2 = T2::transform(state_[im2]) ^
             T3::transform(state_[im3]);
        z3 = z1 ^ z2;
        z4 = T4::transform(z0) ^ T5::transform(z1) ^
             T6::transform(z2) ^ T7::transform(z3);

        state_[i] = z3; // v[i+1,1]
        state_[k] = z4; // v[i+1,0]

        index_ = k;

        return Tempering::template apply<UIntType, r>(z4, state_, im2);
    }

    /**
     * Returns the smallest value that the generator can produce.
     */
    static result_type min BOOST_PREVENT_MACRO_SUBSTITUTION ()
    {
        return 0U;
    }

    /**
     * Returns the largest value that the generator can produce.
     */
    static result_type max BOOST_PREVENT_MACRO_SUBSTITUTION ()
    {
        return ~0U >> (std::numeric_limits<UIntType>::digits - w);
    }

    /**
     * Fills a range with random values.
     */
    template<class InputIterator>
    void generate(InputIterator first, InputIterator last)
    {
        detail::generate_from_int(*this, first, last);
    }

    /**
    * Advances the state of the generator by @c z steps.  Equivalent to
    *
    * @code
    * for(unsigned long long i = 0; i < z; ++i) {
    *     gen();
    * }
    * @endcode
    */
    void discard(uintmax_t z)
    {
        while (z-- > 0) {
            (*this)();
        }
    }

    /**
     * @brief Compares the state of two generators for equality.
     */
    friend bool operator==(const well_engine& lhs, const well_engine& rhs)
    {
        for (std::size_t i = 0; i != state_size; ++i)
            if (lhs.compute(i) != rhs.compute(i))
                return false;

        return true;
    }

    /**
     * @brief Compares the state of two generators for inequality.
     */
    friend bool operator!=(const well_engine& lhs, const well_engine& rhs)
    {
        return !(lhs == rhs);
    }

    /**
     * @brief Writes the state to the specified stream.
     */
    template<class E, class T>
    friend std::basic_ostream<E, T>&
        operator<<(std::basic_ostream<E, T>& out, const well_engine& gen)
    {
        const E space = out.widen(' ');

        for (std::size_t i = 0; i != state_size; ++i)
            out << gen.compute(i) << space;

        return out;
    }

    /**
     * @brief Reads the generator state from the specified input stream.
     */
    template<class E, class T>
    friend std::basic_istream<E, T>&
        operator>>(std::basic_istream<E, T>& in, well_engine& gen)
    {
        for (std::size_t i = 0; i != state_size; ++i)
            in >> gen.state_[i] >> std::ws;

        gen.index_ = state_size;

        return in;
    }

private:
    //! @cond show_private

    template<class T>
    static T mod(T value)
    {
        return detail::selective_modulo<T, r>::compute(value);
    }

    UIntType compute(std::size_t index) const
    {
        return state_[mod(index_ + index + r)];
    }

    UIntType state_[r];
    std::size_t index_;

    //! @endcond
};

//! @cond show_private

#ifndef BOOST_NO_INCLASS_MEMBER_INITIALIZATION
template<class UIntType, std::size_t w, std::size_t r, std::size_t p,
    std::size_t m1, std::size_t m2, std::size_t m3, class T0, class T1,
    class T2, class T3, class T4, class T5, class T6, class T7, class Tempering>
const std::size_t well_engine<UIntType, w, r, p, m1, m2, m3, T0, T1, T2, T3, T4, T5,
      T6, T7, Tempering>::word_size;
template<class UIntType, std::size_t w, std::size_t r, std::size_t p,
    std::size_t m1, std::size_t m2, std::size_t m3, class T0, class T1,
    class T2, class T3, class T4, class T5, class T6, class T7, class Tempering>
const std::size_t well_engine<UIntType, w, r, p, m1, m2, m3, T0, T1, T2, T3, T4, T5,
      T6, T7, Tempering>::state_size;
template<class UIntType, std::size_t w, std::size_t r, std::size_t p,
    std::size_t m1, std::size_t m2, std::size_t m3, class T0, class T1,
    class T2, class T3, class T4, class T5, class T6, class T7, class Tempering>
const std::size_t well_engine<UIntType, w, r, p, m1, m2, m3, T0, T1, T2, T3, T4, T5,
      T6, T7, Tempering>::mask_bits;
template<class UIntType, std::size_t w, std::size_t r, std::size_t p,
    std::size_t m1, std::size_t m2, std::size_t m3, class T0, class T1,
    class T2, class T3, class T4, class T5, class T6, class T7, class Tempering>
const UIntType well_engine<UIntType, w, r, p, m1, m2, m3, T0, T1, T2, T3, T4, T5, T6,
      T7, Tempering>::default_seed;
template<class UIntType, std::size_t w, std::size_t r, std::size_t p,
    std::size_t m1, std::size_t m2, std::size_t m3, class T0, class T1,
    class T2, class T3, class T4, class T5, class T6, class T7, class Tempering>
const bool well_engine<UIntType, w, r, p, m1, m2, m3, T0, T1, T2, T3, T4, T5, T6,
      T7, Tempering>::has_fixed_range;
#endif // BOOST_NO_INCLASS_MEMBER_INITIALIZATION

//! @endcond


//! @cond hide_private

namespace detail {

/**
 * @name Base WELL definitions with pluggable tempering method
 * @{
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
    class T7
>
struct well_quoted
{
    template<class T>
    struct apply
    {
        typedef well_engine<UIntType, w, r, p, m1, m2, m3, T0, T1, T2, T3, T4, T5, T6,
            T7, T> type;
    };
};

namespace well_params {

const boost::uint_least32_t a1 = 0xda442d24U;
const boost::uint_least32_t a2 = 0xd3e43ffdU;
const boost::uint_least32_t a3 = 0x8bdcb91eU;
const boost::uint_least32_t a4 = 0x86a9d87eU;
const boost::uint_least32_t a5 = 0xa8c296d1U;
const boost::uint_least32_t a6 = 0x5d6b45ccU;
const boost::uint_least32_t a7 = 0xb729fcecU;

} // namespace well_params

typedef well_quoted<uint32_t, 32, 16, 0, 13, 9, 5,
    M3<-16>, M3<-15>, M3<11>, M0, M3<-2>, M3<-18>, M2<-28>,
    M5<-5, well_params::a1> > well512a_base;

typedef well_quoted<uint32_t, 32, 17, 23, 13, 11, 10,
    M3<-13>, M3<-15>, M1, M2<-21>, M3<-13>, M2<1>, M0, M3<11> >
    well521a_base;

typedef well_quoted<uint32_t, 32, 17, 23, 11, 10, 7,
    M3<-21>, M3<6>, M0, M3<-13>, M3<13>, M2<-10>, M2<-5>, M3<13> >
    well521b_base;

typedef well_quoted<uint32_t, 32, 19, 1, 16, 15, 14,
    M3<19>, M3<11>, M3<-14>, M1, M3<18>, M1, M0, M3<-5> > well607a_base;

typedef well_quoted<uint32_t, 32, 19, 1, 16, 18, 13,
    M3<-18>, M3<-14>, M0, M3<18>, M3<-24>, M3<5>, M3<-1>, M0>
    well607b_base;

typedef well_quoted<uint32_t, 32, 25, 0, 14, 18, 17,
    M1, M3<-15>, M3<10>, M3<-11>, M3<16>, M2<20>, M1, M3<-28> >
    well800a_base;

typedef well_quoted<uint32_t, 32, 25, 0, 9, 4, 22,
    M3<-29>, M2<-14>, M1, M2<19>, M1, M3<10>, M4<well_params::a2>, M3<-25> >
    well800b_base;

typedef well_quoted<uint32_t, 32, 32, 0, 3, 24, 10,
    M1, M3<8>, M3<-19>, M3<-14>, M3<-11>, M3<-7>, M3<-13>, M0>
    well1024a_base;

typedef well_quoted<uint32_t, 32, 32, 0, 22, 25, 26,
    M3<-21>, M3<17>, M4<well_params::a3>, M3<15>, M3<-14>, M3<-21>, M1, M0>
    well1024b_base;

typedef well_quoted<uint32_t, 32, 624, 31, 70, 179, 449,
    M3<-25>, M3<27>, M2<9>, M3<1>, M1, M3<-9>, M3<-21>, M3<21> >
    well19937a_base;

typedef well_quoted<uint32_t, 32, 624, 31, 203, 613, 123,
    M3<7>, M1, M3<12>, M3<-10>, M3<-19>, M2<-11>, M3<4>, M3<-10> >
    well19937b_base;

typedef well_quoted<uint32_t, 32, 679, 27, 151, 327, 84,
    M1, M3<-26>, M3<19>, M0, M3<27>, M3<-11>,
    M6<32, 15, 27, 10, well_params::a4>, M3<-16> >
    well21701a_base;

typedef well_quoted<uint32_t, 32, 726, 23, 667, 43, 462,
    M3<28>, M1, M3<18>, M3<3>, M3<21>, M3<-17>, M3<-28>, M3<-1> >
    well23209a_base;

typedef well_quoted<uint32_t, 32, 726, 23, 610, 175, 662,
    M4<well_params::a5>, M1, M6<32, 15, 15, 30, well_params::a6>,
    M3<-24>, M3<-26>, M1, M0, M3<16> > well23209b_base;

typedef well_quoted<uint32_t, 32, 1391, 15, 23, 481, 229,
    M3<-24>, M3<30>, M3<-10>, M2<-26>, M1, M3<20>,
    M6<32, 9, 5, 14, well_params::a7>, M1> well44497a_base;

//! @}

} // namespace detail

//! @endcond

typedef mpl::apply1<detail::well512a_base,
    detail::no_tempering>::type well512a;
typedef mpl::apply1<detail::well521a_base,
    detail::no_tempering>::type well521a;
typedef mpl::apply1<detail::well521b_base,
    detail::no_tempering>::type well521b;
typedef mpl::apply1<detail::well607a_base,
    detail::no_tempering>::type well607a;
typedef mpl::apply1<detail::well607b_base,
    detail::no_tempering>::type well607b;
typedef mpl::apply1<detail::well800a_base,
    detail::no_tempering>::type well800a;
typedef mpl::apply1<detail::well800b_base,
    detail::no_tempering>::type well800b;
typedef mpl::apply1<detail::well1024a_base,
    detail::no_tempering>::type well1024a;
typedef mpl::apply1<detail::well1024b_base,
    detail::no_tempering>::type well1024b;
typedef mpl::apply1<detail::well19937a_base,
    detail::no_tempering>::type well19937a;
typedef mpl::apply1<detail::well19937b_base,
    detail::no_tempering>::type well19937b;
typedef mpl::apply1<detail::well19937a_base,
    detail::matsumoto_kurita_tempering<0xe46e1700U, 0x9b868000U> >::type well19937c;
typedef mpl::apply1<detail::well21701a_base,
    detail::no_tempering>::type well21701a;
typedef mpl::apply1<detail::well23209a_base,
    detail::no_tempering>::type well23209a;
typedef mpl::apply1<detail::well23209b_base,
    detail::no_tempering>::type well23209b;
typedef mpl::apply1<detail::well44497a_base,
    detail::no_tempering>::type well44497a;
typedef mpl::apply1<detail::well44497a_base,
    detail::matsumoto_kurita_tempering<0x93dd1400U, 0xfa118000U> >::type well44497b;

/**
 * @name Maximally equidistributed WELL versions using Harase's tempering method
 * @{
 */

typedef mpl::apply1<detail::well800a_base,
    detail::harase_tempering<0x4880U> >::type maxeqdist_well800a;
typedef mpl::apply1<detail::well800b_base,
    detail::harase_tempering<0x17030806U> >::type maxeqdist_well800b;
typedef mpl::apply1<detail::well19937a_base,
    detail::harase_tempering<0x4118000U> >::type maxeqdist_well19937a;
typedef mpl::apply1<detail::well19937b_base,
    detail::harase_tempering<0x30200010U> >::type maxeqdist_well19937b;
typedef mpl::apply1<detail::well21701a_base,
    detail::harase_tempering<0x1002U> >::type maxeqdist_well21701a;
typedef mpl::apply1<detail::well23209a_base,
    detail::harase_tempering<0x5100000U> >::type maxeqdist_well23209a;
typedef mpl::apply1<detail::well23209b_base,
    detail::harase_tempering<0x34000300U> >::type maxeqdist_well23209b;
typedef mpl::apply1<detail::well44497a_base,
    detail::harase_tempering<0x48000000U> >::type maxeqdist_well44497a;

//! @}

} // namespace random
} // namespace boost

#endif // !defined(BOOST_RANDOM_WELL_HPP)
