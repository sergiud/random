// Copyright (c) Sergiu Dotenco 2010
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @brief WELL PRNG implementation unit test.
 * @file welltest.cpp
 */

#include <algorithm>

#define BOOST_TEST_MODULE WELL

#include <boost/range/size.hpp>
#include <boost/test/unit_test.hpp>

#include "well.hpp"

/**
 * @brief Generic WELL test case.
 *
 * The test case performs the following checks:
 * -# The last generated value is equal to the value generate by the reference
 *    implementation after @f$10^9@f$ iterations. The generator is seeded using
 *    an array filled with 1s.
 * -# The @c min and @c max methods of the @ref Well generator return 0 and
 *    @f$2^{32}-1@f$ respectively.
 *
 * @tparam RandomNumberGenerator WELL PRNG implementation type.
 * @tparam Expected The expected result after @f$10^9@f$ iterations.
 */
template
<
    class RandomNumberGenerator,
    typename RandomNumberGenerator::result_type Expected
>
class WellTestCase
{
    RandomNumberGenerator rng;

    typedef typename RandomNumberGenerator::result_type result_type;

    result_type generate()
    {
        unsigned state[RandomNumberGenerator::state_size];
        std::uninitialized_fill_n(state, boost::size(state), 1);

        unsigned* p = state;
        rng.seed(p, boost::end(state));

        result_type x;

        int iterations = 1000000000;

        while (iterations-- > 0)
            x = rng();

        return x;
    }

public:
    static void run()
    {
        WellTestCase c;

        BOOST_CHECK_EQUAL(c.generate(), Expected);
        BOOST_CHECK_EQUAL(c.rng.min(), 0U);
        BOOST_CHECK_EQUAL(c.rng.max(), ~0U);
        BOOST_CHECK_EQUAL(c.rng, c.rng);
        BOOST_CHECK(c.rng == c.rng);
    }
};

/**
 * @brief Defines the actual test case.
 *
 * @param name The name of the test case.
 * @param type WELL pseudo-random generator type.
 * @param expected The expected result after @f$10^9@f$ iterations.
 *
 * @hideinitializer
 */
#define DEFINE_TEST_CASE(name, type, expected) \
    BOOST_AUTO_TEST_CASE(name) { WellTestCase<type, expected>::run(); }

DEFINE_TEST_CASE(WELL512a, boost::random::well512a, 0x2b3fe99e)
DEFINE_TEST_CASE(WELL521a, boost::random::well521a, 0xc9878363)
DEFINE_TEST_CASE(WELL521b, boost::random::well521b, 0xb75867f6)
DEFINE_TEST_CASE(WELL607a, boost::random::well607a, 0x7b5043ea)
DEFINE_TEST_CASE(WELL607b, boost::random::well607b, 0xaedee7da)
DEFINE_TEST_CASE(WELL800a, boost::random::well800a, 0x2bfe686f)
DEFINE_TEST_CASE(WELL800b, boost::random::well800b, 0xf009e1bd)
DEFINE_TEST_CASE(WELL1024a, boost::random::well1024a, 0xd07f528c)
DEFINE_TEST_CASE(WELL1024b, boost::random::well1024b, 0x867f7993)
DEFINE_TEST_CASE(WELL19937a, boost::random::well19937a, 0xb33a2cd5)
DEFINE_TEST_CASE(WELL19937b, boost::random::well19937b, 0x191de86a)
DEFINE_TEST_CASE(WELL19937c, boost::random::well19937c, 0x243eaed5)
DEFINE_TEST_CASE(WELL21701a, boost::random::well21701a, 0x7365a269)
DEFINE_TEST_CASE(WELL23209a, boost::random::well23209a, 0x807dacb)
DEFINE_TEST_CASE(WELL23209b, boost::random::well23209b, 0xf1a77751)
DEFINE_TEST_CASE(WELL44497a, boost::random::well44497a, 0xfdd7c07b)
DEFINE_TEST_CASE(WELL44497b, boost::random::well44497b, 0x9406547b)

BOOST_AUTO_TEST_CASE(abc)
{
    typedef boost::random::well512a RandomNumberGenerator;
    RandomNumberGenerator rng;

    unsigned state[RandomNumberGenerator::state_size] = {};
    state[RandomNumberGenerator::state_size - 1] = 1;

    unsigned* p = state;
    rng.seed(p, boost::end(state));

    unsigned x;

    int iterations = 9999;

    while (iterations-- > 0)
        x = rng();

    BOOST_CHECK_EQUAL(x, 0x4df08652U);
}