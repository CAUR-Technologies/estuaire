#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE nvect
#include <boost/test/unit_test.hpp>
#include <nvect.hpp>

using namespace agsis;

BOOST_AUTO_TEST_SUITE(v)

#define DOUBLE_EPS 1e-10

BOOST_AUTO_TEST_CASE (TestPlusOperation)
{

    vect<double, 2> a2(1, 2);
    vect<double, 2> b2(0, 0);
    vect<double, 2> c2(-1, 5);

    vect<double, 2> c3(-1, 5);

    // NORM VERIFICATION
    BOOST_CHECK_CLOSE(norm(a2), 2.2360679774997898, DOUBLE_EPS);
    BOOST_CHECK_CLOSE(norm(b2), 0, DOUBLE_EPS);
    BOOST_CHECK_CLOSE(norm(c3), 5.0990195135927845, DOUBLE_EPS);
}

BOOST_AUTO_TEST_SUITE_END()
