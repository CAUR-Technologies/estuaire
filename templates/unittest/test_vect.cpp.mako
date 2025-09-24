#include <nvect.hpp>
#define BOOST_TEST_MODULE nvect
#include <boost/test/unit_test.hpp>

using namespace agsis;

BOOST_AUTO_TEST_SUITE(v)

BOOST_AUTO_TEST_CASE (TestOP)
{
% for op in ['+', '-', '/', '*']
    vect<double, 2> a(1, 2);
    BOOST_CHECK_EQUAL((a ${op} a).x, a.x);
% endfor
}


BOOST_AUTO_TEST_SUITE_E
