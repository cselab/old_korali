#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "priors.hpp"
#include "density.hpp"

using namespace priors;

TEST_CASE( "Density", "[density cpp]") {
    auto prior = Prior("test_priors.par");
    prior.print();
    REQUIRE(1==1);
}


TEST_CASE( "Priors", "[prior cpp]") {
    //auto prior = Prior("test_priors.par");
    REQUIRE(1==1);
}


