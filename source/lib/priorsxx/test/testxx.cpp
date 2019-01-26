#define CATCH_CONFIG_MAIN

#include <array>

#include "catch.hpp"
#include "myrand.hpp"
#include "priors.hpp"
#include "density.hpp"


#define EPSILON_TEST 1e-9

#define M 1e5

namespace priors {

TEST_CASE( "Uniform Prior", "[prior cpp]") {
    
    const int DIM = 4;
    
    spmd_gsl_rand_init(1337);
    
    auto prior = Prior("test_uniform_priors.par");

    double x[DIM] = {-2, 0, 2, 4.5};
    double expected = pow(1.0/12.0,DIM);

    Approx expectedApprox = Approx(expected).margin(EPSILON_TEST);

    double res = prior.eval_pdf(x);
    double logres = prior.eval_logpdf(x);

    REQUIRE( res == expectedApprox );
    
    REQUIRE( exp(logres) == expectedApprox );

    bool ok = true;
    for(int i = 0; i < M; ++i) {    
        double r = prior.rand(i%DIM);
        if( r < -6.0 || r > 6.0) ok = false;;
    }

    REQUIRE( ok );
}


TEST_CASE( "Gaussian Prior", "[prior cpp]") {

    const int DIM = 4;
    
    spmd_gsl_rand_init(1337);
    
    auto prior = Prior("test_gaussian_priors.par");

    double x[DIM] = {-2, 0, 2, 4.5};
    double expected = pow(1.0/12.0,DIM);

    Approx expectedApprox = Approx(expected).margin(EPSILON_TEST);

    double res = prior.eval_pdf(x);
    double logres = prior.eval_logpdf(x);

    REQUIRE( res == expectedApprox );
    
    REQUIRE( exp(logres) == expectedApprox );

    bool ok = true;
    for(int i = 0; i < M; ++i) {    
        double r = prior.rand(i%DIM);
        if( r < -6.0 || r > 6.0) ok = false;;
    }
    REQUIRE( ok );
}

} //namespace priors
