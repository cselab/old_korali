#define CATCH_CONFIG_MAIN
#include <vector>

#include "catch.hpp"
#include "engine_cmaes.hpp"

extern "C" {
//    #include "torc.h"
    #include "fitfun_tests.h"
}

#define EPSILON_TEST 1e-5

void compare_best_with_ref(CmaesEngine* engine, std::vector<double> refBestEver, 
                            double refFunVal, double eps);


TEST_CASE( "Ackley", "[engine cmaes cpp]") {
    
    fitfun_initialize_simple("Ackley");
    
    auto engine = CmaesEngine(&fitfun, "./ackley/");
    engine.run();

    compare_best_with_ref( &engine, 
                           { -2.00967e-12,  -7.28353e-13, -3.85381e-13, -2.02323e-14 }, 
                           9.93963, 
                           EPSILON_TEST);
}


TEST_CASE( "Dixon Price", "[engine cmaes cpp]") {
    
    fitfun_initialize_simple("Dixon_Price");
    
    auto engine = CmaesEngine( &fitfun, "./dixon_price/");
    engine.run();

    compare_best_with_ref( &engine, 
                           { 1, 0.707107, 0.594604, -0.545254 },
                           9.93963, 
                           EPSILON_TEST);
}


TEST_CASE( "Griewank", "[engine cmaes cpp]") {
    
    fitfun_initialize_simple("Griewank");
    
    auto engine = CmaesEngine( &fitfun, "./griewank/");
    engine.run();

    compare_best_with_ref( &engine, 
                           { -2.90498e-08, -1.56991e-09, 5.67367e-09, -8.38087e-09 },
                           9.93963,
                           EPSILON_TEST);
}


TEST_CASE( "Levy", "[engine cmaes cpp]") {
    
    fitfun_initialize_simple("Levy");
    
    auto engine = CmaesEngine( &fitfun, "./levy/");
    engine.run();

    compare_best_with_ref( &engine, 
                           { 1.0, 1.0, 1.0, 1.0 },
                           9.93963,
                           EPSILON_TEST);
}


TEST_CASE( "Multivariate Gaussian", "[engine cmaes cpp]") {
    
    fitfun_initialize_simple("multivariate_gaussian");
    
    auto engine = CmaesEngine( &fitfun, "./multivariate_G/");
    engine.run();

    compare_best_with_ref( &engine, 
                           { -2.55317e-06, 1.0, 2.0, 3.0 }, 
                           8.50505, 
                           EPSILON_TEST);
}


TEST_CASE( "Rastrigin", "[engine cmaes cpp]") {
    
    fitfun_initialize_simple("Rastrigin");
    
    auto engine = CmaesEngine( &fitfun, "./rastrigin/");
    engine.run();

    compare_best_with_ref( &engine,
                           { 1.74349e-09, 1.26991e-09, 2.85422e-09, -4.20577e-09 },
                           9.93963,
                           EPSILON_TEST);
}


TEST_CASE( "Rosenbrock", "[engine cmaes cpp]") {
    
    fitfun_initialize_simple("Rosenbrock");
    
    auto engine = CmaesEngine( &fitfun, "./rosenbrock/");
    engine.run();

    compare_best_with_ref( &engine, 
                           { 1.0, 1.0, 1.0, 1.0 },
                           9.93963,
                           EPSILON_TEST);
}


void compare_best_with_ref(CmaesEngine* engine, std::vector<double>
refBestEver, double refFunVal, double eps) {

    auto evo = engine->getEvo();

    Approx refBestFunValApprox = Approx(refFunVal).margin(eps);
   
    SECTION ("compare best function evaluation ever") {
        double bestFunVal = cmaes_Get(evo,"fbestever");
        REQUIRE(bestFunVal == refBestFunValApprox);
    }
    
    SECTION ("compare vector of best point") {
        for(int i = 0; i < evo->sp.N; ++i) {
            Approx refBestXApprox = Approx(refBestEver[i]).margin(eps); 
            REQUIRE( evo->rgxbestever[i] == refBestXApprox);
        }
    }
}
