
#define CATCH_CONFIG_MAIN

#include <catch.hpp>

#include "solver/linear/LinearSolverFactory.hpp"

#include <stdexcept>

#include <json/json.h>

using namespace neon;

/** Create a SPD matrix for solver testing */
SparseMatrix create_sparse_matrix()
{
    std::vector<Eigen::Triplet<double>> triplets = {{0, 0, 2.0},
                                                    {0, 1, -1.0},
                                                    {0, 2, 0.0},
                                                    {1, 0, -1.0},
                                                    {1, 1, 2.0},
                                                    {1, 2, -1.0},
                                                    {2, 0, 0.0},
                                                    {2, 1, -1.0},
                                                    {2, 2, 2.0}};

    SparseMatrix A(3, 3);
    A.setFromTriplets(std::begin(triplets), std::end(triplets));
    A.finalize();
    return A;
}

/** Create an example right hand side */
Vector create_right_hand_side()
{
    Eigen::VectorXd b(3);
    b(0) = 1.0;
    b(1) = 0.5;
    b(2) = 1.5;
    return b;
}

/** The known solution */
Vector solution()
{
    Vector x(3);
    x << 1.375, 1.75, 1.625;
    return x;
}

TEST_CASE("Linear solver test suite")
{
    SparseMatrix A = create_sparse_matrix();
    Vector b = create_right_hand_side();
    Vector x = b;

    SECTION("Preconditioned Conjugate Gradient Default")
    {
        std::string solver_str = "{\"Solver\":\"ConjugateGradient\"}";

        Json::Value solver_data;
        Json::Reader solver_file;
        REQUIRE(solver_file.parse(solver_str.c_str(), solver_data));

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0));
        REQUIRE((A * x - b).norm() == Approx(0.0));
    }
    SECTION("Preconditioned Conjugate Gradient Tolerance")
    {
        std::string solver_str = "{\"Solver\":\"ConjugateGradient\",\"Tolerance\":1e-8}";

        Json::Value solver_data;
        Json::Reader solver_file;
        REQUIRE(solver_file.parse(solver_str.c_str(), solver_data));

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0));
        REQUIRE((A * x - b).norm() == Approx(0.0));
    }
    SECTION("Preconditioned Conjugate Gradient Iterations")
    {
        std::string
            solver_str = "{\"Solver\":\"ConjugateGradient\",\"MaxIterations\":100}";

        Json::Value solver_data;
        Json::Reader solver_file;
        REQUIRE(solver_file.parse(solver_str.c_str(), solver_data));

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0));
        REQUIRE((A * x - b).norm() == Approx(0.0));
    }
    SECTION("Preconditioned Conjugate Gradient Iterations and Tolerance")
    {
        std::string solver_str = "{\"Solver\":\"ConjugateGradient\",\"MaxIterations\":"
                                 "100,\"Tolerance\":1e-8}";

        Json::Value solver_data;
        Json::Reader solver_file;
        REQUIRE(solver_file.parse(solver_str.c_str(), solver_data));

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0));
        REQUIRE((A * x - b).norm() == Approx(0.0));
    }
    SECTION("Preconditioned Bi-conjugate Gradient Stab Default")
    {
        std::string solver_str = "{\"Solver\":\"BiCGStab\"}";

        Json::Value solver_data;
        Json::Reader solver_file;
        REQUIRE(solver_file.parse(solver_str.c_str(), solver_data));

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0));
        REQUIRE((A * x - b).norm() == Approx(0.0));
    }
    SECTION("Preconditioned Bi-conjugate Gradient Stab Tolerance")
    {
        std::string solver_str = "{\"Solver\":\"BiCGStab\",\"Tolerance\":1e-8}";

        Json::Value solver_data;
        Json::Reader solver_file;
        REQUIRE(solver_file.parse(solver_str.c_str(), solver_data));

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0));
        REQUIRE((A * x - b).norm() == Approx(0.0));
    }
    SECTION("Preconditioned Bi-conjugate Gradient Stab Iterations")
    {
        std::string solver_str = "{\"Solver\":\"BiCGStab\",\"MaxIterations\":100}";

        Json::Value solver_data;
        Json::Reader solver_file;
        REQUIRE(solver_file.parse(solver_str.c_str(), solver_data));

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0));
        REQUIRE((A * x - b).norm() == Approx(0.0));
    }
    SECTION("Preconditioned Bi-conjugate Gradient Stab Iterations and Tolerance")
    {
        std::string solver_str = "{\"Solver\":\"BiCGStab\",\"MaxIterations\":100,"
                                 "\"Tolerance\":1e-8}";

        Json::Value solver_data;
        Json::Reader solver_file;
        REQUIRE(solver_file.parse(solver_str.c_str(), solver_data));

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0));
        REQUIRE((A * x - b).norm() == Approx(0.0));
    }
    SECTION("PaStiX")
    {
        std::string solver_str = "{\"Solver\":\"Pastix\"}";

        Json::Value solver_data;
        Json::Reader solver_file;
        REQUIRE(solver_file.parse(solver_str.c_str(), solver_data));

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0));
        REQUIRE((A * x - b).norm() == Approx(0.0));
    }
    SECTION("MUMPS")
    {
        std::string solver_str = "{\"Solver\":\"MUMPS\"}";

        Json::Value solver_data;
        Json::Reader solver_file;
        REQUIRE(solver_file.parse(solver_str.c_str(), solver_data));

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0));
        REQUIRE((A * x - b).norm() == Approx(0.0));
    }
    SECTION("SparseLU")
    {
        std::string solver_str = "{\"Solver\":\"SparseLU\"}";

        Json::Value solver_data;
        Json::Reader solver_file;
        REQUIRE(solver_file.parse(solver_str.c_str(), solver_data));

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0));
        REQUIRE((A * x - b).norm() == Approx(0.0));
    }
#ifdef ENABLE_CUDA
    SECTION("GPU Preconditioned Conjugate Gradient Default")
    {
        std::string solver_str = "{\"Solver\":\"ConjugateGradientGPU\"}";

        Json::Value solver_data;
        Json::Reader solver_file;
        REQUIRE(solver_file.parse(solver_str.c_str(), solver_data));

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0));
        REQUIRE((A * x - b).norm() == Approx(0.0));
    }
    SECTION("GPU Preconditioned Conjugate Gradient Tolerance")
    {
        std::string
            solver_str = "{\"Solver\":\"ConjugateGradientGPU\",\"Tolerance\":1e-8}";

        Json::Value solver_data;
        Json::Reader solver_file;
        REQUIRE(solver_file.parse(solver_str.c_str(), solver_data));

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0));
        REQUIRE((A * x - b).norm() == Approx(0.0));
    }
    SECTION("GPU Preconditioned Conjugate Gradient Iterations")
    {
        std::string
            solver_str = "{\"Solver\":\"ConjugateGradientGPU\",\"MaxIterations\":100}";

        Json::Value solver_data;
        Json::Reader solver_file;
        REQUIRE(solver_file.parse(solver_str.c_str(), solver_data));

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0));
        REQUIRE((A * x - b).norm() == Approx(0.0));
    }
    SECTION("GPU Preconditioned Conjugate Gradient Iterations and Tolerance")
    {
        std::string solver_str = "{\"Solver\":\"ConjugateGradientGPU\",\"MaxIterations\":"
                                 "100,\"Tolerance\":1e-8}";

        Json::Value solver_data;
        Json::Reader solver_file;
        REQUIRE(solver_file.parse(solver_str.c_str(), solver_data));

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0));
        REQUIRE((A * x - b).norm() == Approx(0.0));
    }
#endif
    SECTION("Error")
    {
        std::string solver_str = "{\"Solver\":\"PurpleMonkey\"}";

        Json::Value solver_data;
        Json::Reader solver_file;
        REQUIRE(solver_file.parse(solver_str.c_str(), solver_data));

        REQUIRE_THROWS_AS(make_linear_solver(solver_data), std::runtime_error);
    }
}
