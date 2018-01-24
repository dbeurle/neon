
#define CATCH_CONFIG_MAIN

#include <catch.hpp>

#include "solver/linear/LinearSolverFactory.hpp"

#include <stdexcept>

#include "io/json.hpp"

using namespace neon;

constexpr auto ZERO_MARGIN = 1.0e-5;

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
vector create_right_hand_side()
{
    Eigen::VectorXd b(3);
    b(0) = 1.0;
    b(1) = 0.5;
    b(2) = 1.5;
    return b;
}

/** The known solution */
vector solution()
{
    vector x(3);
    x << 1.375, 1.75, 1.625;
    return x;
}

TEST_CASE("Linear solver test suite")
{
    SparseMatrix A = create_sparse_matrix();
    vector b = create_right_hand_side();
    vector x = b;

    SECTION("Preconditioned Conjugate Gradient Default")
    {
        json solver_data{{"Solver", "Iterative"}};

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Preconditioned Conjugate Gradient Tolerance")
    {
        json solver_data{{"Solver", "Iterative"}, {"Tolerance", 1.0e-8}};

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Preconditioned Conjugate Gradient Iterations")
    {
        json solver_data{{"Solver", "Iterative"}, {"MaxIterations", 100}};

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Preconditioned Conjugate Gradient Iterations and Tolerance")
    {
        json solver_data{{"Solver", "Iterative"}, {"MaxIterations", 100}, {"Tolerance", 1.0e-8}};

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Preconditioned Bi-conjugate Gradient Stab Default")
    {
        json solver_data{{"Solver", "Iterative"}};

        auto linear_solver = make_linear_solver(solver_data, false);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Preconditioned Bi-conjugate Gradient Stab Tolerance")
    {
        json solver_data{{"Solver", "Iterative"}, {"Tolerance", 1.0e-8}};

        auto linear_solver = make_linear_solver(solver_data, false);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Preconditioned Bi-conjugate Gradient Stab Iterations")
    {
        json solver_data{{"Solver", "Iterative"}, {"MaxIterations", 100}};

        auto linear_solver = make_linear_solver(solver_data, false);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Preconditioned Bi-conjugate Gradient Stab Iterations and Tolerance")
    {
        json solver_data{{"Solver", "Iterative"}, {"MaxIterations", 100}, {"Tolerance", 1.0e-8}};

        auto linear_solver = make_linear_solver(solver_data, false);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("PaStiX SPD")
    {
        json solver_data{{"Solver", "PaStiX"}};

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("PaStiX LU")
    {
        json solver_data{{"Solver", "PaStiX"}};

        auto linear_solver = make_linear_solver(solver_data, false);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("MUMPS SPD")
    {
        json solver_data{{"Solver", "MUMPS"}};

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("MUMPS LU")
    {
        json solver_data{{"Solver", "MUMPS"}};

        auto linear_solver = make_linear_solver(solver_data, false);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("SparseLU")
    {
        json solver_data{{"Solver", "Direct"}};

        auto linear_solver = make_linear_solver(solver_data, false);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("SparseLDLT")
    {
        json solver_data{{"Solver", "Direct"}};

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
#ifdef ENABLE_CUDA
    SECTION("GPU Preconditioned Conjugate Gradient Default")
    {
        json solver_data{{"Solver", "IterativeGPU"}};

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("GPU Preconditioned Conjugate Gradient Tolerance")
    {
        json solver_data{{"Solver", "Iterative"}, {"Tolerance", 1.0e-8}};

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("GPU Preconditioned Conjugate Gradient Iterations")
    {
        json solver_data{{"Solver", "Iterative"}, {"MaxIterations", 100}};

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("GPU Preconditioned Conjugate Gradient Iterations and Tolerance")
    {
        json solver_data{{"Solver", "Iterative"}, {"MaxIterations", 100}, {"Tolerance", 1.0e-8}};

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
#endif
    SECTION("Error")
    {
        json solver_data{{"Solver", "PurpleMonkey"}};

        REQUIRE_THROWS_AS(make_linear_solver(solver_data), std::runtime_error);
    }
}
