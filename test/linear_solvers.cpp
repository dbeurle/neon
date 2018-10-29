
#include <catch2/catch.hpp>

#include "solver/linear/linear_solver_factory.hpp"

#include <stdexcept>

#include "io/json.hpp"

using namespace neon;

constexpr auto ZERO_MARGIN = 1.0e-5;

/** Create a SPD matrix for solver testing */
sparse_matrix create_sparse_matrix()
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

    sparse_matrix A(3, 3);
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
    sparse_matrix A = create_sparse_matrix();
    vector b = create_right_hand_side();
    vector x = b;

    SECTION("Preconditioned Conjugate Gradient Default")
    {
        json solver_data{{"type", "iterative"}};

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Preconditioned Conjugate Gradient Tolerance")
    {
        json solver_data{{"type", "iterative"}, {"tolerance", 1.0e-8}};

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Preconditioned Conjugate Gradient Iterations")
    {
        json solver_data{{"type", "iterative"}, {"maximum_iterations", 100}};

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Preconditioned Conjugate Gradient Iterations and Tolerance")
    {
        json solver_data{{"type", "iterative"}, {"maximum_iterations", 100}, {"tolerance", 1.0e-8}};

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Preconditioned Bi-conjugate Gradient Stab Default")
    {
        json solver_data{{"type", "iterative"}};

        auto linear_solver = make_linear_solver(solver_data, false);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Preconditioned Bi-conjugate Gradient Stab Tolerance")
    {
        json solver_data{{"type", "iterative"}, {"tolerance", 1.0e-8}};

        auto linear_solver = make_linear_solver(solver_data, false);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Preconditioned Bi-conjugate Gradient Stab Iterations")
    {
        json solver_data{{"type", "iterative"}, {"maximum_iterations", 100}};

        auto linear_solver = make_linear_solver(solver_data, false);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Preconditioned Bi-conjugate Gradient Stab Iterations and Tolerance")
    {
        json solver_data{{"type", "iterative"}, {"maximum_iterations", 100}, {"tolerance", 1.0e-8}};

        auto linear_solver = make_linear_solver(solver_data, false);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("PaStiX SPD")
    {
        json solver_data{{"type", "PaStiX"}};

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("PaStiX LU")
    {
        json solver_data{{"type", "PaStiX"}};

        auto linear_solver = make_linear_solver(solver_data, false);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("MUMPS SPD")
    {
        json solver_data{{"type", "MUMPS"}};

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("MUMPS LU")
    {
        json solver_data{{"type", "MUMPS"}};

        auto linear_solver = make_linear_solver(solver_data, false);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("SparseLU")
    {
        json solver_data{{"type", "direct"}};

        auto linear_solver = make_linear_solver(solver_data, false);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("SparseLDLT")
    {
        json solver_data{{"type", "direct"}};

        auto linear_solver = make_linear_solver(solver_data);

        linear_solver->solve(A, x, b);

        REQUIRE((x - solution()).norm() == Approx(0.0).margin(ZERO_MARGIN));
        REQUIRE((A * x - b).norm() == Approx(0.0).margin(ZERO_MARGIN));
    }
    SECTION("Error")
    {
        json solver_data{{"type", "PurpleMonkey"}};

        REQUIRE_THROWS_AS(make_linear_solver(solver_data), std::domain_error);
    }
}
