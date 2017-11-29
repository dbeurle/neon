
#pragma once

#include <queue>
#include <tuple>
#include <vector>

#include <json/forwards.h>

namespace neon
{
/**
 * AdaptiveLoadStep responsibility is to handle the pseudo time step for each
 * SimulationCase (see input file).
 *
 * The load factor is always between zero and one for a given load case, which
 * is defined as time_n and time_n+1.  However, to retain convergence during
 * the non-linear step, cutbacks in the load factor can be required.
 *
 * These cutbacks are determined from the initial, minimum and maximum time
 * increments specified by the user.  The convergence state of the non-linear
 * solution scheme is given at the end and this is used to update the internal
 * state of the adaptive load step algorithm to decide when to reduce the load
 * factor.
 *
 * Consider the simpliest case of a ramped load:
 *
 *                        applied load
 *                             ^
 *                             |
 *                         1 --|        *
 *                             |      *
 *                             |    *
 *                             |  *
 *                         0 --|*______________> pseudo time
 *                             0        1
 *
 * Where the load factor, α, is between 0 and 1, pseudo time is between 0 and 1
 * and the applied load is between 0 and 1 (\sa Dirichlet).
 * Cutbacks in this case are handled by reducing the load factor.  The current
 * factor is given by the load factor.  The decision on the ramp or instantaneous
 * is handled by the relevant boundary condition.
 *
 * Consider a load-relaxation-unload case, where each boundary condition
 * is applied, held and unapplied regardless of internal load incrementation.
 *
 *         applied load
 *              ^
 *              |
 *          1 --|        * * * * * * * * * *
 *              |      *                     *
 *              |    *                         *
 *              |  *                             *
 *          0 --|*_________________________________*_______> pseudo time
 *              0       1                 3        4
 *
 * In this case, the load factor gives the progress through the time step.
 * The boundary conditions are responsible for determining their application
 * method (ramp or instantaneous), and are given a load factor for interpolation
 * purposes.
 *
 * The unloading step is handled by the boundary conditions as they will be
 * responsible for interpolation of their current and final states and will
 * require the load factor from this algorithm.  If g is the applied value
 * for a Dirichlet condition, then g depends on α (load factor) such that g(α).
 *
 * TODO This class can be extended by including a smart prediction algorithm
 * for the determination of the best next step.
 *
 */
class AdaptiveLoadStep
{
public:
    AdaptiveLoadStep(Json::Value const& increment_data, std::vector<double> mandatory_time_history);

    /** Check if the load increment is finalised */
    [[nodiscard]] bool is_fully_applied() const { return is_applied; }

    /** Get the global time (including past load cases) */
    [[nodiscard]] double time() const { return total_time + last_converged_time; }

    /** Get the time only for the current load case */
    [[nodiscard]] double step_time() const { return current_time; }

    /** Get the time only for the last converged load case */
    [[nodiscard]] double last_step_time() const { return last_converged_time; }

    /** Get the pseudo time step size */
    [[nodiscard]] double increment() const { return current_time - last_converged_time; }

    /** The number of steps taken for all time */
    [[nodiscard]] auto step() const { return successful_increments; }

    /** Update the convergence state to determine the next increment */
    void update_convergence_state(bool const is_converged);

    void reset(Json::Value const& new_increment_data);

protected:
    void parse_input(Json::Value const& increment_data, double const maximum_mandatory_time);

    void check_increment_data(Json::Value const& increment_data);

    [[nodiscard]] bool is_highly_nonlinear() const
    {
        return consecutive_unconverged > 0 || consecutive_converged < 4;
    }

protected:
    int const increment_limit = 10; //!< Maximum allowable increments
    int successful_increments = 0;  //!< Number of converged steps

    int consecutive_converged = 0;   //!< Number of consecutive successful attempts
    int consecutive_unconverged = 0; //!< Number of consecutive unsuccessful attempts

    double initial_time = 1.0;
    double final_time = 1.0;
    double current_time = 1.0;

    double total_time = 0.0; //!< Time history for multi-step simulations

    double last_converged_time = 0.0;           //!< Last time when convergence was reached
    double last_converged_time_step_size = 0.0; //!< Last good timestep size

    double minimum_increment; //!< Minimum increment allowed by the algorithm
    double maximum_increment; //!< Maximum increment allowed by the algorithm

    bool is_applied = false;

    std::priority_queue<double, std::vector<double>, std::greater<double>> time_queue;
};
}
