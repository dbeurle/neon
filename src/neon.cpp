
#include "SimulationControl.hpp"

#include <iostream>

int main(int argc, char* argv[])
{
    using namespace neon;

    if (argc <= 1)
    {
        std::cerr << "No input file was provided.  Use <filename>.json\n";
        return 1;
    }

    SimulationControl simulation(argv[1]);

    simulation.start();

    return 0;
}
