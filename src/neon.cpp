
#include "simulation_parser.hpp"
#include "exceptions.hpp"

#include <termcolor/termcolor.hpp>

int main(int argc, char* argv[])
{
    if (argc <= 1)
    {
        std::cerr << "No input file was provided.  Use <filename>.json\n";
        return 1;
    }
    try
    {
        neon::simulation_parser simulation(argv[1]);

        simulation.start();
    }
    catch (std::exception& error)
    {
        std::cout << std::endl
                  << std::string(2, ' ') << termcolor::red << termcolor::bold << error.what()
                  << termcolor::reset << std::flush << std::endl;
        return 1;
    }
    return 0;
}
