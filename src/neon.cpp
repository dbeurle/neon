
#include "simulation_parser.hpp"

#include "exceptions.hpp"

#include <stdexcept>
#include <termcolor/termcolor.hpp>

int main(int argc, char* argv[])
{
    using namespace neon;

    if (argc <= 1)
    {
        std::cerr << "No input file was provided.  Use <filename>.json\n";
        return 1;
    }

    try
    {
        simulation_parser simulation(argv[1]);

        simulation.start();
    }
    catch (std::runtime_error& error)
    {
        std::cout << std::endl
                  << std::string(2, ' ') << termcolor::red << termcolor::bold << error.what()
                  << termcolor::reset << std::flush << std::endl;
        return 1;
    }
    catch (std::domain_error& error)
    {
        std::cout << std::endl
                  << std::string(2, ' ') << termcolor::red << termcolor::bold << error.what()
                  << termcolor::reset << std::flush << std::endl;
        return 1;
    }
    return 0;
}
