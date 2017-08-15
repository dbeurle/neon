
#include "Exceptions.hpp"

namespace neon
{
const char* UnknownBoundaryTypeException::what() const noexcept
{
    std::cout << "\n!! Error: Boundary type \"" << boundaryType << "\" in " << input_file
              << " not recognized.  The recognized boudary types are: \n"
                 "Traction \n"
                 "Pressure \n"
                 "DistributedForce \n"
                 "Robin \n"
                 "Temperature \n"
                 "HeatFlux\n";
    // todo Add the solver context to which each of these boundary conditions apply
    return nullptr;
}

const char* MaterialPropertyException::what() const noexcept
{
    std::cout << "\n!! Error: " << property << " has not been provided in " << input_file;
    return nullptr;
}

const char* PartNameException::what() const noexcept
{
    std::cout << "\n!! Error: Part name \"" << partName << "\" not found in gmsh file.";
    std::cout << "Please check that part name in " << input_file
              << ".neon file matches highest dimension physical name in .msh file.\n";
    return nullptr;
}

const char* UnknownOutputRequestException::what() const noexcept
{
    std::cout << "\n!! Error: Output request " << outputRequest << " in " << input_file
              << " is not recognized.  The recognized output requests are: \n"
                 "\t \"Displacement\" \n"
                 "\t \"ModeShapes\" \n"
                 "\t \"VonMises\" \n"
                 "\t \"CauchyStress\" \n"
                 "\t \"CauchyStrain\" \n"
                 "\t \"PrincipalStress\" \n"
                 "\t \"PrincipalStrain\" \n"
                 "\t \"Tresca\" \n"
                 "\t \"NodalTemperature\" \n"
                 "\t \"HeatFlux\" \n";
    return nullptr;
}

const char* NoInputException::what() const noexcept
{
    std::cout << "\n!! Error: No input file found.  An input file needs to be provided: "
              << "\"<filename>.neon\"\n";
    return nullptr;
}

const char* InvalidExtensionException::what() const noexcept
{
    std::cout << "\n!! Error: Extension \"" << extension << "\" is not supported.\n"
              << "Supported extension is \".neon\"";
    return nullptr;
}

const char* DuplicateNameException::what() const noexcept
{
    std::cout << "\n!! Error: Duplicate " << duplParameter
              << " names not allowed!  Please check " << duplParameter << " names in "
              << input_file;
    return nullptr;
}
}
