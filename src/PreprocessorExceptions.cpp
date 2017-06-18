/*
 * neon - A finite element solver.
 *
 * For licensing please refer to the LICENSE.md file
 *
 * Copyright Darcy Beurle, 2016.
 */

#include "PreprocessorExceptions.hpp"

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

const char* NegativeElementAreaException::what() const noexcept
{
    std::cout << "\n!! Error: Element area is less than zero.  Please check mesh.\n";
    return nullptr;
}

DistortedElement::DistortedElement(int element, int quadraturePoint)
    : element(element), quadraturePoint(quadraturePoint)
{
}

const char* DistortedElement::what() const noexcept
{
    std::cout << "Element " + std::to_string(element) + ", quadrature point " +
                     std::to_string(quadraturePoint) + " is exhibiting excessive distortion\n";
    return nullptr;
}

// todo write list of valid element codes? For each analysis type? Are they even
// different?
const char* InvalidElementCodeException::what() const noexcept
{
    std::cout << "\n!! Error: Element code \"" << baseElementCode << "\" not supported for "
              << analysisType
              << " analysis.  Please file a bug if you think this should be supported.";
    return nullptr;
}

// todo write a list of valid element types?
const char* UnknownElementTypeException::what() const noexcept
{
    std::cout << "\n!! Error: Unable to determine element type\n";
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

const char* UnknownScenarioException::what() const noexcept
{
    std::cout << "\n!! Error: Scenario type \"" << scenario << "\" in " << input_file
              << " is not recognized.  The recognized scenario types are: \n"
                 "\t \"Heat\" \n"
                 "\t \"Elastostatics\" \n"
                 "\t \"PlaneStress\" \n"
                 "\t \"PlaneStrain\" \n"
                 "\t \"NaturalFrequency\" \n";
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

const char* JsonFileParseException::what() const noexcept
{
    std::cout << "\n!! Error: Please check input file: \n";
    std::cout << jsonError;
    return nullptr;
}

// todo SH do I need an if statement for "Part" (preprocessor) and "Scenario"
//(preprocessor) empty field?
const char* EmptyFieldException::what() const noexcept
{
    std::cout << "\n!! Error: " << emptyField << " not specified in " << input_file << ".  "
              << "Please specify using \"" << emptyField << "\" : \"<FieldDetails>\"";

    if (emptyField == "AnalysisName" || emptyField == "AnalysisType")
    {
        std::cout << "in \"Analysis\" section.";
    }
    else if (emptyField == "Material")
    {
        std::cout << "in \"Part\" section.";
    }
    else if (emptyField == "Solver")
    {
        std::cout << "in \"Analysis\" section using the format:";
        std::cout << "\"" << emptyField << "\" : [";
        std::cout << "\t {";
        std::cout << "\t\t \"Type\" : \"<SolverType>\"";
        std::cout << "\t }";
        std::cout << "\t ],";

        std::cout << "Accepted solver types are:";
        std::cout << "\t \"pCG\" - preconditioned (Jacobi) parallel symmetric Conjugate "
                     "Gradient";
        std::cout << "\t \"Pastix\" - parallel multifrontal direct LU solver with AMD "
                     "reordering";
        std::cout << "\t \"BiCGStab\" - parallel BiConjugate Gradient Stabilised";
    }
    else if (emptyField == "Output")
    {
        std::cout << "in \"Output\" section of simulation section of \"" << input_file
                  << "\" using the format:";
        std::cout << "\"" << emptyField << "\" : [";
        std::cout << "\t {";
        std::cout << "\t\t \"Fields\" : \"<OutputRequest1>\"";
        std::cout << "\t },";
        std::cout << "\t {";
        std::cout << "\t\t \"Fields\" : \"<OutputRequest2>\"";
        std::cout << "\t }";
        std::cout << "\t ]";

        std::cout << "Accepted output requests are: \n";
        std::cout << "\t \"Displacement\"";
        std::cout << "\t \"ModeShapes\"";
        std::cout << "\t \"VonMises\"";
        std::cout << "\t \"Stress\"";
        std::cout << "\t \"Strain\"";
        std::cout << "\t \"PrincipalStress\"";
        std::cout << "\t \"PrincipalStrain\"";
        std::cout << "\t \"Tresca\"";
        std::cout << "\t \"NodalTemperature\"\n";
        std::cout << "\t \"HeatFlux\"\n";
    }
    else if (emptyField == "BoundaryCondition")
    {
        std::cout << "in \"BoundaryCondition\" section of simulation section of " << input_file
                  << ".  The recognized boudary types are: \n"
                     "\t \"Traction\" \n"
                     "\t \"Pressure\" \n"
                     "\t \"DistributedForce\" \n"
                     "\t \"Robin\" \n"
                     "\t \"Temperature\" \n"
                     "\t \"HeatFlux\" \n";
        // todo Add the solver context to which each of these boundary conditions apply
    }
    else
    {
    }
    return nullptr;
}

const char* DuplicateNameException::what() const noexcept
{
    std::cout << "\n!! Error: Duplicate " << duplParameter << " names not allowed!  Please check "
              << duplParameter << " names in " << input_file;
    return nullptr;
}
