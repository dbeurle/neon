
#pragma once

#include <exception>
#include <iostream>
#include <string>

namespace neon
{
/*
 * Exception header file for handling all exceptions thrown
 */
class Exception : public std::exception
{
public:
    Exception() = default;

    void setInputFileName(std::string const& fn) { input_file = fn; }

protected:
    std::string input_file;
};

class NoInputException : public Exception
{
public:
    NoInputException() = default;

    const char* what() const noexcept;
};

class InvalidExtensionException : public Exception
{
public:
    InvalidExtensionException(std::string extension) : extension(extension) {}

    const char* what() const noexcept;

protected:
    std::string extension;
};

class JsonFileParseException : public Exception
{
public:
    JsonFileParseException(std::string jsonError) : jsonError(jsonError) {}

    const char* what() const noexcept;

protected:
    std::string jsonError;
};

class EmptyFieldException : public Exception
{
public:
    EmptyFieldException(std::string emptyField) : emptyField(emptyField) {}

    const char* what() const noexcept;

protected:
    std::string emptyField;
};

class DuplicateNameException : public Exception
{
public:
    DuplicateNameException(std::string duplParameter) : duplParameter(duplParameter) {}

    const char* what() const noexcept;

protected:
    std::string duplParameter;
};

class UnknownOutputRequestException : public Exception
{
public:
    UnknownOutputRequestException(std::string outputRequest) : outputRequest(outputRequest) {}

    const char* what() const noexcept;

protected:
    std::string outputRequest;
};

class UnknownScenarioException : public Exception
{
public:
    UnknownScenarioException(std::string scenario) : scenario(scenario) {}

    const char* what() const noexcept;

protected:
    std::string scenario;
};

class UnknownBoundaryTypeException : public Exception
{
public:
    UnknownBoundaryTypeException(const std::string& boundaryType) : boundaryType(boundaryType) {}

    const char* what() const noexcept;

protected:
    std::string boundaryType;
};

class MaterialPropertyException : public Exception
{
public:
    MaterialPropertyException(std::string property) : property(property) {}

    const char* what() const noexcept;

protected:
    std::string property;
};

class NegativeElementAreaException : public Exception
{
public:
    NegativeElementAreaException() {}

    const char* what() const noexcept;
};

class InvalidElementCodeException : public Exception
{
public:
    InvalidElementCodeException(short baseElementCode, std::string analysisType)
        : baseElementCode(baseElementCode), analysisType(analysisType)
    {
    }

    const char* what() const noexcept;

private:
    short baseElementCode;
    std::string analysisType;
};

class DistortedElement : public Exception
{
public:
    DistortedElement(int element, int quadrature_point);

    const char* what() const noexcept;

private:
    int element, quadrature_point;
};

class UnknownElementTypeException : public Exception
{
public:
    UnknownElementTypeException() = default;

    const char* what() const noexcept;
};

class PartNameException : public Exception
{
public:
    PartNameException(std::string partName) : partName(partName) {}

    const char* what() const noexcept;

protected:
    std::string partName;
};

template <typename KeyTp>
class KeyNotFoundInMap : public Exception
{
public:
    KeyNotFoundInMap(KeyTp missing_key) : missing_key(missing_key) {}

    const char* what() const noexcept
    {
        std::cout << "\n!! Error: Key " << missing_key << " not found in an internal datastructure\n";
        return nullptr;
    }

protected:
    KeyTp missing_key;
};
}
