
#pragma once

#include <iostream>
#include <stdexcept>

namespace neon
{
/**
 * step_runtime_error should be thrown when something is wrong with the
 * computation during a time / load step that can potentially be corrected
 * with a change in load step.  These occur during the linear solution phase,
 * with element distortion and violations of internal variable updates as
 * examples
 */
struct computational_error : public std::runtime_error
{
    using std::runtime_error::runtime_error;
};

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

class DuplicateNameException : public Exception
{
public:
    DuplicateNameException(std::string duplParameter) : duplParameter(duplParameter) {}

    const char* what() const noexcept;

protected:
    std::string duplParameter;
};

class MaterialPropertyException : public Exception
{
public:
    MaterialPropertyException(std::string property) : property(property) {}

    const char* what() const noexcept;

protected:
    std::string property;
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
