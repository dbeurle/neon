
#pragma once

#include <iostream>
#include <stdexcept>

namespace neon
{
/// computational_error should be thrown when something is wrong with the
/// computation that can potentially be corrected with a change in load step.
/// These occur during the linear solution phase, with element distortion and
/// violations of internal variable updates as examples
struct computational_error : public std::domain_error
{
    using std::domain_error::domain_error;
};

class Exception : public std::exception
{
public:
    Exception() = default;

    void setInputFileName(std::string const& fn) { input_file = fn; }

protected:
    std::string input_file;
};

template <typename KeyTp>
class KeyNotFoundInMap : public Exception
{
public:
    KeyNotFoundInMap(KeyTp missing_key) : missing_key(missing_key) {}

    char const* what() const noexcept
    {
        std::cout << "\n!! Error: Key " << missing_key << " not found in an internal datastructure\n";
        return nullptr;
    }

protected:
    KeyTp missing_key;
};
}
