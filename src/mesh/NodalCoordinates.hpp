
#pragma once

#include "numeric/DenseMatrix.hpp"
#include "numeric/IndexTypes.hpp"

#include <json/forwards.h>

namespace neon
{
/** NodalCoordinates stores the list of coordinates of a discretized geometry */
class NodalCoordinates
{
public:
    NodalCoordinates() = default;

    /** Construct with a list of coordinates */
    explicit NodalCoordinates(matrix3x const coordinates);

    /** Construct with a list of coordinates in json format */
    explicit NodalCoordinates(Json::Value const& mesh_file);

    [[nodiscard]] auto size() const { return X.cols(); }

    [[nodiscard]] matrix3x const& coordinates() const;

    /** @return the coordinates using fancy indexing */
    [[nodiscard]] matrix3x coordinates(std::vector<int64> const& local_node_list) const;

protected:
    matrix3x X; //!< Reference configuration encoded as (x1, y1, z1, x2, y2, z2)
};
}
