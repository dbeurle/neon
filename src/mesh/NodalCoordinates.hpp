
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

    /** Construct with a list of coordinates (x, y, z, x, y, z, ...) */
    explicit NodalCoordinates(vector coordinates);

    /** Construct with a list of coordinates in json format */
    explicit NodalCoordinates(Json::Value const& mesh_file);

    [[nodiscard]] auto size() const { return X.rows() / 3; }

    [[nodiscard]] vector const& coordinates() const { return X; }

    /** @return the coordinates using fancy indexing */
    [[nodiscard]] vector coordinates(List const& local_node_list) const;

protected:
    vector X; //!< Reference configuration encoded as (x1, y1, z1, x2, y2, z2)
};
}
