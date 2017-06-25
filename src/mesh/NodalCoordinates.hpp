
#pragma once

#include <json/forwards.h>

#include "numeric/DenseTypes.hpp"

namespace neon
{
/** NodalCoordinates stores the list of coordinates of a discretized geometry */
class NodalCoordinates
{
public:
    NodalCoordinates() = default;

    /** Construct with a list of coordinates (x, y, z, x, y, z, ...) */
    explicit NodalCoordinates(Vector coordinates);

    /** Construct with a list of coordinates in json format */
    explicit NodalCoordinates(Json::Value const& mesh_file);

    auto size() const { return X.rows() / 3; }

    Vector const& coordinates() const { return X; }

    /** @return the coordinates using fancy indexing */
    Vector operator[](List const& local_node_list) const;

protected:
    Vector X; //!< Reference configuration encoded as (x1, y1, z1, x2, y2, z2)
};
}
