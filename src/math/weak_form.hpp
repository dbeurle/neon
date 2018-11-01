
#pragma once

namespace neon::fem
{
/// Allowable options of bilinear form
enum class form {
    /// standard linear form (w, h)
    linear,
    /// gradient linear form (gradient(w), gradient(h))
    linear_gradient,
    /// standard bilinear form (w, u)
    bilinear,
    /// gradient form (gradient(w), gradient(u))
    bilinear_gradient,
};
}
