
#pragma once

namespace neon::fem
{
/// Classification of bilinear form where
/// standard -> (w, u)
/// gradient -> (gradient(w), gradient(u))
enum class bilinear { standard, gradient };

/// Classification of linear form where
/// standard -> (w, h)
/// gradient -> (gradient(w), gradient(h))
enum class linear { standard, gradient };
}
