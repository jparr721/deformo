#include "Rve.h"
#include <memory>

auto Rve::ToImplicitSurface() -> Tensor3r {
    const auto implicit_surface_generator =
        std::make_unique<ImplicitSurfaceGenerator>(width, height, depth);
    const auto inclusion = MakeVolumeAwareBinaryInclusion();

    return implicit_surface_generator->Generate(inclusion);
}

auto Rve::MakeVolumeAwareBinaryInclusion() -> BinaryInclusion {
    // TODO(@jparr721) - Actually claculate the volume fraction here
    return BinaryInclusion{15, 4, 5, 5, 5};
}

auto Rve::SetInclusionDimenions(const Vector3<unsigned int>& dimensions)
    -> void {
    inclusion_width = dimensions.x();
    inclusion_height = dimensions.y();
    inclusion_depth = dimensions.z();
}
