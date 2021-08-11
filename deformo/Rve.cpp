#include "Rve.h"
#include <memory>

auto Rve::ToImplicitSurface() -> Tensor3r {
    //if (homogenous) {
        Tensor3r output(width, height, depth);
        output.SetConstant(material_1.number);
        return output;
    //}

    //const auto implicit_surface_generator =
    //    std::make_unique<ImplicitSurfaceGenerator>(width, height, depth);
    //const auto inclusion = MakeVolumeAwareBinaryInclusion();

    //return implicit_surface_generator->Generate();
}
