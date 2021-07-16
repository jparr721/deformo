#include "Rve.h"

auto Material::Empty() -> bool { return E == -1. && v == -1. && name.empty(); }

auto Rve::AssignSections(const std::string& material_name) -> void {}

auto Rve::SetVolumeFractionForMaterial(const std::string& material_name,
                                       Real volume_fraction) -> void {
    
}
