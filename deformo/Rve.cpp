#include "Rve.h"

Rve::Rve(const Vector3<int>& size, const Material& material)
    : material_1_(material) {
    height_ = size.x();
    width_ = size.y();
    depth_ = size.z();
}

Rve::Rve(const Vector3<int>& size, const Material& material_1,
         const Material& material_2)
    : material_1_(material_1), material_2_(material_2) {
    height_ = size.x();
    width_ = size.y();
    depth_ = size.z();
}

auto Rve::Homogenize() -> void {
    deformo_assertion_.Assert(contains_surface_mesh_, __FUNCTION__, __FILE__,
                              __LINE__, "No surface mesh found");
    if (material_2_.IsInit()) {
		homogenization_ = std::make_unique<Homogenization>(
			surface_mesh_, material_1_, material_2_);
    } else {
		homogenization_ = std::make_unique<Homogenization>(
			surface_mesh_, material_1_);
    }

    homogenization_->Solve();
    C_ = homogenization_->Stiffness();
    is_homogenized_ = true;
}

auto Rve::ComputeRenderableMesh(MatrixXr& V, MatrixX<int>& F) -> void {
    deformo_assertion_.Assert(contains_surface_mesh_, __FUNCTION__, __FILE__,
                              __LINE__, "No surface mesh found");
    // Need to add the padding layers so that way marching cubes works properly
    surface_mesh_ = generator_->AddSquarePaddingLayers();
    marching_cubes_ = std::make_unique<MarchingCubes>(
        material_1_.number, kCellLength, surface_mesh_.Instance().data());
    marching_cubes_->GenerateGeometry(V, F, height_ + 1, width_ + 1,
                                      depth_ + 1);
}

auto Rve::ComputeSurfaceMesh() -> void {
    generator_ = std::make_unique<ImplicitSurfaceGenerator<Real>>(
        height_, width_, depth_, material_1_.number);

    surface_mesh_ = generator_->Generate();

    contains_surface_mesh_ = true;
}

auto Rve::ComputeSurfaceMesh(const Vector3<int>& inclusion_size,
                             const int n_inclusions, const bool is_isotropic)
    -> void {
    const ImplicitSurfaceGenerator<Real>::Inclusion inclusion{
        n_inclusions,       inclusion_size.x(), inclusion_size.x(),
        inclusion_size.y(), inclusion_size.z(),
    };

    const ImplicitSurfaceGenerator<Real>::ImplicitSurfaceMicrostructure
        microstructure = ImplicitSurfaceGenerator<
            Real>::ImplicitSurfaceMicrostructure::kComposite;

    const ImplicitSurfaceGenerator<Real>::ImplicitSurfaceCharacteristics
        characteristics =
            is_isotropic
                ? ImplicitSurfaceGenerator<
                      Real>::ImplicitSurfaceCharacteristics::kIsotropic
                : ImplicitSurfaceGenerator<
                      Real>::ImplicitSurfaceCharacteristics::kAnisotropic;

    const int material_2_number = material_2_.IsInit() ? material_2_.number : 0;

    generator_ = std::make_unique<ImplicitSurfaceGenerator<Real>>(
        height_, width_, depth_, characteristics, microstructure, inclusion,
        material_1_.number, material_2_number);

    surface_mesh_ = generator_->Generate();

    contains_surface_mesh_ = true;
}
