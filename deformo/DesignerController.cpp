#include "DesignerController.h"

void DesignerController::SetImplicitSurfaceHeight(int value) {
    implicit_surface_height_ = value;
}

void DesignerController::SetImplicitSurfaceWidth(int value) {
    implicit_surface_width_ = value;
}

void DesignerController::SetImplicitSurfaceDepth(int value) {
    implicit_surface_depth_ = value;
}

void DesignerController::SetNumberOfInclusions(int value) {
    material_2_number_of_inclusions_ = value;
}

void DesignerController::SetMaterialOneName(const std::string& value) {
    material_1.name = value;
}

void DesignerController::SetMaterialOnePoissonsRatio(double value) {
    material_1.v = value;
}

void DesignerController::SetMaterialOneYoungsModulus(double value) {
    material_1.E = value;
}

void DesignerController::SetMaterialTwoName(const std::string& value) {
    material_2.name = value;
}

void DesignerController::SetMaterialTwoPoissonsRatio(double value) {
    material_2.v = value;
}

void DesignerController::SetMaterialTwoYoungsModulus(double value) {
    material_2.E = value;
}

void DesignerController::SetSquareShapedMaterial(const bool checked) {
    is_square_ = checked;
}

void DesignerController::SetUniformMaterial(const bool checked,
                                            const Ui::deformoClass& ui) {
    is_uniform_ = checked;
    if (is_uniform_) {
        ui.designer_material_2_name_line_edit->setDisabled(true);
        ui.designer_material_2_poissons_ratio_double_spinbox->setDisabled(true);
        ui.designer_material_2_youngs_modulus_double_spinbox->setDisabled(true);
        ui.designer_isotropic_material_checkbox->setDisabled(true);
    } else {
        ui.designer_material_2_name_line_edit->setDisabled(false);
        ui.designer_material_2_poissons_ratio_double_spinbox->setDisabled(
            false);
        ui.designer_material_2_youngs_modulus_double_spinbox->setDisabled(
            false);
        ui.designer_isotropic_material_checkbox->setDisabled(false);
    }
}

void DesignerController::SetIsotropicMaterial(const bool checked) {
    is_isotropic_ = checked;
}

void DesignerController::ComputeDesignedShapeButtonPressed() {
    material_1 =
        MaterialFromEandv(1, material_1.name, material_1.E, material_1.v);

    if (!is_uniform_) {
        material_2 =
            MaterialFromEandv(1, material_2.name, material_2.E, material_2.v);
    }
}

void DesignerController::SetInclusionHeight(int value) {
    inclusion_height_ = value;
}

void DesignerController::SetInclusionWidth(int value) {
    inclusion_width_ = value;
}

void DesignerController::SetInclusionDepth(int value) {
    inclusion_depth_ = value;
}

void DesignerController::SetSquareShapedInclusion(bool checked) {
    inclusion_is_square_ = checked;
}
