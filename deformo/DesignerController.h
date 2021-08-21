#pragma once

#include "CsvFile.h"
#include "HomogenizationDatasetGeneratorDialog.h"
#include "Material.h"
#include "Mesh.h"
#include "Numerics.h"
#include "ui_deformo.h"
#include <QObject>
#include <string>

class DesignerController : public QObject {
    Q_OBJECT
  public:
    Material material_1;
    Material material_2;

    DesignerController(const std::shared_ptr<Mesh> mesh,
                       const Ui::deformoClass& ui);
    ~DesignerController();

    // Implicit surface
    void SetImplicitSurfaceHeight(int value);
    void SetImplicitSurfaceWidth(int value);
    void SetImplicitSurfaceDepth(int value);

    // Designer -- Size of inclusion
    void SetNumberOfInclusions(int value);

    // Designer -- Material 1 Specifications
    void SetMaterialOneName(const std::string& value);
    void SetMaterialOnePoissonsRatio(double value);
    void SetMaterialOneYoungsModulus(double value);

    // Designer -- Material 2 Specifications
    void SetMaterialTwoName(const std::string& value);
    void SetMaterialTwoPoissonsRatio(double value);
    void SetMaterialTwoYoungsModulus(double value);

    void SetSquareShapedMaterial(bool checked);
    void SetUniformMaterial(bool checked, const Ui::deformoClass& ui);
    void SetIsotropicMaterial(bool checked);

    void ComputeDesignedShapeButtonPressed();

    void SetInclusionHeight(int value);
    void SetInclusionWidth(int value);
    void SetInclusionDepth(int value);

    void SetSquareShapedInclusion(bool checked);

    // Dataset Generator
    void SetCSVPathButtonClicked(const Ui::deformoClass& ui);
    void SetOutputCSVFileName(const std::string& value);
    void SetGenerator(const std::string& value);
    void ComputeDatasetButtonPressed();
    void SetupAndDisplayHomogenizationDialog();
    void SetNumberOfEntries(int value);

  public slots:
    void HomogenizationDialogSubmitted();

  private:
    enum class Generator {
        kHomogenization = 0,
    };

    // Implicit Surface Options
    int implicit_surface_height_ = 0;
    int implicit_surface_width_ = 0;
    int implicit_surface_depth_ = 0;

    // Inclusion Options
    int inclusion_height_ = 0;
    int inclusion_width_ = 0;
    int inclusion_depth_ = 0;

    // Dataset Generator
    int number_of_dataset_entries_ = 1;

    // Square shaped object
    bool is_square_ = false;

    // Is material made of only material 1
    bool is_uniform_ = false;

    // Isotropic Material Generator
    bool is_isotropic_ = false;

    // Square shaped inclusion
    bool inclusion_is_square_ = false;

    // Inclusion ratio for material 2
    int material_2_number_of_inclusions_ = 0;

    std::string output_csv_path_;
    std::string output_csv_filename_;
    Generator generator_ = Generator::kHomogenization;

    std::shared_ptr<Mesh> mesh_;
    std::unique_ptr<Homogenization> homogenization_;

    HomogenizationDatasetGeneratorDialog* homogenization_dialog_;

    Ui::deformoClass ui_;
};
