#include "DesignerController.h"
#include "MarchingCubes.h"
#include "QTUtils.h"
#include "Rve.h"
#include <QFileDialog>
#include <QPushButton>
#include <QString>
#include <random>
#include <igl/writeOBJ.h>

DesignerController::DesignerController(const std::shared_ptr<Mesh> mesh,
                                       const Ui::deformoClass& ui)
    : mesh_(mesh), ui_(ui) {
    homogenization_dialog_ =
        new HomogenizationDatasetGeneratorDialog(ui_.centralWidget);

    ui_.designer_dataset_generator_progressbar->setVisible(false);
}

DesignerController::~DesignerController() { delete homogenization_dialog_; }

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
        ui.designer_isotropic_material_checkbox->setChecked(true);
    } else {
        ui.designer_material_2_name_line_edit->setDisabled(false);
        ui.designer_material_2_poissons_ratio_double_spinbox->setDisabled(
            false);
        ui.designer_material_2_youngs_modulus_double_spinbox->setDisabled(
            false);
        ui.designer_isotropic_material_checkbox->setDisabled(false);
        ui.designer_isotropic_material_checkbox->setChecked(false);
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
            MaterialFromEandv(0, material_2.name, material_2.E, material_2.v);
    }

    const auto rve = std::make_unique<Rve>(
        Vector3<int>(implicit_surface_height_, implicit_surface_width_,
                     implicit_surface_depth_),
        material_1, material_2);

    if (is_uniform_) {
        rve->ComputeSurfaceMesh();
    } else {
        rve->ComputeSurfaceMesh(
            Vector3<int>(inclusion_height_, inclusion_width_, inclusion_depth_),
            material_2_number_of_inclusions_, is_isotropic_);
    }

    MatrixXr V;
    MatrixX<int> F;
    rve->ComputeRenderableMesh(V, F);
    igl::writeOBJ("cuboid.obj", V, F);
    mesh_->RefreshData(V, F);
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

void DesignerController::SetCSVPathButtonClicked(const Ui::deformoClass& ui) {
    const QString q_dirname =
        QFileDialog::getExistingDirectory(nullptr, "Choose CSV File");
    ui.designer_dataset_generator_csv_path_line_edit->setText(q_dirname);
    const std::string dirname = utils::qt::QStringToString(q_dirname);
    output_csv_path_ = dirname + "/";
}

void DesignerController::SetOutputCSVFileName(const std::string& value) {
    output_csv_filename_ = value;
}

void DesignerController::SetGenerator(const std::string& value) {
    if (value == "Homogenization") {
        generator_ = Generator::kHomogenization;
    }
}

void DesignerController::ComputeDatasetButtonPressed() {
    switch (generator_) {
    case Generator::kHomogenization:
        SetupAndDisplayHomogenizationDialog();
    }
}

void DesignerController::HomogenizationDialogSubmitted() {
    ui_.designer_dataset_generator_progressbar->setVisible(true);
    ui_.designer_dataset_generator_progressbar->setValue(0);
    homogenization_dialog_->accept();

    assert(!output_csv_path_.empty());
    assert(!output_csv_filename_.empty());
    const std::string filename = output_csv_path_ + output_csv_filename_;

    CsvFile<std::string> csv_file(
        filename,
        std::vector<std::string>({"Voxel", "Rows", "Cols", "Layers", "Coefficients"}));

    for (int i = 0; i < number_of_dataset_entries_; ++i) {
        ui_.designer_dataset_generator_progressbar->setValue(
            ((i + 1) / number_of_dataset_entries_) * 100);

        const Material m =
            MaterialFromEandv(1, "material", homogenization_dialog_->E_,
                              homogenization_dialog_->v_);
        const auto rve = std::make_unique<Rve>(
            Vector3<int>{homogenization_dialog_->cube_dimensions_,
                         homogenization_dialog_->cube_dimensions_,
                         homogenization_dialog_->cube_dimensions_},
            m);

        // Select random inclusion dimension in range
        std::default_random_engine generator;
        const std::uniform_int_distribution<int> inclusion_dims_distribution(
            homogenization_dialog_->min_inclusion_dimensions_,
            homogenization_dialog_->max_inclusion_dimensions_);

        const int dim = inclusion_dims_distribution(generator);

        const std::uniform_int_distribution<int> inclusions_distribution(
            homogenization_dialog_->min_inclusions_,
            homogenization_dialog_->max_inclusions_);

        const int inclusions = inclusions_distribution(generator);

        const Vector3<int> inclusion_size(dim, dim, dim);
        rve->ComputeSurfaceMesh(inclusion_size, inclusions, false);
        rve->Homogenize();

        // Pull Params from homogenized system
        // Voxel
        const Tensor3r _voxel = rve->Homogenized()->Voxel();
        const std::string voxel = rve->Homogenized()->Voxel().ToString();
        const std::string rows = std::to_string(_voxel.Dimension(0));
        const std::string cols = std::to_string(_voxel.Dimension(1));
        const std::string layers = std::to_string(_voxel.Dimension(2));
        const VectorXr _coefficients = rve->Homogenized()->CoefficientVector();
        std::stringstream ss;
        ss << _coefficients.transpose();
        const std::string coefficients = ss.str();
        const std::vector<std::string> entries(
            {voxel, rows, cols, layers, coefficients});

        csv_file << entries;
    }
}

void DesignerController::SetupAndDisplayHomogenizationDialog() {
    connect(homogenization_dialog_->submit_button, &QPushButton::pressed, this,
            &DesignerController::HomogenizationDialogSubmitted);

    connect(homogenization_dialog_->min_inclusion_dimensions,
            QOverload<int>::of(&QSpinBox::valueChanged), homogenization_dialog_,
            &HomogenizationDatasetGeneratorDialog::SetMinInclusionDimensions);
    connect(homogenization_dialog_,
            &HomogenizationDatasetGeneratorDialog::OnSetMinInclusionDimensions,
            homogenization_dialog_->min_inclusion_dimensions,
            &QSpinBox::setValue);

    connect(homogenization_dialog_->max_inclusion_dimensions,
            QOverload<int>::of(&QSpinBox::valueChanged), homogenization_dialog_,
            &HomogenizationDatasetGeneratorDialog::SetMaxInclusionDimensions);
    connect(homogenization_dialog_,
            &HomogenizationDatasetGeneratorDialog::OnSetMaxInclusionDimensions,
            homogenization_dialog_->max_inclusion_dimensions,
            &QSpinBox::setValue);

    connect(homogenization_dialog_->min_inclusions,
            QOverload<int>::of(&QSpinBox::valueChanged), homogenization_dialog_,
            &HomogenizationDatasetGeneratorDialog::SetMinInclusions);
    connect(homogenization_dialog_,
            &HomogenizationDatasetGeneratorDialog::OnSetMinInclusions,
            homogenization_dialog_->min_inclusions, &QSpinBox::setValue);

    connect(homogenization_dialog_->max_inclusions,
            QOverload<int>::of(&QSpinBox::valueChanged), homogenization_dialog_,
            &HomogenizationDatasetGeneratorDialog::SetMaxInclusions);
    connect(homogenization_dialog_,
            &HomogenizationDatasetGeneratorDialog::OnSetMaxInclusions,
            homogenization_dialog_->max_inclusions, &QSpinBox::setValue);

    connect(homogenization_dialog_->cube_dimensions,
            QOverload<int>::of(&QSpinBox::valueChanged), homogenization_dialog_,
            &HomogenizationDatasetGeneratorDialog::SetCubeDimensions);
    connect(homogenization_dialog_,
            &HomogenizationDatasetGeneratorDialog::OnSetCubeDimensions,
            homogenization_dialog_->cube_dimensions, &QSpinBox::setValue);

    connect(homogenization_dialog_->E,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            homogenization_dialog_,
            &HomogenizationDatasetGeneratorDialog::SetMaterialYoungsModulus);
    connect(homogenization_dialog_->v,
            QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            homogenization_dialog_,
            &HomogenizationDatasetGeneratorDialog::SetMaterialPoissonsRatio);

    homogenization_dialog_->exec();
}

void DesignerController::SetNumberOfEntries(int value) {
    number_of_dataset_entries_ = value;
}
