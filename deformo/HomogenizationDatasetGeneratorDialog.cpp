#include "HomogenizationDatasetGeneratorDialog.h"
#include <QDialogButtonBox>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QSpinBox>
#include <QVBoxLayout>
#include <iostream>

HomogenizationDatasetGeneratorDialog::HomogenizationDatasetGeneratorDialog(
    QWidget* parent)
    : QDialog(parent) {
    // Inclusions ================================
    min_inclusions_label = new QLabel("Min number of inclusions:");
    min_inclusions = MakeSpinBox();
    min_inclusions_label->setBuddy(min_inclusions);

    max_inclusions_label = new QLabel("Max number of inclusions:");
    max_inclusions = MakeSpinBox();
    max_inclusions->setValue(20);
    max_inclusions_label->setBuddy(max_inclusions);

    // Dims =========================================
    min_inclusion_dimensions_label = new QLabel("Min inclusion size (nxnxn):");
    min_inclusion_dimensions = MakeSpinBox();
    min_inclusion_dimensions_label->setBuddy(min_inclusion_dimensions);

    max_inclusion_dimensions_label = new QLabel("Max inclusion size (nxnxn):");
    max_inclusion_dimensions = MakeSpinBox();
    max_inclusion_dimensions->setValue(5);
    max_inclusion_dimensions_label->setBuddy(max_inclusion_dimensions);

    // Cube =========================================
    cube_dimensions_label = new QLabel("Cuboid Size (nxnxn):");
    cube_dimensions = MakeSpinBox();
    cube_dimensions->setValue(50);
    cube_dimensions_label->setBuddy(cube_dimensions);

    // Material =========================================
    material_E_label = new QLabel("Young's Modulus");
    E = new QDoubleSpinBox;
    E->setMinimum(1000);
    E->setMaximum(1000000000000000);
    material_E_label->setBuddy(E);

    material_v_label = new QLabel("Poisson's Ratio");
    v = new QDoubleSpinBox;
    v->setMinimum(0.00000001);
    v->setMaximum(0.5);
    material_v_label->setBuddy(v);

    // Submission Buttons ================================
    submit_button = new QPushButton("Submit");
    // Automatically submit when the user presses enter.
    submit_button->setDefault(true);

    // Layout ============================================
    QVBoxLayout* main_layout = new QVBoxLayout;
    QGridLayout* parameters_layout = new QGridLayout;
    parameters_layout->addWidget(min_inclusions_label, 0, 0, 1, 1);
    parameters_layout->addWidget(max_inclusions_label, 0, 1, 1, 1);
    parameters_layout->addWidget(min_inclusions, 1, 0, 1, 1);
    parameters_layout->addWidget(max_inclusions, 1, 1, 1, 1);

    parameters_layout->addWidget(min_inclusion_dimensions_label, 2, 0, 1, 1);
    parameters_layout->addWidget(max_inclusion_dimensions_label, 2, 1, 1, 1);
    parameters_layout->addWidget(min_inclusion_dimensions, 3, 0, 1, 1);
    parameters_layout->addWidget(max_inclusion_dimensions, 3, 1, 1, 1);

    parameters_layout->addWidget(cube_dimensions_label, 4, 0, 1, 1);
    parameters_layout->addWidget(cube_dimensions, 5, 0, 1, 1);

    parameters_layout->addWidget(material_E_label, 6, 0, 1, 1);
    parameters_layout->addWidget(material_v_label, 6, 1, 1, 1);
    parameters_layout->addWidget(E, 7, 0, 1, 1);
    parameters_layout->addWidget(v, 7, 1, 1, 1);

    parameters_layout->addWidget(submit_button, 8, 1);

    setLayout(parameters_layout);
    setMinimumSize(300, 300);
    setWindowTitle("Homogenization Generator");
}

auto HomogenizationDatasetGeneratorDialog::MakeSpinBox() -> QSpinBox* {
    auto spin_box = new QSpinBox;
    spin_box->setMinimum(1);
    spin_box->setMaximum(10000000);
    return spin_box;
}

void HomogenizationDatasetGeneratorDialog::SetMinInclusionDimensions(
    int value) {
    min_inclusion_dimensions_ = value;
    OnSetMinInclusionDimensions(value);
}

void HomogenizationDatasetGeneratorDialog::SetMaxInclusionDimensions(
    int value) {
    max_inclusion_dimensions_ = value;
    OnSetMaxInclusionDimensions(value);
}

void HomogenizationDatasetGeneratorDialog::SetMinInclusions(int value) {
    min_inclusions_ = value;
    OnSetMinInclusions(value);
}

void HomogenizationDatasetGeneratorDialog::SetMaxInclusions(int value) {
    max_inclusions_ = value;
    OnSetMaxInclusions(value);
}

void HomogenizationDatasetGeneratorDialog::SetCubeDimensions(int value) {
    cube_dimensions_ = value;
    OnSetCubeDimensions(value);
}

void HomogenizationDatasetGeneratorDialog::SetMaterialYoungsModulus(
    Real value) {
    E_ = value;
}

void HomogenizationDatasetGeneratorDialog::SetMaterialPoissonsRatio(
    Real value) {
    v_ = value;
}
