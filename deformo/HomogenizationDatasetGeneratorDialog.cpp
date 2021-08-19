#include "HomogenizationDatasetGeneratorDialog.h"
#include <QDialogButtonBox>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QSpinBox>
#include <QVBoxLayout>

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

    parameters_layout->addWidget(submit_button, 6, 1);

    setLayout(parameters_layout);
    setWindowTitle("Homogenization Generator");
}

auto HomogenizationDatasetGeneratorDialog::MakeSpinBox() -> QSpinBox* {
    auto spin_box = new QSpinBox;
    spin_box->setMinimum(1);
    spin_box->setMaximum(10000000);
    return spin_box;
}

auto HomogenizationDatasetGeneratorDialog::ConnectWidgets() -> void {
    connect(min_inclusion_dimensions,
            QOverload<int>::of(&QSpinBox::valueChanged), this,
            &HomogenizationDatasetGeneratorDialog::SetMinInclusionDimensions);
    connect(this,
            &HomogenizationDatasetGeneratorDialog::OnSetMinInclusionDimensions,
            min_inclusion_dimensions, &QSpinBox::setValue);

    connect(max_inclusion_dimensions,
            QOverload<int>::of(&QSpinBox::valueChanged), this,
            &HomogenizationDatasetGeneratorDialog::SetMaxInclusionDimensions);
    connect(this,
            &HomogenizationDatasetGeneratorDialog::OnSetMaxInclusionDimensions,
            max_inclusion_dimensions, &QSpinBox::setValue);

    connect(min_inclusions, QOverload<int>::of(&QSpinBox::valueChanged), this,
            &HomogenizationDatasetGeneratorDialog::SetMinInclusions);
    connect(this, &HomogenizationDatasetGeneratorDialog::SetMinInclusions,
            min_inclusions, &QSpinBox::setValue);

    connect(max_inclusions, QOverload<int>::of(&QSpinBox::valueChanged), this,
            &HomogenizationDatasetGeneratorDialog::SetMaxInclusions);
    connect(this, &HomogenizationDatasetGeneratorDialog::SetMaxInclusions,
            max_inclusions, &QSpinBox::setValue);

    connect(cube_dimensions, QOverload<int>::of(&QSpinBox::valueChanged), this,
            &HomogenizationDatasetGeneratorDialog::SetCubeDimensions);
    connect(this, &HomogenizationDatasetGeneratorDialog::SetCubeDimensions,
            cube_dimensions, &QSpinBox::setValue);

    connect(submit_button, &QPushButton::clicked, this,
            &HomogenizationDatasetGeneratorDialog::done);
}

void HomogenizationDatasetGeneratorDialog::SetMinInclusionDimensions(
    int value) {
    min_inclusion_dimensions_ = value;
    OnSetMinInclusionDimensions(value);
}
//
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
