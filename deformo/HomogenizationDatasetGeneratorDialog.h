#pragma once

#include "Numerics.h"
#include <QDialog>

QT_BEGIN_NAMESPACE
class QDialogButtonBox;
class QPushButton;
class QLineEdit;
class QSpinBox;
class QDoubleSpinBox;
class QLabel;
QT_END_NAMESPACE

class HomogenizationDatasetGeneratorDialog : public QDialog {
    Q_OBJECT

  public:
    explicit HomogenizationDatasetGeneratorDialog(QWidget* parent = nullptr);

    // Labels
    QLabel* min_inclusions_label;
    QLabel* max_inclusions_label;

    QLabel* min_inclusion_dimensions_label;
    QLabel* max_inclusion_dimensions_label;

    QLabel* cube_dimensions_label;

    QLabel* material_E_label;
    QLabel* material_v_label;

    // Buttons
    QPushButton* submit_button;

    // Button Box
    QDialogButtonBox* bottom_button_box;

    // Inclusion Height
    QSpinBox* min_inclusion_dimensions;
    QSpinBox* max_inclusion_dimensions;

    // N inclusions
    QSpinBox* min_inclusions;
    QSpinBox* max_inclusions;

    // Main cube dimensions
    QSpinBox* cube_dimensions;

    // Material coeff for mesh
    QDoubleSpinBox* E;
    QDoubleSpinBox* v;

    int min_inclusion_dimensions_ = 1;
    int max_inclusion_dimensions_ = 5;

    int min_inclusions_ = 1;
    int max_inclusions_ = 20;

    int cube_dimensions_ = 50;

    Real E_;
    Real v_;

  public slots:
    void SetMinInclusionDimensions(int value);
    void SetMaxInclusionDimensions(int value);

    void SetMinInclusions(int value);
    void SetMaxInclusions(int value);

    void SetCubeDimensions(int value);

    void SetMaterialYoungsModulus(Real value);
    void SetMaterialPoissonsRatio(Real value);

  signals:
    void OnSetMinInclusionDimensions(int value);
    void OnSetMaxInclusionDimensions(int value);

    void OnSetMinInclusions(int value);
    void OnSetMaxInclusions(int value);

    void OnSetCubeDimensions(int value);

  private:
    auto MakeSpinBox() -> QSpinBox*;
};
