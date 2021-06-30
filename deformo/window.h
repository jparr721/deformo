#pragma once

#include "WindowController.h"
#include <QtWidgets/QMainWindow>
#include <memory>

#include "ui_deformo.h"

class Window : public QMainWindow {
    Q_OBJECT

  public:
    explicit Window(QWidget* parent = Q_NULLPTR);

  private:
    Ui::deformoClass ui;
    std::shared_ptr<WindowController> controller_;
};
