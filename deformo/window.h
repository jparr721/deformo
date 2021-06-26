#pragma once

#include <QtWidgets/QMainWindow>
#include "WindowController.h"

#include "ui_deformo.h"

class Window : public QMainWindow {
    Q_OBJECT

  public:
    Window(QWidget* parent = Q_NULLPTR);
    ~Window();

  private:
    Ui::deformoClass ui;
    WindowController* window_controller_;
};
