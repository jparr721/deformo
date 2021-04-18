#pragma once

#include <QtWidgets/QMainWindow>

#include "ui_deformo.h"

class Window : public QMainWindow {
  Q_OBJECT

 public:
  Window(QWidget *parent = Q_NULLPTR);

 private:
  Ui::deformoClass ui;
};
