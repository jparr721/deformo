#include <Eigen/Dense>
#include <QtWidgets/QApplication>
#include <iostream>

#include "Window.h"

int main(int argc, char *argv[]) {
  QApplication a(argc, argv);
  Window w;
  w.show();
  return a.exec();
}
