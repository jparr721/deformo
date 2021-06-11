#include <QtWidgets/QApplication>

#include "Window.h"
#include <QtQuickControls2/QtQuickControls2>
#include <QDebug>

int main(int argc, char* argv[]) {
    QApplication a(argc, argv);
    Window w;
    w.show();
    return a.exec();
}
