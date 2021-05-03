#pragma once

#include <QMatrix4x4>
#include <QVector3D>

class Camera {
 public:
  const QVector3D kForward = QVector3D(0., 0., -1.);
  const QVector3D kUp = QVector3D(0., 1., 0.);
  const QVector3D kRight = QVector3D(1., 0., 0.);

  void Translate(const QVector3D& translation);
  void Rotate(const QQuaternion& rotation);

  const QMatrix4x4 Matrix() const;

 private:
  QMatrix4x4 world;
};
