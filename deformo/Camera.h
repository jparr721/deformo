#pragma once

#include <QMatrix4x4>
#include <QVector3D>

class Camera {
 public:
  constexpr static float kTranslateSpeed = 0.05;
  const QVector3D kForward = QVector3D(0., 0., -1.);
  const QVector3D kUp = QVector3D(0., 1., 0.);
  const QVector3D kRight = QVector3D(1., 0., 0.);

  Camera() { Reset(); }
  void Translate(const QVector3D& t);
  void Rotate(float angle, const QVector3D& axis);
  void Rotate(const QQuaternion& rot);
  void Reset();

  const QMatrix4x4 Matrix() const;

 private:
  QMatrix4x4 world;
  QVector3D translation;
  QQuaternion rotation;
};
