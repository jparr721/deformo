#pragma once

#include <QMatrix4x4>
#include <QVector3D>

class Camera {
 public:
  constexpr static float kTranslateSpeed = 0.05;
  constexpr static float kRotationSpeed = 0.005;
  const QVector3D kForward = QVector3D(0., 0., -kTranslateSpeed);
  const QVector3D kUp = QVector3D(0., kTranslateSpeed, 0.);
  const QVector3D kRight = QVector3D(kTranslateSpeed, 0., 0.);

  QMatrix4x4 world;
  QVector3D translation;
  QQuaternion rotation;

  Camera() { Reset(); }
  void Translate(const QVector3D& t);

  void Rotate(float angle, const QVector3D& axis);
  void Rotate(const QQuaternion& rot);
  void Reset();

  void Forward();
  void Backward();
  void Left();
  void Right();
  void Up();
  void Down();

  const QMatrix4x4 Matrix() const;

};
