#include "Camera.h"

void Camera::Translate(const QVector3D& t) {
  world.translate(-t);
  translation = QVector3D(0.f, 0.f, 0.f);
}

void Camera::Forward() {
  translation += kForward;
  Translate(translation);
}

void Camera::Backward() {
  translation -= kForward;
  Translate(translation);
}

void Camera::Left() {
  translation -= kRight;
  Translate(translation);
}

void Camera::Right() {
  translation += kRight;
  Translate(translation);
}

void Camera::Up() {
  translation += kUp;
  Translate(translation);
}

void Camera::Down() {
  translation -= kUp;
  Translate(translation);
}

void Camera::Rotate(float angle, const QVector3D& axis) {
  const auto rot = QQuaternion::fromAxisAndAngle(axis, angle);
  Rotate(rot);
}
void Camera::Rotate(const QQuaternion& rot) {
  rotation *= rot;
  world.rotate(rotation);
}

const QMatrix4x4 Camera::Matrix() const { return world; }

void Camera::Reset() {
  QMatrix4x4 new_world;
  new_world.perspective(45.f, 4.f / 3.f, 0.f, 2000.f);
  world = new_world;
  translation = QVector3D(0.f, 0.f, 0.f);
  rotation = QQuaternion();
}
