#include "Camera.h"

void Camera::Translate(const QVector3D& t) {
  translation += t;
  world.translate(kTranslateSpeed * translation);
}

void Camera::Rotate(float angle, const QVector3D& axis) {
  const auto rot = QQuaternion::fromAxisAndAngle(axis, angle);
  Rotate(rot);
}
void Camera::Rotate(const QQuaternion& rot) {
  rotation *= rot;
  world.rotate(rotation);
}

const QMatrix4x4 Camera::Matrix() const {
  return world;
}

void Camera::Reset() {
  QMatrix4x4 new_world;
  new_world.perspective(45.f, 4.f / 3.f, 0.f, 2000.f);
  world = new_world;
}
