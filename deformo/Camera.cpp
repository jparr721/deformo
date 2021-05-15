#include "Camera.h"

void Camera::Translate(const QVector3D& t) {
  translation += t;
  world.translate(translation);
}

void Camera::Rotate(float angle, const QVector3D& axis) {
  world.rotate(QQuaternion::fromAxisAndAngle(axis, angle));
}
void Camera::Rotate(const QQuaternion& rotation) { world.rotate(rotation); }

const QMatrix4x4 Camera::Matrix() const { return world; }
