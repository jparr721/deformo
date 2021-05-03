#include "Camera.h"

void Camera::Translate(const QVector3D& translation) {
  world.translate(-translation);
}

void Camera::Rotate(const QQuaternion& rotation) { world.rotate(rotation); }

const QMatrix4x4 Camera::Matrix() const { return world; }
