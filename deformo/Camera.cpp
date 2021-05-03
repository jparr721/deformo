#include "Camera.h"

//void Camera::Translate(Eigen::Ref<const Eigen::Vector3d> translation) {
//  const Eigen::Affine3d trans(translation);
//  world *= trans.matrix();
//}
//
//void Camera::Rotate(Eigen::Quaterniond& rotation) {
//  const auto rot = rotation.normalized().toRotationMatrix().eval();
//  world *= rot;
//}
//
//const Eigen::Matrix4d Camera::Matrix() const { return world; }
