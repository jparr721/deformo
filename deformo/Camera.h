#pragma once

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Core>

class Camera {
 public:
  const Eigen::Vector3d kForward = Eigen::Vector3d(0., 0., -1.);
  const Eigen::Vector3d kUp = Eigen::Vector3d(0., 1., 0.);
  const Eigen::Vector3d kRight = Eigen::Vector3d(1., 0., 0.);

  void Translate(Eigen::Ref<const Eigen::Vector3d> translation);
  void Rotate(Eigen::Quaterniond& rotation);

  const Eigen::Matrix4d Matrix() const;

 private:
  Eigen::Matrix4d world;
};
