#pragma once

#include "Camera.h"

#include <Eigen/Dense>
#include <QMouseEvent>
#include <QWheelEvent>

class Input {
  public:
    Eigen::Vector2i current_mouse_position = Eigen::Vector2i::Zero();
    Eigen::Vector2i last_mouse_position = Eigen::Vector2i::Zero();

    // Mouse Wheel
    void WheelEvent(QWheelEvent* e,
                    const std::shared_ptr<Camera<Real>>& camera);

    // Mouse Buttons
    void MouseMoveEvent(QMouseEvent* e,
                        const std::shared_ptr<Camera<Real>>& camera);
    void MousePressEvent(QMouseEvent* e,
                         const std::shared_ptr<Camera<Real>>& camera);
    void MouseReleaseEvent(QMouseEvent* e,
                           const std::shared_ptr<Camera<Real>>& camera);
};
