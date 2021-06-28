#include "Input.h"

void Input::WheelEvent(QWheelEvent* e,
                       const std::shared_ptr<Camera<Real>>& camera) {
    const auto dt = e->delta();
    constexpr Real zoom_inverse = -0.01;
    camera->zoom(dt * zoom_inverse);
}

void Input::MouseMoveEvent(QMouseEvent* e,
                           const std::shared_ptr<Camera<Real>>& camera) {
    current_mouse_position.x() = e->x();
    current_mouse_position.y() = e->y();

    const Real x_diff = last_mouse_position.x() - current_mouse_position.x();
    const Real y_diff = last_mouse_position.y() - current_mouse_position.y();

    if (camera->is_panning) {
        camera->pan(x_diff, y_diff);
    }

    if (camera->is_rotating) {
        camera->rotate(x_diff, y_diff);
    }

    if (camera->is_zooming) {
        camera->zoom(y_diff);
    }

    last_mouse_position = current_mouse_position;
}

void Input::MousePressEvent(QMouseEvent* e,
                            const std::shared_ptr<Camera<Real>>& camera) {
    if (e->button() == Qt::MiddleButton)
        camera->is_panning = true;
    else if (e->button() == Qt::LeftButton)
        camera->is_rotating = true;
    else if (e->button() == Qt::RightButton)
        camera->is_zooming = true;

    last_mouse_position.x() = e->x();
    last_mouse_position.y() = e->y();
}

void Input::MouseReleaseEvent(QMouseEvent* e,
                              const std::shared_ptr<Camera<Real>>& camera) {
    if (e->button() == Qt::MiddleButton)
        camera->is_panning = false;
    else if (e->button() == Qt::LeftButton)
        camera->is_rotating = false;
    else if (e->button() == Qt::RightButton)
        camera->is_zooming = false;

    last_mouse_position.x() = e->x();
    last_mouse_position.y() = e->y();
}
