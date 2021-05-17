#include "Input.h"

// ===========================================

bool Input::KeyPressed(Qt::Key key) {
  return keys.find(key) != keys.end() && keys[key];
}

bool Input::MouseButtonPressed(Qt::MouseButton button) {
  return mouse_buttons.find(button) != mouse_buttons.end() &&
         mouse_buttons[button];
}

void Input::RegisterKeyPress(int key) {
  const auto k = static_cast<Qt::Key>(key);
    keys[k] = true;
}

void Input::RegisterKeyRelease(int key) {
  const auto k = static_cast<Qt::Key>(key);
  keys[k] = false;
}

void Input::RegisterMouseButtonPress(Qt::MouseButton button) {}

void Input::RegisterMouseButtonRelease(Qt::MouseButton button) {}

QPoint Input::MousePosition() { return QCursor::pos(); }

QPoint Input::MouseDelta() { return mouse_pos_delta; }

// ===========================================
