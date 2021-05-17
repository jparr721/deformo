#include "Input.h"
#include <iostream>

void Input::Update() {
  last_mouse_cursor_pos = mouse_cursor_pos;
  mouse_cursor_pos = QCursor::pos();
  mouse_pos_delta = mouse_cursor_pos - last_mouse_cursor_pos;
}

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

void Input::RegisterMouseButtonPress(Qt::MouseButton button) {
  mouse_buttons[button] = true;
}

void Input::RegisterMouseButtonRelease(Qt::MouseButton button) {
  mouse_buttons[button] = false;
}

QPoint Input::MousePosition() { return QCursor::pos(); }

QPoint Input::MouseDelta() { return mouse_pos_delta; }
