#include "Input.h"

InputState Input::KeyState(Qt::Key key) {
  const auto instance = GetKeyInstance(key);
  return instance != std::nullopt ? instance.value().second
                                  : InputState::Invalid;
}

bool Input::KeyTriggered(Qt::Key key) {
  return KeyState(key) == InputState::Triggered;
}

bool Input::KeyPressed(Qt::Key key) {
  return KeyState(key) == InputState::Pressed;
}

bool Input::KeyReleased(Qt::Key key) {
  return KeyState(key) == InputState::Released;
}

// ===========================================

InputState Input::MouseButtonState(Qt::MouseButton button) {
  const auto instance = GetMouseButtonInstance(button);
  return instance != std::nullopt ? instance.value().second
                                  : InputState::Invalid;
}

bool Input::MouseButtonTriggered(Qt::MouseButton key) {
  return MouseButtonState(key) == InputState::Triggered;
}

bool Input::MouseButtonPressed(Qt::MouseButton key) {
  return MouseButtonState(key) == InputState::Pressed;
}

bool Input::MouseButtonReleased(Qt::MouseButton key) {
  return MouseButtonState(key) == InputState::Released;
}

// ===========================================

QPoint Input::MousePosition() { return QCursor::pos(); }

QPoint Input::MouseDelta() { return mouse_pos_delta; }

// ===========================================

void Input::Update() {
  last_mouse_cursor_pos = mouse_cursor_pos;
  mouse_cursor_pos = MousePosition();
  mouse_pos_delta = mouse_cursor_pos - last_mouse_cursor_pos;

  UpdateInputContainer(keys);
  UpdateInputContainer(mouse_buttons);
}

void Input::Reset() {
  keys.clear();
  mouse_buttons.clear();
}

// ===========================================

void Input::RegisterKeyPress(int key) {
  const int key_idx = KeyExists(static_cast<Qt::Key>(key));

  if (key_idx == -1) {
    keys.push_back(
        KeyPressInstance(static_cast<Qt::Key>(key), InputState::Registered));
  }
}

void Input::RegisterKeyRelease(int key) {
  const int key_idx = KeyExists(static_cast<Qt::Key>(key));

  if (key_idx > -1) {
    keys[key_idx].second = InputState::Unregistered;
  }
}

// ===========================================

void Input::RegisterMouseButtonPress(Qt::MouseButton button) {
  const int mouse_button_idx = MouseButtonExists(button);

  if (mouse_button_idx == -1) {
    mouse_buttons.push_back(
        MouseButtonPressInstance(button, InputState::Registered));
  }
}

void Input::RegisterMouseButtonRelease(Qt::MouseButton button) {
  const int mouse_button_idx = MouseButtonExists(button);

  if (mouse_button_idx > -1) {
    mouse_buttons[mouse_button_idx].second = InputState::Unregistered;
  }
}

// ===========================================

int Input::KeyExists(Qt::Key key) {
  const auto it = std::find(keys.begin(), keys.end(), key);

  return it != keys.end() ? it - keys.begin() : -1;
}

int Input::MouseButtonExists(Qt::MouseButton button) {
  const auto it = std::find(mouse_buttons.begin(), mouse_buttons.end(), button);

  return it != mouse_buttons.end() ? it - mouse_buttons.begin() : -1;
}

// ===========================================

std::optional<KeyPressInstance> Input::GetKeyInstance(Qt::Key key) {
  const int key_idx = KeyExists(key);
  return key_idx == -1 ? std::nullopt : std::optional{keys.at(key_idx)};
}

std::optional<MouseButtonPressInstance> Input::GetMouseButtonInstance(
    Qt::MouseButton button) {
  const int mouse_button_idx = MouseButtonExists(button);
  return mouse_button_idx == -1
             ? std::nullopt
             : std::optional{mouse_buttons.at(mouse_button_idx)};
}
