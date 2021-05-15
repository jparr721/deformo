#pragma once

#include <QCursor>
#include <QPoint>
#include <Qt>
#include <algorithm>
#include <optional>
#include <utility>
#include <vector>

enum class InputState {
  Invalid = 0x00,
  Registered,
  Unregistered,
  Triggered,
  Pressed,
  Released
};

template <typename InputType>
struct InputInstance : std::pair<InputType, InputState> {
  inline InputInstance(InputType input_type)
      : std::pair<InputType, InputState>(input_type, InputState::Invalid) {}
  inline InputInstance(InputType input_type, InputState state)
      : std::pair<InputType, InputState>(input_type, state) {}
  inline bool operator==(const InputInstance& rhs) const {
    return first == rhs.first;
  }
};

using KeyPressInstance = InputInstance<Qt::Key>;
using MouseButtonPressInstance = InputInstance<Qt::MouseButton>;

class Input {
 public:
  QPoint mouse_cursor_pos;
  QPoint last_mouse_cursor_pos;
  QPoint mouse_pos_delta;

  std::vector<KeyPressInstance> keys;
  std::vector<MouseButtonPressInstance> mouse_buttons;

  // Keyboard Input Handlers
  InputState KeyState(Qt::Key key);
  [[nodiscard]] bool KeyTriggered(Qt::Key key);
  [[nodiscard]] bool KeyPressed(Qt::Key key);
  [[nodiscard]] bool KeyReleased(Qt::Key key);

  // Mouse Input Handlers
  InputState MouseButtonState(Qt::MouseButton button);
  [[nodiscard]] bool MouseButtonTriggered(Qt::MouseButton button);
  [[nodiscard]] bool MouseButtonPressed(Qt::MouseButton button);
  [[nodiscard]] bool MouseButtonReleased(Qt::MouseButton button);

  // Mouse Position Handlers
  [[nodiscard]] QPoint MousePosition();
  [[nodiscard]] QPoint MouseDelta();

  // Keyboard State Changes
  void RegisterKeyPress(int key);
  void RegisterKeyRelease(int key);

  // Mouse State Changes
  void RegisterMouseButtonPress(Qt::MouseButton button);
  void RegisterMouseButtonRelease(Qt::MouseButton button);

  // State Handlers
  template <typename InputInstanceType>
  inline void UpdateInputState(InputInstanceType& instance) {
    switch (instance.second) {
      case InputState::Registered:
        instance.second = InputState::Triggered;
        break;
      case InputState::Triggered:
        instance.second = InputState::Pressed;
        break;
      case InputState::Unregistered:
        instance.second = InputState::Released;
        break;
      default:
        break;
    }
  }

  // Container Management
  template <typename InputInstanceContainer>
  inline void UpdateInputContainer(InputInstanceContainer& container) {
    const auto should_remove = [](auto instance) -> bool {
      return instance.second == InputState::Released;
    };

    // Remove buttons that aren't being pressed.
    container.erase(
        std::remove_if(container.begin(), container.end(), should_remove),
        container.end());

    // Update states of buttons which have been pressed.
    for (auto& value : container) {
      UpdateInputState(value);
    }
  }

 private:
  void Update();
  void Reset();

  // Utilities for Key management
  [[nodiscard]] int KeyExists(Qt::Key key);
  [[nodiscard]] int MouseButtonExists(Qt::MouseButton button);

  // Mutable reference returns for changing key states
  std::optional<KeyPressInstance> GetKeyInstance(Qt::Key key);
  std::optional<MouseButtonPressInstance> GetMouseButtonInstance(
      Qt::MouseButton button);
};
