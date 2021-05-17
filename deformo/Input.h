#pragma once

#include <QCursor>
#include <QPoint>
#include <Qt>
#include <algorithm>
#include <optional>
#include <utility>
#include <vector>
#include <unordered_map>

class Input {
 public:
  QPoint mouse_cursor_pos;
  QPoint last_mouse_cursor_pos;
  QPoint mouse_pos_delta;

  std::unordered_map<Qt::Key, bool> keys;
  std::unordered_map<Qt::MouseButton, bool> mouse_buttons;

  // State Checks
  [[nodiscard]] bool KeyPressed(Qt::Key key);
  [[nodiscard]] bool MouseButtonPressed(Qt::MouseButton button);

  // Keyboard State Changes
  void RegisterKeyPress(int key);
  void RegisterKeyRelease(int key);

  // Mouse State Changes
  void RegisterMouseButtonPress(Qt::MouseButton button);
  void RegisterMouseButtonRelease(Qt::MouseButton button);

  // Mouse Position Handlers
  [[nodiscard]] QPoint MousePosition();
  [[nodiscard]] QPoint MouseDelta();
};
