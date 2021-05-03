#pragma once
#include <QVector3D>

class VVertex
{
public:
  // Constructors
  Q_DECL_CONSTEXPR VVertex();
  Q_DECL_CONSTEXPR explicit VVertex(const QVector3D &position);
  Q_DECL_CONSTEXPR VVertex(const QVector3D &position, const QVector3D &color);

  // Accessors / Mutators
  Q_DECL_CONSTEXPR const QVector3D& position() const;
  Q_DECL_CONSTEXPR const QVector3D& color() const;
  void setPosition(const QVector3D& position);
  void setColor(const QVector3D& color);

  // OpenGL Helpers
  static const int PositionTupleSize = 3;
  static const int ColorTupleSize = 3;
  static Q_DECL_CONSTEXPR int positionOffset();
  static Q_DECL_CONSTEXPR int colorOffset();
  static Q_DECL_CONSTEXPR int stride();

  QVector3D m_position;
  QVector3D m_color;
};

/*******************************************************************************
 * Inline Implementation
 ******************************************************************************/

// Note: Q_MOVABLE_TYPE means it can be memcpy'd.
Q_DECLARE_TYPEINFO(VVertex, Q_MOVABLE_TYPE);

// Constructors
Q_DECL_CONSTEXPR inline VVertex::VVertex() {}
Q_DECL_CONSTEXPR inline VVertex::VVertex(const QVector3D &position) : m_position(position) {}
Q_DECL_CONSTEXPR inline VVertex::VVertex(const QVector3D &position, const QVector3D &color) : m_position(position), m_color(color) {}

// Accessors / Mutators
Q_DECL_CONSTEXPR inline const QVector3D& VVertex::position() const { return m_position; }
Q_DECL_CONSTEXPR inline const QVector3D& VVertex::color() const { return m_color; }
void inline VVertex::setPosition(const QVector3D& position) { m_position = position; }
void inline VVertex::setColor(const QVector3D& color) { m_color = color; }

// OpenGL Helpers
Q_DECL_CONSTEXPR inline int VVertex::positionOffset() { return offsetof(VVertex, m_position); }
Q_DECL_CONSTEXPR inline int VVertex::colorOffset() { return offsetof(VVertex, m_color); }
Q_DECL_CONSTEXPR inline int VVertex::stride() { return sizeof(VVertex); }
