#include "Mesh.h"

int Vertex::PositionOffset() { return offsetof(Vertex, position); }
int Vertex::ColorOffset() { return offsetof(Vertex, color); }
int Vertex::PositionSize() { return 3; }
int Vertex::ColorSize() { return 3; }
int Vertex::Stride() { return sizeof(Vertex); }
