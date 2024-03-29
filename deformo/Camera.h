#pragma once

#include "Numerics.h"

template <typename Real> class Camera {
  public:
    bool is_panning = false;
    bool is_rotating = false;
    bool is_zooming = false;

    Camera();
    virtual ~Camera();

    auto Resize(std::size_t width, std::size_t height) -> bool;
    auto Resize(std::size_t width, std::size_t height, Real near_plane,
                Real far_plane) -> bool;

    auto Rotate(Real du, Real dv) -> bool;
    auto Pan(Real du, Real dv) -> bool;
    auto Zoom(Real dr) -> bool;

    void SetPerspective(Real fov, std::size_t width, std::size_t height,
                        Real nearPlane, Real farPlane);
    void SetFrustum(Real left, Real right, Real top, Real bottom,
                    Real nearPlane, Real farPlane);

    void SetPosition(const Vector3r& position);
    void SetPosition(Real x, Real y, Real z);
    void AddPosition(const Vector3r& displacement);
    void AddPosition(Real dx, Real dy, Real dz);

    void SetRotation(Real theta, Real phi);
    void SetHorizontalRotation(Real theta);
    void SetVerticalRotation(Real phi);
    void SetSphericalPosition(Real r, Real theta, Real phi);
    void AddRadius(Real dRadius);
    void SetRadius(Real radius);
    void SetFieldOfView(Real fov);
    void SetNearPlane(Real nearPlane);
    void SetFarPlane(Real farPlane);
    void SetMinRadius(Real minRadius);
    void SetMaxRadius(Real maxRadius);

    void SetRotationSensitivity(Real sen);
    void SetPanSensitivity(Real sen);
    void SetZoomSensitivity(Real sen);
    void SetScopeSensitivity(Real sen);

    void Reset();
    void ResetPlanes();
    void ResetView();
    void ResetMatrices();
    void ResetSensitivities();

    auto GetRotationSensitivity() const -> Real;
    auto GetPanSensitivity() const -> Real;
    auto GetZoomSensitivity() const -> Real;
    auto GetScopeSensitivity() const -> Real;

    auto GetViewDirection() -> Vector3r;
    auto GetRightDirection() -> Vector3r;
    auto GetLeftDirection() -> Vector3r;
    auto GetUpDirection() -> Vector3r;
    auto GetDownDirection() -> Vector3r;

    auto GetEye() const -> const Vector3r&;
    auto GetLookAt() const -> const Vector3r&;
    auto GetUp() const -> const Vector3r&;

    auto GetWidth() const -> std::size_t;
    auto GetHeight() const -> std::size_t;
    auto GetNear() const -> const Real&;
    auto GetFar() const -> const Real&;

    auto GetViewMatrix() -> const Matrix4r&;
    auto GetProjectionMatrix() -> const Matrix4r&;

    auto ToViewMatrix() -> Matrix4r;
    auto ToProjectionMatrix() -> Matrix4r;

    static auto LookAt(const Vector3r& eye, const Vector3r& lookAt,
                       const Vector3r& up) -> Matrix4r;
    static auto LookAt(Real eyex, Real eyey, Real eyez, Real atx, Real aty,
                       Real atz, Real upx, Real upy, Real upz) -> Matrix4r;

    static auto PerspectiveMatrix(Real fov, std::size_t width,
                                  std::size_t height, Real nearPlane,
                                  Real farPlane) noexcept -> Matrix4r;
    static auto OrthographicMatrix(Real left, Real right, Real bottom, Real top,
                                   Real nearPlane, Real farPlane) noexcept
        -> Matrix4r;

    static auto SphereicalToCartesian(Real r, Real theta, Real phi) -> Vector3r;
    static auto SphereicalToCartesian_dTheta(Real r, Real theta, Real phi)
        -> Vector3r;
    static auto SphereicalToCartesian_dPhi(Real r, Real theta, Real phi)
        -> Vector3r;
    static auto SphereicalToCartesian_dPhiCrossdTheta(Real r, Real theta,
                                                      Real phi) -> Vector3r;

  protected:
    void Compile() noexcept;

  protected:
    Matrix4r view_matrix_;
    Matrix4r projection_matrix_;
    Real near_plane_, far_plane_;

    std::size_t width_, height_;
    Vector3r eye_;
    Vector3r look_at_;
    Vector3r up_;

    Real fov_, aspect_ratio_;
    Real min_radius_, max_radius_;

    Real r_, theta_, phi_;
    Vector3r displacement_;

    Real pan_sensitivity_;
    Real zoom_sensitivity_;
    Real rotation_sensitivity_;
    Real scope_sensitivity_;
};

constexpr Real kPi = static_cast<Real>(3.14159);
constexpr Real kCameraEpsilon = static_cast<Real>(0.00001);
constexpr Real kDefaultNearPlane = static_cast<Real>(0.1);
constexpr Real kDefaultFarPlane = static_cast<Real>(1000.0);
constexpr Real kDefaultRotationSensitivity = static_cast<Real>(0.01);
constexpr Real kDefaultPanSensitivity = static_cast<Real>(0.01);
constexpr Real kDefaultZoomSensitivity = static_cast<Real>(0.1);
constexpr Real kDefaultScopeSensitivity = static_cast<Real>(0.01);
constexpr Real kDefaultRadius = static_cast<Real>(1);
constexpr Real kDefaultTheta = static_cast<Real>(1.57079632679);
constexpr Real kDefaultPhi = static_cast<Real>(1.57079632679);
constexpr Real kDefaultFov = static_cast<Real>(65.0);
constexpr Real kDefaultAspectRatio = 4.f / 3.f;
constexpr Real kDefaultMinRadius = static_cast<Real>(0.1);
constexpr Real kDefaultMaxRadius = static_cast<Real>(1000.0);
constexpr Real kMaxFov = static_cast<Real>(120.0);
constexpr Real kMinFov = static_cast<Real>(10.0);

template <typename Real> Camera<Real>::Camera() {
    view_matrix_ = Matrix4r::Identity();
    projection_matrix_ = Matrix4r::Identity();

    eye_ = Vector3r::UnitZ();
    look_at_ = Vector3r::Zero();
    up_ = Vector3r::UnitY();

    near_plane_ = Real(kDefaultNearPlane);
    far_plane_ = Real(kDefaultFarPlane);

    rotation_sensitivity_ = Real(kDefaultRotationSensitivity);
    pan_sensitivity_ = Real(kDefaultPanSensitivity);
    zoom_sensitivity_ = Real(kDefaultZoomSensitivity);
    scope_sensitivity_ = Real(kDefaultScopeSensitivity);

    r_ = Real(1);
    theta_ = Real(kDefaultTheta);
    phi_ = Real(kDefaultPhi);
    fov_ = Real(kDefaultFov);
    aspect_ratio_ = Real(kDefaultAspectRatio);
    min_radius_ = Real(kDefaultMinRadius);
    max_radius_ = Real(kDefaultMaxRadius);
    displacement_ = Vector3r::Zero();
}

template <typename Real> Camera<Real>::~Camera() {}

template <typename Real>
auto Camera<Real>::Resize(std::size_t width, std::size_t height) -> bool {
    SetPerspective(fov_, width, height, near_plane_, far_plane_);
    Compile();
    return true;
}

template <typename Real>
auto Camera<Real>::Resize(std::size_t width, std::size_t height,
                          Real near_plane, Real far_plane) -> bool {
    SetPerspective(fov_, width, height, near_plane, far_plane);
    Compile();
    return true;
}

template <typename Real> auto Camera<Real>::Rotate(Real du, Real dv) -> bool {
    theta_ -= (du * rotation_sensitivity_);
    phi_ += (dv * rotation_sensitivity_);
    Compile();
    return true;
}

template <typename Real> auto Camera<Real>::Pan(Real du, Real dv) -> bool {
    Vector3r uDir = GetLeftDirection();
    Vector3r vDir = GetDownDirection();

    Vector3r uDisp = (du * pan_sensitivity_) * uDir;
    Vector3r vDisp = (dv * pan_sensitivity_) * vDir;
    Vector3r panDisp = uDisp + vDisp;

    displacement_ += panDisp;
    Compile();
    return true;
}

template <typename Real> auto Camera<Real>::Zoom(Real dr) -> bool {
    if (r_ + dr > max_radius_) {
        r_ = max_radius_;
        Compile();
        return false;
    }
    if (r_ + dr < min_radius_) {
        r_ = min_radius_;
        Compile();
        return false;
    }
    r_ += dr;
    Compile();
    return true;
}

template <typename Real>
void Camera<Real>::SetPerspective(Real fov, std::size_t width,
                                  std::size_t height, Real nearPlane,
                                  Real farPlane) {
    if (fov > Real(kMaxFov))
        fov_ = Real(kMaxFov);
    else if (fov < Real(kMinFov))
        fov_ = Real(kMinFov);
    else
        fov_ = fov;

    width_ = width;
    height_ = height;
    near_plane_ = nearPlane;
    far_plane_ = farPlane;

    Real aspectRatio = static_cast<Real>(width) / static_cast<Real>(height);
    Real ymax = nearPlane * std::tan(fov * Real(kPi) / Real(360));
    Real xmax = ymax * aspectRatio;
    SetFrustum(-xmax, xmax, ymax, -ymax, nearPlane, farPlane);
}

template <typename Real>
void Camera<Real>::SetFrustum(Real left, Real right, Real top, Real bottom,
                              Real nearPlane, Real farPlane) {
    projection_matrix_.setZero();

    Real temp, temp2, temp3, temp4;
    temp = Real(2) * nearPlane;
    temp2 = right - left;
    temp3 = top - bottom;
    temp4 = farPlane - nearPlane;

    projection_matrix_(0) = temp / temp2;
    projection_matrix_(5) = temp / temp3;
    projection_matrix_(8) = (right + left) / temp2;
    projection_matrix_(9) = (top + bottom) / temp3;
    projection_matrix_(10) = (-farPlane - nearPlane) / temp4;
    projection_matrix_(11) = Real(-1);
    projection_matrix_(14) = (-temp * farPlane) / temp4;
}

template <typename Real> void Camera<Real>::SetPosition(const Vector3r& pos) {
    position.x() = pos.x();
    position.y() = pos.y();
    position.z() = pos.z();
    Compile();
}

template <typename Real>
void Camera<Real>::SetPosition(Real x, Real y, Real z) {
    position.x() = x;
    position.y() = y;
    position.z() = z;
    Compile();
}

template <typename Real>
void Camera<Real>::AddPosition(const Vector3r& displacement) {
    position += displacement;
    Compile();
}

template <typename Real>
void Camera<Real>::AddPosition(Real dx, Real dy, Real dz) {
    displacement_.x() += dx;
    displacement_.y() += dy;
    displacement_.z() += dz;
    Compile();
}

template <typename Real> void Camera<Real>::SetRotation(Real theta, Real phi) {
    theta_ = theta;
    phi_ = phi;
    Compile();
}

template <typename Real> void Camera<Real>::SetHorizontalRotation(Real theta) {
    theta_ = theta;
    Compile();
}

template <typename Real> void Camera<Real>::SetVerticalRotation(Real phi) {
    phi_ = phi;
    Compile();
}

template <typename Real>
void Camera<Real>::SetSphericalPosition(Real r, Real theta, Real phi) {
    r_ = r;
    theta_ = theta;
    phi_ = phi;
    Compile();
}

template <typename Real> void Camera<Real>::AddRadius(Real dRadius) {
    r_ += dRadius;
    Compile();
}

template <typename Real> void Camera<Real>::SetRadius(Real radius) {
    r_ = radius;
    Compile();
}

template <typename Real> void Camera<Real>::SetFieldOfView(Real fov) {
    if (fov > Real(kMaxFov))
        fov_ = Real(kMaxFov);
    else if (fov < Real(kMinFov))
        fov_ = Real(kMinFov);
    else
        fov_ = fov;
    SetPerspective(fov_, aspect_ratio_, near_plane_, far_plane_);
}

template <typename Real> void Camera<Real>::SetNearPlane(Real nearPlane) {
    near_plane_ = nearPlane;
    SetPerspective(fov_, aspect_ratio_, nearPlane, far_plane_);
}

template <typename Real> void Camera<Real>::SetFarPlane(Real farPlane) {
    far_plane_ = farPlane;
    SetPerspective(fov_, aspect_ratio_, near_plane_, farPlane);
}

template <typename Real> void Camera<Real>::SetMinRadius(Real minRadius) {
    min_radius_ = minRadius;
}

template <typename Real> void Camera<Real>::SetMaxRadius(Real maxRadius) {
    max_radius_ = maxRadius;
}

template <typename Real> void Camera<Real>::SetRotationSensitivity(Real sen) {
    rotation_sensitivity_ = sen;
}

template <typename Real> void Camera<Real>::SetPanSensitivity(Real sen) {
    pan_sensitivity_ = sen;
}

template <typename Real> void Camera<Real>::SetZoomSensitivity(Real sen) {
    zoom_sensitivity_ = sen;
}

template <typename Real> void Camera<Real>::SetScopeSensitivity(Real sen) {
    scope_sensitivity_ = sen;
}

template <typename Real> void Camera<Real>::Reset() {
    ResetMatrices();
    ResetPlanes();
    ResetView();
    ResetSensitivities();
    Compile();
}

template <typename Real> void Camera<Real>::ResetPlanes() {
    near_plane_ = Real(kDefaultNearPlane);
    far_plane_ = Real(kDefaultFarPlane);
}

template <typename Real> void Camera<Real>::ResetView() {
    eye_ = Vector3r::UnitZ();
    look_at_ = Vector3r::Zero();
    up_ = Vector3r::UnitY();
}

template <typename Real> void Camera<Real>::ResetMatrices() {
    view_matrix_ = Matrix4r::Identity();
    projection_matrix_ = Matrix4r::Identity();
}

template <typename Real> void Camera<Real>::ResetSensitivities() {
    rotation_sensitivity_ = Real(kDefaultRotationSensitivity);
    pan_sensitivity_ = Real(kDefaultPanSensitivity);
    zoom_sensitivity_ = Real(kDefaultZoomSensitivity);
    scope_sensitivity_ = Real(kDefaultScopeSensitivity);
}

template <typename Real>
auto Camera<Real>::GetRotationSensitivity() const -> Real {
    return rotation_sensitivity_;
}

template <typename Real> auto Camera<Real>::GetPanSensitivity() const -> Real {
    return pan_sensitivity_;
}

template <typename Real> auto Camera<Real>::GetZoomSensitivity() const -> Real {
    return zoom_sensitivity_;
}

template <typename Real>
auto Camera<Real>::GetScopeSensitivity() const -> Real {
    return scope_sensitivity_;
}

template <typename Real>
inline auto Camera<Real>::GetViewDirection() -> Vector3r {
    Compile();
    return (look_at_ - eye_).normalized();
}

template <typename Real>
inline auto Camera<Real>::GetRightDirection() -> Vector3r {
    Compile();
    Vector3r dir = (look_at_ - eye_).normalized();
    return (up_.cross(dir)).normalized();
}

template <typename Real>
inline auto Camera<Real>::GetLeftDirection() -> Vector3r {
    Compile();
    Vector3r dir = (look_at_ - eye_).normalized();
    return (dir.cross(up_)).normalized();
}

template <typename Real> inline auto Camera<Real>::GetUpDirection() -> Vector3r {
    Compile();
    return up_;
}

template <typename Real>
inline auto Camera<Real>::GetDownDirection() -> Vector3r {
    Compile();
    return -up_;
}

template <typename Real> auto Camera<Real>::GetEye() const -> const Vector3r& {
    return eye_;
}

template <typename Real>
auto Camera<Real>::GetLookAt() const -> const Vector3r& {
    return look_at_;
}

template <typename Real> auto Camera<Real>::GetUp() const -> const Vector3r& {
    return up_;
}

template <typename Real> auto Camera<Real>::GetWidth() const -> std::size_t {
    return width_;
}

template <typename Real> auto Camera<Real>::GetHeight() const -> std::size_t {
    return height_;
}

template <typename Real> auto Camera<Real>::GetNear() const -> const Real& {
    return near_plane_;
}

template <typename Real> auto Camera<Real>::GetFar() const -> const Real& {
    return far_plane_;
}

template <typename Real> auto Camera<Real>::GetViewMatrix() -> const Matrix4r& {
    Compile();
    return view_matrix_;
}

template <typename Real>
auto Camera<Real>::GetProjectionMatrix() -> const Matrix4r& {
    Compile();
    return projection_matrix_;
}

template <typename Real> auto Camera<Real>::ToViewMatrix() -> Matrix4r {
    Compile();
    return view_matrix_;
}

template <typename Real> auto Camera<Real>::ToProjectionMatrix() -> Matrix4r {
    Compile();
    return projection_matrix_;
}

template <typename Real>
auto Camera<Real>::LookAt(const Vector3r& eye, const Vector3r& lookAt,
                          const Vector3r& up) -> Matrix4r {
    return Camera<Real>::LookAt(eye.x(), eye.y(), eye.z(), lookAt.x(),
                                lookAt.y(), lookAt.z(), up.x(), up.y(), up.z());
}

template <typename Real>
auto Camera<Real>::LookAt(Real eyex, Real eyey, Real eyez, Real atx, Real aty,
                          Real atz, Real upx, Real upy, Real upz) -> Matrix4r {
    Matrix4r matrix;
    Vector3r x, y, z;
    Vector3r eye(eyex, eyey, eyez);

    y = Vector3r(upx, upy, upz);
    z = Vector3r(atx - eyex, aty - eyey, atz - eyez);
    x = y.cross(z).normalized();
    y = z.cross(x).normalized();
    z.normalize();

    matrix(0, 0) = -x.x();
    matrix(0, 1) = -x.y();
    matrix(0, 2) = -x.z();
    matrix(0, 3) = x.dot(eye);

    matrix(1, 0) = y.x();
    matrix(1, 1) = y.y();
    matrix(1, 2) = y.z();
    matrix(1, 3) = -y.dot(eye);

    matrix(2, 0) = -z.x();
    matrix(2, 1) = -z.y();
    matrix(2, 2) = -z.z();
    matrix(2, 3) = z.dot(eye);

    matrix(3, 0) = Real(0);
    matrix(3, 1) = Real(0);
    matrix(3, 2) = Real(0);
    matrix(3, 3) = Real(1);

    return matrix;
}

template <typename Real>
auto Camera<Real>::PerspectiveMatrix(Real fov, std::size_t width,
                                     std::size_t height, Real nearPlane,
                                     Real farPlane) noexcept -> Matrix4r {
    Matrix4r proj;
    proj.setZero();

    Real aspectRatio = static_cast<Real>(width) / static_cast<Real>(height);
    Real ymax = nearPlane * std::tan(fov * Real(kPi) / Real(360));
    Real xmax = ymax * aspectRatio;
    Real left = -xmax;
    Real right = xmax;
    Real top = ymax;
    Real bottom = -ymax;

    Real temp, temp2, temp3, temp4;
    temp = Real(2) * nearPlane;
    temp2 = right - left;
    temp3 = top - bottom;
    temp4 = farPlane - nearPlane;

    proj(0) = temp / temp2;
    proj(5) = temp / temp3;
    proj(8) = (right + left) / temp2;
    proj(9) = (top + bottom) / temp3;
    proj(10) = (-farPlane - nearPlane) / temp4;
    proj(11) = Real(-1);
    proj(14) = (-temp * farPlane) / temp4;

    return proj;
}

template <typename Real>
auto Camera<Real>::OrthographicMatrix(Real left, Real right, Real bottom,
                                      Real top, Real nearPlane,
                                      Real farPlane) noexcept -> Matrix4r {
    Matrix4r proj;

    Real x = Real(2) / (right - left);
    Real y = Real(2) / (top - bottom);
    Real z = Real(-2) / (farPlane - nearPlane);
    Real tx = (right + left) / (right - left);
    Real ty = (top + bottom) / (top - bottom);
    Real tz = (farPlane + nearPlane) / (farPlane - nearPlane);

    proj.setZero();

    proj(0, 0) = x;
    proj(1, 1) = y;
    proj(2, 2) = z;
    proj(3, 3) = Real(1);

    proj(0, 3) = -tx;
    proj(1, 3) = -ty;
    proj(2, 3) = -tz;

    return proj;
}

/*
 * Cartesian from Spherical (form)
 * x = r * cos(theta) * sin(phi)
 * y = r * sin(theta) * sin(phi)
 * z = r * cos(phi)
 *
 * Spherical from Cartesian (form)
 * r = sqrt(x^2 + y^2 + z^2)
 * theta = atan(y / x)
 * phi = acos(z / r);
 *
 * Unswapped: R(r, t, p) = rcos(theta)sin(phi)i + rsin(phi)sin(theta)j +
 * rcos(phi)k
 */
template <typename Real>
inline auto Camera<Real>::SphereicalToCartesian(Real r, Real theta, Real phi)
    -> Vector3r {
    Vector3r result;

    Real sinPhi = std::sin(phi);
    Real cosPhi = std::cos(phi);
    Real sinTheta = std::sin(theta);
    Real cosTheta = std::cos(theta);

    result.x() = r * (cosTheta * sinPhi);
    result.z() = r * (sinTheta * sinPhi);
    result.y() = r * cosPhi;
    return result;
}

/* Unswapped: Rt(r, t, p) = -rsin(phi)sin(theta)i + rsin(phi)cos(theta)j + 0k */
template <typename Real>
inline auto Camera<Real>::SphereicalToCartesian_dTheta(Real r, Real theta,
                                                       Real phi) -> Vector3r {
    Vector3r result;

    Real sinPhi = std::sin(phi);
    Real cosPhi = std::cos(phi);
    Real sinTheta = std::sin(theta);
    Real cosTheta = std::cos(theta);

    result.x() = -r * (sinPhi * sinTheta);
    result.z() = r * (sinPhi * cosTheta);
    result.y() = Real(0);
    return result;
}

/* Unswapped: Rp(r, t, p) = rcos(phi)cos(theta)i + rcos(phi)sin(theta)j -
 * rsin(phi)k */
template <typename Real>
inline auto Camera<Real>::SphereicalToCartesian_dPhi(Real r, Real theta,
                                                     Real phi) -> Vector3r {
    Vector3r result;

    Real sinPhi = std::sin(phi);
    Real cosPhi = std::cos(phi);
    Real sinTheta = std::sin(theta);
    Real cosTheta = std::cos(theta);

    result.x() = r * (cosPhi * cosTheta);
    result.z() = r * (cosPhi * sinTheta);
    result.y() = -r * sinPhi;
    return result;
}

/* Rp X Rt = r^2 * sin^2(phi)cos(theta)i + r^2 * sin^2(phi)sin(theta)j + r^2 *
 * sin(phi)cos(phi)k */
template <typename Real>
inline auto Camera<Real>::SphereicalToCartesian_dPhiCrossdTheta(Real r,
                                                                Real theta,
                                                                Real phi)
    -> Vector3r {
    Vector3r result;

    Real rs = (r * r);
    Real sinPhi = std::sin(phi);
    Real cosPhi = std::cos(phi);
    Real sinTheta = std::sin(theta);
    Real cosTheta = std::cos(theta);

    result.x() = -rs * ((sinPhi * sinPhi) * cosTheta);
    result.y() = -rs * ((sinPhi * sinPhi) * sinTheta);
    result.z() = -rs * sinPhi * cosPhi;
    return result;
}

template <typename Real> inline void Camera<Real>::Compile() noexcept {
    look_at_ = Vector3r::Zero();
    eye_ = Camera<Real>::SphereicalToCartesian(r_, theta_, phi_);
    up_ =
        Camera<Real>::SphereicalToCartesian_dPhi(r_, theta_, phi_).normalized();

    //--------------------------------------------------------------------------------
    // Invert the up direction (since the spherical coordinates have phi
    // increasing downwards. Therefore we would like to have the (vector)
    // direction of the derivative inversed.
    //--------------------------------------------------------------------------------
    up_ *= Real(-1.0);

    look_at_ += displacement_;
    eye_ += displacement_;

    view_matrix_ = Camera<Real>::LookAt(eye_, look_at_, up_);
}
