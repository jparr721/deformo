#pragma once
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Geometry>

typedef float Real;

template <typename Real> using Vector2 = Eigen::Matrix<Real, 2, 1>;

template <typename Real> using Vector3 = Eigen::Matrix<Real, 3, 1>;

template <typename Real> using Matrix4 = Eigen::Matrix<Real, 4, 4>;

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

    auto rotate(Real du, Real dv) -> bool;
    auto pan(Real du, Real dv) -> bool;
    auto zoom(Real du) -> bool;

    void setPerspective(Real fov, std::size_t width, std::size_t height,
                        Real nearPlane, Real farPlane);
    void setFrustum(Real left, Real right, Real top, Real bottom,
                    Real nearPlane, Real farPlane);

    void setPosition(const Vector3<Real>& position);
    void setPosition(Real x, Real y, Real z);
    void addPosition(const Vector3<Real>& displacement);
    void addPosition(Real dx, Real dy, Real dz);

    void setRotation(Real theta, Real phi);
    void setHorizontalRotation(Real theta);
    void setVerticalRotation(Real phi);
    void setSphericalPosition(Real r, Real theta, Real phi);
    void addRadius(Real dRadius);
    void setRadius(Real radius);
    void setFieldOfView(Real fov);
    void setNearPlane(Real nearPlane);
    void setFarPlane(Real farPlane);
    void setMinRadius(Real minRadius);
    void setMaxRadius(Real maxRadius);

    void setRotationSensitivity(Real sen);
    void setPanSensitivity(Real sen);
    void setZoomSensitivity(Real sen);
    void setScopeSensitivity(Real sen);

    void reset();
    void resetPlanes();
    void resetView();
    void resetMatrices();
    void resetSensitivities();

    auto getRotationSensitivity() const -> Real;
    auto getPanSensitivity() const -> Real;
    auto getZoomSensitivity() const -> Real;
    auto getScopeSensitivity() const -> Real;

    auto getViewDirection() -> Vector3<Real>;
    auto getRightDirection() -> Vector3<Real>;
    auto getLeftDirection() -> Vector3<Real>;
    auto getUpDirection() -> Vector3<Real>;
    auto getDownDirection() -> Vector3<Real>;

    auto getEye() const -> const Vector3<Real>&;
    auto getLookAt() const -> const Vector3<Real>&;
    auto getUp() const -> const Vector3<Real>&;

    auto getWidth() const -> std::size_t;
    auto getHeight() const -> std::size_t;
    auto getNear() const -> const Real&;
    auto getFar() const -> const Real&;

    auto getViewMatrix() -> const Matrix4<Real>&;
    auto getProjectionMatrix() -> const Matrix4<Real>&;

    auto toViewMatrix() -> Matrix4<Real>;
    auto toProjectionMatrix() -> Matrix4<Real>;

    static auto LookAt(const Vector3<Real>& eye, const Vector3<Real>& lookAt,
                       const Vector3<Real>& up) -> Matrix4<Real>;
    static auto LookAt(Real eyex, Real eyey, Real eyez, Real atx, Real aty,
                       Real atz, Real upx, Real upy, Real upz) -> Matrix4<Real>;

    static auto PerspectiveMatrix(Real fov, std::size_t width,
                                  std::size_t height, Real nearPlane,
                                  Real farPlane) noexcept -> Matrix4<Real>;
    static auto OrthographicMatrix(Real left, Real right, Real bottom, Real top,
                                   Real nearPlane, Real farPlane) noexcept
        -> Matrix4<Real>;

    static auto SphereicalToCartesian(Real r, Real theta, Real phi)
        -> Vector3<Real>;
    static auto SphereicalToCartesian_dTheta(Real r, Real theta, Real phi)
        -> Vector3<Real>;
    static auto SphereicalToCartesian_dPhi(Real r, Real theta, Real phi)
        -> Vector3<Real>;
    static auto SphereicalToCartesian_dPhiCrossdTheta(Real r, Real theta,
                                                      Real phi)
        -> Vector3<Real>;

  protected:
    void compile() noexcept;

  protected:
    Matrix4<Real> view_matrix_;
    Matrix4<Real> projection_matrix_;
    Real near_plane_, far_plane_;

    std::size_t width_, height_;
    Vector3<Real> eye_;
    Vector3<Real> look_at_;
    Vector3<Real> up_;

    Real fov_, aspect_ratio_;
    Real min_radius_, max_radius_;

    Real r_, theta_, phi_;
    Vector3<Real> displacement_;

    Real pan_sensitivity_;
    Real zoom_sensitivity_;
    Real rotation_sensitivity_;
    Real scope_sensitivity_;
};

constexpr Real kPi = static_cast<Real>(3.14159);
constexpr Real kCameraEpsilon = static_cast<Real>(0.00001);
constexpr Real kDefaultNearPlane = static_cast<Real>(-1000.0);
constexpr Real kDefaultFarPlane = static_cast<Real>(10000.0);
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
    view_matrix_ = Matrix4<Real>::Identity();
    projection_matrix_ = Matrix4<Real>::Identity();

    eye_ = Vector3<Real>::UnitZ();
    look_at_ = Vector3<Real>::Zero();
    up_ = Vector3<Real>::UnitY();

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
    displacement_ = Vector3<Real>::Zero();
}

template <typename Real> Camera<Real>::~Camera() {}

template <typename Real>
auto Camera<Real>::Resize(std::size_t width, std::size_t height) -> bool {
    std::cout << "Near: " << near_plane_ << " Far: " << far_plane_ << std::endl;
    setPerspective(fov_, width, height, near_plane_, far_plane_);
    compile();
    return true;
}

template <typename Real>
auto Camera<Real>::Resize(std::size_t width, std::size_t height,
                          Real near_plane, Real far_plane) -> bool {
    setPerspective(fov_, width, height, near_plane, far_plane);
    compile();
    return true;
}

template <typename Real> auto Camera<Real>::rotate(Real du, Real dv) -> bool {
    theta_ -= (du * rotation_sensitivity_);
    phi_ += (dv * rotation_sensitivity_);
    compile();
    return true;
}

template <typename Real> auto Camera<Real>::pan(Real du, Real dv) -> bool {
    Vector3<Real> uDir = getLeftDirection();
    Vector3<Real> vDir = getDownDirection();

    Vector3<Real> uDisp = (du * pan_sensitivity_) * uDir;
    Vector3<Real> vDisp = (dv * pan_sensitivity_) * vDir;
    Vector3<Real> panDisp = uDisp + vDisp;

    displacement_ += panDisp;
    compile();
    return true;
}

template <typename Real> auto Camera<Real>::zoom(Real dr) -> bool {
    if (r_ + dr > max_radius_) {
        r_ = max_radius_;
        compile();
        return false;
    }
    if (r_ + dr < min_radius_) {
        r_ = min_radius_;
        compile();
        return false;
    }
    r_ += dr;
    compile();
    return true;
}

template <typename Real>
void Camera<Real>::setPerspective(Real fov, std::size_t width,
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
    setFrustum(-xmax, xmax, ymax, -ymax, nearPlane, farPlane);
}

template <typename Real>
void Camera<Real>::setFrustum(Real left, Real right, Real top, Real bottom,
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

template <typename Real>
void Camera<Real>::setPosition(const Vector3<Real>& pos) {
    position.x() = pos.x();
    position.y() = pos.y();
    position.z() = pos.z();
    compile();
}

template <typename Real>
void Camera<Real>::setPosition(Real x, Real y, Real z) {
    position.x() = x;
    position.y() = y;
    position.z() = z;
    compile();
}

template <typename Real>
void Camera<Real>::addPosition(const Vector3<Real>& displacement) {
    position += displacement;
    compile();
}

template <typename Real>
void Camera<Real>::addPosition(Real dx, Real dy, Real dz) {
    displacement_.x() += dx;
    displacement_.y() += dy;
    displacement_.z() += dz;
    compile();
}

template <typename Real> void Camera<Real>::setRotation(Real theta, Real phi) {
    theta_ = theta;
    phi_ = phi;
    compile();
}

template <typename Real> void Camera<Real>::setHorizontalRotation(Real theta) {
    theta_ = theta;
    compile();
}

template <typename Real> void Camera<Real>::setVerticalRotation(Real phi) {
    phi_ = phi;
    compile();
}

template <typename Real>
void Camera<Real>::setSphericalPosition(Real r, Real theta, Real phi) {
    r_ = r;
    theta_ = theta;
    phi_ = phi;
    compile();
}

template <typename Real> void Camera<Real>::addRadius(Real dRadius) {
    r_ += dRadius;
    compile();
}

template <typename Real> void Camera<Real>::setRadius(Real radius) {
    r_ = radius;
    compile();
}

template <typename Real> void Camera<Real>::setFieldOfView(Real fov) {
    if (fov > Real(kMaxFov))
        fov_ = Real(kMaxFov);
    else if (fov < Real(kMinFov))
        fov_ = Real(kMinFov);
    else
        fov_ = fov;
    setPerspective(fov_, aspect_ratio_, near_plane_, far_plane_);
}

template <typename Real> void Camera<Real>::setNearPlane(Real nearPlane) {
    near_plane_ = nearPlane;
    setPerspective(fov_, aspect_ratio_, nearPlane, far_plane_);
}

template <typename Real> void Camera<Real>::setFarPlane(Real farPlane) {
    far_plane_ = farPlane;
    setPerspective(fov_, aspect_ratio_, near_plane_, farPlane);
}

template <typename Real> void Camera<Real>::setMinRadius(Real minRadius) {
    min_radius_ = minRadius;
}

template <typename Real> void Camera<Real>::setMaxRadius(Real maxRadius) {
    max_radius_ = maxRadius;
}

template <typename Real> void Camera<Real>::setRotationSensitivity(Real sen) {
    rotation_sensitivity_ = sen;
}

template <typename Real> void Camera<Real>::setPanSensitivity(Real sen) {
    pan_sensitivity_ = sen;
}

template <typename Real> void Camera<Real>::setZoomSensitivity(Real sen) {
    zoom_sensitivity_ = sen;
}

template <typename Real> void Camera<Real>::setScopeSensitivity(Real sen) {
    scope_sensitivity_ = sen;
}

template <typename Real> void Camera<Real>::reset() {
    resetMatrices();
    resetPlanes();
    resetView();
    resetSensitivities();
    compile();
}

template <typename Real> void Camera<Real>::resetPlanes() {
    near_plane_ = Real(kDefaultNearPlane);
    far_plane_ = Real(kDefaultFarPlane);
}

template <typename Real> void Camera<Real>::resetView() {
    eye_ = Vector3<Real>::UnitZ();
    look_at_ = Vector3<Real>::Zero();
    up_ = Vector3<Real>::UnitY();
}

template <typename Real> void Camera<Real>::resetMatrices() {
    view_matrix_ = Matrix4<Real>::Identity();
    projection_matrix_ = Matrix4<Real>::Identity();
}

template <typename Real> void Camera<Real>::resetSensitivities() {
    rotation_sensitivity_ = Real(kDefaultRotationSensitivity);
    pan_sensitivity_ = Real(kDefaultPanSensitivity);
    zoom_sensitivity_ = Real(kDefaultZoomSensitivity);
    scope_sensitivity_ = Real(kDefaultScopeSensitivity);
}

template <typename Real>
auto Camera<Real>::getRotationSensitivity() const -> Real {
    return rotation_sensitivity_;
}

template <typename Real> auto Camera<Real>::getPanSensitivity() const -> Real {
    return pan_sensitivity_;
}

template <typename Real> auto Camera<Real>::getZoomSensitivity() const -> Real {
    return zoom_sensitivity_;
}

template <typename Real>
auto Camera<Real>::getScopeSensitivity() const -> Real {
    return scope_sensitivity_;
}

template <typename Real>
inline auto Camera<Real>::getViewDirection() -> Vector3<Real> {
    compile();
    return (look_at_ - eye_).normalized();
}

template <typename Real>
inline auto Camera<Real>::getRightDirection() -> Vector3<Real> {
    compile();
    Vector3<Real> dir = (look_at_ - eye_).normalized();
    return (up_.cross(dir)).normalized();
}

template <typename Real>
inline auto Camera<Real>::getLeftDirection() -> Vector3<Real> {
    compile();
    Vector3<Real> dir = (look_at_ - eye_).normalized();
    return (dir.cross(up_)).normalized();
}

template <typename Real>
inline auto Camera<Real>::getUpDirection() -> Vector3<Real> {
    compile();
    return up_;
}

template <typename Real>
inline auto Camera<Real>::getDownDirection() -> Vector3<Real> {
    compile();
    return -up_;
}

template <typename Real>
auto Camera<Real>::getEye() const -> const Vector3<Real>& {
    return eye_;
}

template <typename Real>
auto Camera<Real>::getLookAt() const -> const Vector3<Real>& {
    return look_at_;
}

template <typename Real>
auto Camera<Real>::getUp() const -> const Vector3<Real>& {
    return up_;
}

template <typename Real> auto Camera<Real>::getWidth() const -> std::size_t {
    return width_;
}

template <typename Real> auto Camera<Real>::getHeight() const -> std::size_t {
    return height_;
}

template <typename Real> auto Camera<Real>::getNear() const -> const Real& {
    return near_plane_;
}

template <typename Real> auto Camera<Real>::getFar() const -> const Real& {
    return far_plane_;
}

template <typename Real>
auto Camera<Real>::getViewMatrix() -> const Matrix4<Real>& {
    compile();
    return view_matrix_;
}

template <typename Real>
auto Camera<Real>::getProjectionMatrix() -> const Matrix4<Real>& {
    compile();
    return projection_matrix_;
}

template <typename Real> auto Camera<Real>::toViewMatrix() -> Matrix4<Real> {
    compile();
    return view_matrix_;
}

template <typename Real>
auto Camera<Real>::toProjectionMatrix() -> Matrix4<Real> {
    compile();
    return projection_matrix_;
}

template <typename Real>
auto Camera<Real>::LookAt(const Vector3<Real>& eye, const Vector3<Real>& lookAt,
                          const Vector3<Real>& up) -> Matrix4<Real> {
    return Camera<Real>::LookAt(eye.x(), eye.y(), eye.z(), lookAt.x(),
                                lookAt.y(), lookAt.z(), up.x(), up.y(), up.z());
}

template <typename Real>
auto Camera<Real>::LookAt(Real eyex, Real eyey, Real eyez, Real atx, Real aty,
                          Real atz, Real upx, Real upy, Real upz)
    -> Matrix4<Real> {
    Matrix4<Real> matrix;
    Vector3<Real> x, y, z;
    Vector3<Real> eye(eyex, eyey, eyez);

    y = Vector3<Real>(upx, upy, upz);
    z = Vector3<Real>(atx - eyex, aty - eyey, atz - eyez);
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
                                     Real farPlane) noexcept -> Matrix4<Real> {
    Matrix4<Real> proj;
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
                                      Real farPlane) noexcept -> Matrix4<Real> {
    Matrix4<Real> proj;

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
    -> Vector3<Real> {
    Vector3<Real> result;

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
                                                       Real phi)
    -> Vector3<Real> {
    Vector3<Real> result;

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
                                                     Real phi)
    -> Vector3<Real> {
    Vector3<Real> result;

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
    -> Vector3<Real> {
    Vector3<Real> result;

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

template <typename Real> inline void Camera<Real>::compile() noexcept {
    look_at_ = Vector3<Real>::Zero();
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
