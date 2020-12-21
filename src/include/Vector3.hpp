#pragma once
#include <cmath>

namespace mathUtl {
    class Vec3 {
    public:
        double x;
        double y;
        double z;

        Vec3() = default;

        constexpr Vec3(const Vec3& v_) : x(v_.x), y(v_.y), z(v_.z) {}

        constexpr Vec3(const double x_, const double y_, const double z_) : x(x_), y(y_), z(z_) {}

        constexpr double Elem(const size_t index_) const {
            return index_ == 0 ? x : index_ == 1 ? y : index_ == 2 ? z : 0;
        }

        constexpr Vec3 operator+() const {
            return *this;
        }

        constexpr Vec3 operator-() const {
            return {-x, -y, -z};
        }

        constexpr Vec3 operator+(const Vec3& v_) const {
            return {x + v_.x, y + v_.y, z + v_.z};
        }

        constexpr Vec3 operator-(const Vec3& v_) const {
            return {x - v_.x, y - v_.y, z - v_.z};
        }

        constexpr Vec3 operator*(const double s_) const {
            return {x * s_, y * s_, z * s_};
        }

        constexpr Vec3 operator*(const Vec3& v_) const {
            return {x * v_.x, y * v_.y, z * v_.z};
        }

        constexpr Vec3 operator/(const double s_) const {
            return *this * (1.0 / s_);
        }

        constexpr Vec3 operator/(const Vec3& v_) const {
            return {x / v_.x, y / v_.y, z / v_.z};
        }

        constexpr Vec3& operator+=(const Vec3& v_) {
            x += v_.x;
            y += v_.y;
            z += v_.z;
            return *this;
        }

        constexpr Vec3& operator-=(const Vec3& v_) {
            x -= v_.x;
            y -= v_.y;
            z -= v_.z;
            return *this;
        }

        constexpr Vec3& operator*=(const double s_) {
            x *= s_;
            y *= s_;
            z *= s_;
            return *this;
        }

        constexpr Vec3& operator*=(const Vec3& v_) {
            x *= v_.x;
            y *= v_.y;
            z *= v_.z;
            return *this;
        }

        constexpr Vec3& operator/=(const double s_) {
            return *this *= (1.0 / s_);
        }

        constexpr Vec3& operator/=(const Vec3& v_) {
            x /= v_.x;
            y /= v_.y;
            z /= v_.z;
            return *this;
        }

        constexpr bool operator==(const Vec3& v_) const {
            return v_.x == x && v_.y == y && v_.z == z;
        }

        constexpr bool operator!=(const Vec3& v_) const {
            return v_.x != x || v_.y != y || v_.z != z;
        }

        constexpr Vec3& Set(const double x_, const double y_, const double z_) {
            x = x_;
            y = y_;
            z = z_;
            return *this;
        }

        constexpr Vec3& Set(const Vec3& v_) {
            return *this = v_;
        }

        constexpr Vec3 MovedBy(const double x_, const double y_, const double z_) const {
            return {x + x_, y + y_, z + z_};
        }

        constexpr Vec3 MovedBy(const Vec3& v_) const {
            return {x + v_.x, y + v_.y, z + v_.z};
        }

        constexpr Vec3& MoveBy(const double x_, const double y_, const double z_) {
            x += x_;
            y += y_;
            z += z_;
            return *this;
        }

        constexpr Vec3& MoveBy(const Vec3& v_) {
            return *this += v_;
        }

        constexpr bool IsZero() const {
            return x == 0.0 && y == 0.0 && z == 0.0;
        }

        bool HasNaN() const {
            return std::isnan(x) || std::isnan(y) || std::isnan(z);
        }

        constexpr double Dot(const Vec3& v_) const {
            return x * v_.x + y * v_.y + z * v_.z;
        }

        constexpr Vec3 Cross(const Vec3& v_) const {
            return {y * v_.z - z * v_.y, z * v_.x - x * v_.z, x * v_.y - y * v_.x};
        }

        double Length() const {
            return std::sqrt(LengthSq());
        }

        constexpr double LengthSq() const {
            return Dot(*this);
        }

        double LengthInv() const {
            return 1.0 / Length();
        }

        Vec3& SetLength(const double length_) {
            const double len = Length();

            if (len == 0.0) {
                return *this;
            }

            return *this *= (length_ / len);
        }

        Vec3 ClampLength(const double maxLength_) const {
            const double len = Length();

            if (len <= maxLength_) {
                return *this;
            }
            else {
                return *this * (maxLength_ / len);
            }
        }

        double DistanceFrom(const double x_, const double y_, const double z_) const {
            return (*this - Vec3(x_, y_, z_)).Length();
        }

        double DistanceFrom(const Vec3& v_) const {
            return (*this - v_).Length();
        }

        constexpr double DistanceFromSq(const double x_, const double y_, const double z_) const {
            return (*this - Vec3(x_, y_, z_)).LengthSq();
        }

        constexpr double DistanceFromSq(const Vec3& v_) const {
            return (*this - v_).LengthSq();
        }

        Vec3 Normalized() const {
            return *this * LengthInv();
        }

        Vec3& Normalize() {
            return *this *= LengthInv();
        }

        constexpr Vec3 Lerp(const Vec3& other_, double f_) const {
            return Vec3(x + (other_.x - x) * f_, y + (other_.y - y) * f_, z + (other_.z - z) * f_);
        }

        /// <summary>
        /// Vec3{ x, x, x }
        /// </summary>
        constexpr Vec3 XXX() const {
            return {x, x, x};
        }

        /// <summary>
        /// Vec3{ y, y, y }
        /// </summary>
        constexpr Vec3 YYY() const {
            return {y, y, y};
        }

        /// <summary>
        /// Vec3{ z, z, z }
        /// </summary>
        constexpr Vec3 ZZZ() const {
            return {z, z, z};
        }

        /// <summary>
        /// Vec3{ x, y, z }
        /// </summary>
        constexpr Vec3 XYZ() const {
            return {x, y, z};
        }

        /// <summary>
        /// Vec3{ x, z, y }
        /// </summary>
        constexpr Vec3 XZY() const {
            return {x, z, y};
        }

        /// <summary>
        /// Vec3{ y, x, z }
        /// </summary>
        constexpr Vec3 YXZ() const {
            return {y, x, z};
        }

        /// <summary>
        /// Vec3{ y, z, x }
        /// </summary>
        constexpr Vec3 YZX() const {
            return {y, z, x};
        }

        /// <summary>
        /// Vec3{ z, x, y }
        /// </summary>
        constexpr Vec3 ZXY() const {
            return {z, x, y};
        }

        /// <summary>
        /// Vec3{ z, y, x }
        /// </summary>
        constexpr Vec3 ZYX() const {
            return {z, y, x};
        }

        /// <summary>
        /// Vec3{ 0, 0, 0 }
        /// </summary>
        static constexpr Vec3 Zero() {
            return {0, 0, 0};
        }

        /// <summary>
        /// Vec3{ 1, 1, 1 }
        /// </summary>
        static constexpr Vec3 One() {
            return {1, 1, 1};
        }

        /// <summary>
        /// Vec3{ value_, value_, value_ }
        /// </summary>
        static constexpr Vec3 All(const double value_ = 1) {
            return {value_, value_, value_};
        }

        /// <summary>
        /// Vec3{ 1, 0, 0 }
        /// </summary>
        static constexpr Vec3 UnitX() {
            return {1, 0, 0};
        }

        /// <summary>
        /// Vec3{ 0, 1, 0 }
        /// </summary>
        static constexpr Vec3 UnitY() {
            return {0, 1, 0};
        }

        /// <summary>
        /// Vec3{ 0, 0, 1 }
        /// </summary>
        static constexpr Vec3 UnitZ() {
            return {0, 0, 1};
        }

        /// <summary>
        /// Vec3{ -length_, 0, 0 }
        /// </summary>
        static constexpr Vec3 Left(const double length_ = 1) {
            return {-length_, 0, 0};
        }

        /// <summary>
        /// Vec3{ length_, 0, 0 }
        /// </summary>
        static constexpr Vec3 Right(const double length_ = 1) {
            return {length_, 0, 0};
        }

        /// <summary>
        /// Vec3{ 0, length_, 0 }
        /// </summary>
        static constexpr Vec3 Up(const double length_ = 1) {
            return {0, length_, 0};
        }

        /// <summary>
        /// Vec3{ 0, -length_, 0 }
        /// </summary>
        static constexpr Vec3 Down(const double length_ = 1) {
            return {0, -length_, 0};
        }

        /// <summary>
        /// Vec3{ 0, 0, length_ }
        /// </summary>
        static constexpr Vec3 Forward(const double length_ = 1) {
            return {0, 0, length_};
        }

        /// <summary>
        /// Vec3{ 0, 0, -length_ }
        /// </summary>
        static constexpr Vec3 Backward(const double length_ = 1) {
            return {0, 0, -length_};
        }
    };

    inline constexpr Vec3 operator*(const double s_, const Vec3& v_) {
        return v_ * s_;
    }
}  // namespace mathUtl