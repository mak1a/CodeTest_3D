#pragma once
#include "MathUtil.hpp"
#include "Vector3.hpp"
#include <cassert>

namespace mathUtl {
    class Quaternion;
    class Mat4x3;
    class RotationMatrix;

    /// <summary>
    /// オイラー角3つ組で表す
    /// </summary>
    class EulerAngles {
    public:
        /// <summary>
        /// 3つの角度をラジアンで保存する
        /// </summary>
        double heading;
        double pitch;
        double bank;

        EulerAngles() = default;

        constexpr EulerAngles(const double heading_, const double pitch_, const double bank_) : heading(heading_), pitch(pitch_), bank(bank_) {}

        /// <summary>
        /// 恒等化する
        /// 全て0を設定する
        /// </summary>
        void Identity() {
            heading = 0.0;
            pitch = 0.0;
            bank = 0.0;
        }

        /// <summary>
        /// 正準オイラー角の3組を決定する
        /// </summary>
        void Canonize() {
            // 最初に、ピッチを範囲 -π...πにラップする
            pitch = WrapPi(pitch);

            // ここで、正準範囲 -π/2...π/2の外側で行列pitchの裏側をチェックする
            if (pitch < -PiOverTwo) {
                pitch = -Pi - pitch;
                heading += Pi;
                bank += Pi;
            }
            else if (pitch > PiOverTwo) {
                pitch = Pi - pitch;
                heading += Pi;
                bank += Pi;
            }

            // ここでジンバルロックのケースをチェックする
            if (std::fabs(pitch) > PiOverTwo - 1e-4) {
                // ジンバルロック内にいる
                // 垂直軸に関する全ての回転をヘディングに割り当てる
                heading += bank;
                bank = 0.0;
            }
            else {
                // ジンバルロック外
                // バンク角を正準備範囲にラップする
                bank = WrapPi(bank);
            }

            // ヘディングを正準範囲にラップする
            heading = WrapPi(heading);
        }

        /// <summary>
        /// 四元数をオイラー角形式に変換する
        /// 入力される四元数は、その名前が示すようにオブジェクト空間から慣性空間、
        /// または、慣性空間からオブジェクト空間への回転を実行するものとする
        /// </summary>
        void FromObjectToInertialQuaternion(const Quaternion& q_);
        void FromInertialToObjectQuaternion(const Quaternion& q_);

        /// <summary>
        /// 座標変換行列をオイラー角形式に変換する
        /// 入力される行列は、それが示すようにオブジェクト空間からワールド空間、
        /// ワールド空間からオブジェクト空間への座標変換を実行するものとする
        /// この行列は、直交しているものとする
        /// </summary>
        void FromObjectToWorldMatrix(const Mat4x3& m_);
        void FromWorldToObjectMatrix(const Mat4x3& m_);

        /// <summary>
        /// 回転行列をオイラー角形式に変換する
        /// </summary>
        void FromRotationMatrix(const RotationMatrix& m_);
    };

    /// <summary>
    /// 3Dにおける各変位(方向)を表す為の四元数を実装する
    /// </summary>
    class Quaternion {
    public:
        /// <summary>
        /// 四元数の4つの値
        /// </summary>
        double w;
        double x;
        double y;
        double z;

        Quaternion() = default;

        constexpr Quaternion(const double w_, const double x_, const double y_, const double z_) : w(w_), x(x_), y(y_), z(z_) {}

        /// <summary>
        /// 恒等四元数にする
        /// </summary>
        void Identity() {
            w = 1.0;
            x = 0.0;
            y = 0.0;
            z = 0.0;
        }

        /// <summary>
        /// 四元数を特定の角度にセットアップする
        /// </summary>
        void SetToRotateAboutX(const double theta_) {
            // 半分の角度を計算
            const double thetaOverTwo{theta_ * 0.5};

            w = std::cos(thetaOverTwo);
            x = std::sin(thetaOverTwo);
            y = 0.0;
            z = 0.0;
        }
        void SetToRotateAboutY(const double theta_) {
            const double thetaOverTwo{theta_ * 0.5};

            w = std::cos(thetaOverTwo);
            x = 0.0;
            y = std::sin(thetaOverTwo);
            z = 0.0;
        }
        void SetToRotateAboutZ(const double theta_) {
            const double thetaOverTwo{theta_ * 0.5};

            w = std::cos(thetaOverTwo);
            x = 0.0;
            y = 0.0;
            z = std::sin(thetaOverTwo);
        }
        void SetToRotateAboutAxis(const Vec3& axis_, const double theta_) {
            assert(std::fabs(axis_.Length() - 1.0) < 0.01);

            const double thetaOverTwo{theta_ * 0.5};
            const double sinThetaOverTwo{std::sin(thetaOverTwo)};

            w = std::cos(thetaOverTwo);
            x = axis_.x * sinThetaOverTwo;
            y = axis_.y * sinThetaOverTwo;
            z = axis_.z * sinThetaOverTwo;
        }

        /// <summary>
        /// オブジェクト空間<->慣性空間の回転を実行するようにセットアップする
        /// 方向はオイラー角形式で与えられる
        /// </summary>
        void SetToRotateObjectToInertial(const EulerAngles& orientation_);
        void SetToRotateInertialToObject(const EulerAngles& orientation_);

        /// <summary>
        /// 外積
        /// </summary>
        constexpr Quaternion operator*(const Quaternion& q_) const {
            return {w * q_.w - x * q_.x - y * q_.y - z * q_.z, w * q_.x + x * q_.w + z * q_.y - y * q_.z, w * q_.y + y * q_.w + x * q_.z - z * q_.x,
                    w * q_.z + z * q_.w + y * q_.x - x * q_.y};
        }

        constexpr Quaternion& operator*=(const Quaternion& q_) {
            *this = *this * q_;
            return *this;
        }

        /// <summary>
        /// 四元数を正規化する
        /// </summary>
        void Normalize() {
            const double mag{std::sqrt(w * w + x * x + y * y + z * z)};

            if (mag > 0.0) {
                const double oneOverMag{1.0 / mag};
                w *= oneOverMag;
                x *= oneOverMag;
                y *= oneOverMag;
                z *= oneOverMag;
                return;
            }

            // 問題発生
            assert(false);

            // リリースビルドでは、何かに割り当てる
            Identity();
        }

        /// <summary>
        /// 回転角と軸を取り出し返す
        /// </summary>
        constexpr double GetRotationAngle() const {
            return SafeAcos(w) * 2.0;
        }
        constexpr Vec3 GetRotationAxis() const {
            const double sinThetaOverTwoSq{1.0 - w * w};

            if (sinThetaOverTwoSq <= 0.0) {
                return {1.0, 0.0, 0.0};
            }

            const double oneOverSinThetaOverTwo{1.0 / std::sqrt(sinThetaOverTwoSq)};
            return {x * oneOverSinThetaOverTwo, y * oneOverSinThetaOverTwo, z * oneOverSinThetaOverTwo};
        }

        constexpr double DotProduct(const Quaternion& q_) const {
            return w * q_.w + x * q_.x + y * q_.y + z * q_.z;
        }

        constexpr Quaternion Slerp(const Quaternion& q_, const double t_) const {
            if (t_ <= 0.0) {
                return *this;
            }
            if (t_ >= 1.0) {
                return q_;
            }

            double cosOmega{DotProduct(q_)};

            double qw{q_.w};
            double qx{q_.x};
            double qy{q_.y};
            double qz{q_.z};

            if (cosOmega < 0.0) {
                qw = -qw;
                qx = -qx;
                qy = -qy;
                qz = -qz;
                cosOmega = -cosOmega;
            }

            assert(cosOmega < 1.1);

            double k0{}, k1{};
            if (cosOmega > 0.9999) {
                k0 = 1.0 - t_;
                k1 = t_;
            }
            else {
                const double sinOmega{std::sqrt(1.0 - cosOmega * cosOmega)};

                // 角度を計算
                const double omega{std::atan2(sinOmega, cosOmega)};

                const double oneOverSinOmega{1.0 / sinOmega};

                k0 = std::sin((1.0 - t_) * omega) * oneOverSinOmega;
                k1 = std::sin(t_ * omega) * oneOverSinOmega;
            }

            return {k0 * w + k1 * qw, k0 * x + k1 * qx, k0 * y + k1 * qy, k0 * z + k1 * qz};
        }

        /// <summary>
        /// 四元数の共役
        /// 元の四元数の反対の回転を持つ四元数
        /// </summary>
        constexpr Quaternion Conjugate() const {
            return {w, -x, -y, -z};
        }

        /// <summary>
        /// 四元数の累乗
        /// </summary>
        Quaternion Pow(const double exponent_) const {
            if (std::fabs(w) > 0.9999) {
                return *this;
            }

            double alpha{std::acos(w)};
            double newAlpha{alpha * exponent_};
            double mult{std::sin(newAlpha) / std::sin(alpha)};

            return {std::cos(newAlpha), x * mult, y * mult, z * mult};
        }
    };

    /// <summary>
    /// 回転でのみ用いる3×3行列を実装する
    /// この行列は直交行列
    /// 座標変換の向きは、座標変換時に指定される
    /// </summary>
    class RotationMatrix {
    public:
        double m11, m12, m13;
        double m21, m22, m23;
        double m31, m32, m33;

        RotationMatrix() = default;
        constexpr RotationMatrix(const double m11_,
                                 const double m12_,
                                 const double m13_,
                                 const double m21_,
                                 const double m22_,
                                 const double m23_,
                                 const double m31_,
                                 const double m32_,
                                 const double m33_)
        : m11(m11_), m12(m12_), m13(m13_), m21(m21_), m22(m22_), m23(m23_), m31(m31_), m32(m32_), m33(m33_) {}

        void Identity() {
            m11 = 1.0;
            m12 = 0.0;
            m13 = 0.0;
            m21 = 0.0;
            m22 = 1.0;
            m23 = 0.0;
            m31 = 0.0;
            m32 = 0.0;
            m33 = 1.0;
        }

        void Setup(const EulerAngles& orientation_);

        void FromInertialToObjectQuaternion(const Quaternion& q_);
        void FromObjectToInertialQuaternion(const Quaternion& q_);

        constexpr Vec3 InertialToObject(const Vec3& v_) const {
            return {m11 * v_.x + m21 * v_.y + m31 * v_.z, m12 * v_.x + m22 * v_.y + m32 * v_.z, m13 * v_.x + m23 * v_.y + m33 * v_.z};
        }

        constexpr Vec3 ObjectToInertial(const Vec3& v_) const {
            return {m11 * v_.x + m12 * v_.y + m13 * v_.z, m21 * v_.x + m22 * v_.y + m23 * v_.z, m31 * v_.x + m32 * v_.y + m33 * v_.z};
        }
    };

    class Mat4x3 {
    public:
        double m11, m12, m13;
        double m21, m22, m23;
        double m31, m32, m33;
        // 平行移動部分
        double tx, ty, tz;

        Mat4x3() = default;

        constexpr Mat4x3(const double m11_,
                         const double m12_,
                         const double m13_,
                         const double m21_,
                         const double m22_,
                         const double m23_,
                         const double m31_,
                         const double m32_,
                         const double m33_,
                         const double tx_,
                         const double ty_,
                         const double tz_)
        : m11(m11_), m12(m12_), m13(m13_), m21(m21_), m22(m22_), m23(m23_), m31(m31_), m32(m32_), m33(m33_), tx(tx_), ty(ty_), tz(tz_) {}

        void Identity() {
            m11 = 1.0;
            m12 = 0.0;
            m13 = 0.0;
            m21 = 0.0;
            m22 = 1.0;
            m23 = 0.0;
            m31 = 0.0;
            m32 = 0.0;
            m33 = 1.0;
            tx = 0.0;
            ty = 0.0;
            tz = 0.0;
        }

        void ZeroTranslation() {
            tx = 0.0;
            ty = 0.0;
            tz = 0.0;
        }

        void SetTranslation(const Vec3& v_) {
            tx = v_.x;
            ty = v_.y;
            tz = v_.z;
        }

        void SetupTranslation(const Vec3& v_) {
            m11 = 1.0;
            m12 = 0.0;
            m13 = 0.0;
            m21 = 0.0;
            m22 = 1.0;
            m23 = 0.0;
            m31 = 0.0;
            m32 = 0.0;
            m33 = 1.0;
            tx = v_.x;
            ty = v_.y;
            tz = v_.z;
        }

        void SetupLocalToParent(const Vec3& pos_, const EulerAngles& orient_);
        void SetupLocalToParent(const Vec3& pos_, const RotationMatrix& orient_);
        void SetupParentToLocal(const Vec3& pos_, const EulerAngles& orient_);
        void SetupParentToLocal(const Vec3& pos_, const RotationMatrix& orient_);

        /// <summary>
        /// 基本軸の周りの回転を実行する行列をセットアップ
        /// </summary>
        /// <param name=axis_>1:x, 2:y, 3:z, それ以外:エラー</param>
        void SetupRotate(const int32 axis_, const float theta_) {
            double s{}, c{};
            SinCos(&s, &c, theta_);

            switch (axis_) {
            case 1:
                m11 = 1.0;
                m12 = 0.0;
                m13 = 0.0;
                m21 = 0.0;
                m22 = c;
                m23 = s;
                m31 = 0.0;
                m32 = -s;
                m33 = c;
                break;
            case 2:
                m11 = c;
                m12 = 0.0;
                m13 = -s;
                m21 = 0.0;
                m22 = 1.0;
                m23 = 0.0;
                m31 = s;
                m32 = 0.0;
                m33 = c;
                break;
            case 3:
                m11 = c;
                m12 = s;
                m13 = 0.0;
                m21 = -s;
                m22 = c;
                m23 = 0.0;
                m31 = 0.0;
                m32 = 0.0;
                m33 = 1.0;
                break;
            default:
                assert(false);
            }

            tx = 0.0;
            ty = 0.0;
            tz = 0.0;
        }

        void SetupRotate(const Vec3& axis_, const double theta_) {
            assert(std::fabs(axis_.Dot(axis_) - 1.0) < 0.01);

            double s{}, c{};
            SinCos(&s, &c, theta_);

            const double a{1.0 - c};
            const double ax{a * axis_.x};
            const double ay{a * axis_.y};
            const double az{a * axis_.z};

            m11 = ax * axis_.x + c;
            m12 = ax * axis_.y + axis_.z * s;
            m13 = ax * axis_.z - axis_.y * s;

            m21 = ay * axis_.x - axis_.z * s;
            m22 = ay * axis_.y + c;
            m23 = ay * axis_.z + axis_.x * s;

            m31 = az * axis_.x + axis_.y * s;
            m32 = az * axis_.y - axis_.x * s;
            m33 = az * axis_.z + c;

            tx = 0.0;
            ty = 0.0;
            tz = 0.0;
        }

        void FromQuaternion(const Quaternion& q_);

        void SetupScale(const Vec3& v_) {
            m11 = v_.x;
            m12 = 0.0;
            m13 = 0.0;
            m21 = 0.0;
            m22 = v_.y;
            m23 = 0.0;
            m31 = 0.0;
            m32 = 0.0;
            m33 = v_.z;

            tx = 0.0;
            ty = 0.0;
            tz = 0.0;
        }

        void SetupScaleAlongAxis(const Vec3& axis_, const float k_) {
            assert(std::fabs(axis_.Dot(axis_) - 1.0) < 0.01);

            const double a{k_ - 1.0};
            const double ax{a * axis_.x};
            const double ay{a * axis_.y};
            const double az{a * axis_.z};

            m11 = ax * axis_.x + 1.0;
            m22 = ay * axis_.y + 1.0;
            m33 = az * axis_.z + 1.0;

            m12 = ax * axis_.y;
            m21 = ax * axis_.y;
            m13 = ax * axis_.z;
            m31 = ax * axis_.z;
            m23 = ay * axis_.z;
            m32 = ay * axis_.z;

            tx = 0.0;
            ty = 0.0;
            tz = 0.0;
        }

        /// <summary>
        /// せん断を実行する行列をセットアップ
        /// </summary>
        /// <param name=axis_>1:x, 2:y, 3:z, それ以外:エラー</param>
        void SetupShear(const int32 axis_, const double s_, const double t_) {
            switch (axis_) {
            case 1:  // xを用いてyとzをせん断
                m11 = 1.0;
                m12 = s_;
                m13 = t_;
                m21 = 0.0;
                m22 = 1.0;
                m23 = 0.0;
                m31 = 0.0;
                m32 = 0.0;
                m33 = 1.0;
                break;
            case 2:  // yを用いてxとzをせん断
                m11 = 1.0;
                m12 = 0.0;
                m13 = 0.0;
                m21 = s_;
                m22 = 1.0;
                m23 = t_;
                m31 = 0.0;
                m32 = 0.0;
                m33 = 1.0;
                break;
            case 3:  // zを用いてxとyをせん断
                m11 = 1.0;
                m12 = 0.0;
                m13 = 0.0;
                m21 = 0.0;
                m22 = 1.0;
                m23 = 0.0;
                m31 = s_;
                m32 = t_;
                m33 = 1.0;
                break;
            default:
                assert(false);
            }

            tx = 0.0;
            ty = 0.0;
            tz = 0.0;
        }

        void SetupProject(const Vec3& n_) {
            assert(std::fabs(n_.Dot(n_) - 1.0) < 0.01);

            m11 = 1.0 - n_.x * n_.x;
            m22 = 1.0 - n_.y * n_.y;
            m33 = 1.0 - n_.z * n_.z;

            m12 = -n_.x * n_.y;
            m21 = -n_.x * n_.y;
            m13 = -n_.x * n_.z;
            m31 = -n_.x * n_.z;
            m23 = -n_.y * n_.z;
            m32 = -n_.y * n_.z;

            tx = 0.0;
            ty = 0.0;
            tz = 0.0;
        }

        void SetupReflect(const int32 axis_, const double k_) {
            switch (axis_) {
            case 1:  // 面x=kに関するリフレクション
                m11 = -1.0;
                m12 = 0.0;
                m13 = 0.0;
                m21 = 0.0;
                m22 = 1.0;
                m23 = 0.0;
                m31 = 0.0;
                m32 = 0.0;
                m33 = 1.0;

                tx = 2.0 * k_;
                ty = 0.0;
                tz = 0.0;
                break;
            case 2:  // 面x=kに関するリフレクション
                m11 = 1.0;
                m12 = 0.0;
                m13 = 0.0;
                m21 = 0.0;
                m22 = -1.0;
                m23 = 0.0;
                m31 = 0.0;
                m32 = 0.0;
                m33 = 1.0;

                tx = 0.0;
                ty = 2.0 * k_;
                tz = 0.0;
                break;
            case 3:  // 面x=kに関するリフレクション
                m11 = 1.0;
                m12 = 0.0;
                m13 = 0.0;
                m21 = 0.0;
                m22 = 1.0;
                m23 = 0.0;
                m31 = 0.0;
                m32 = 0.0;
                m33 = -1.0;

                tx = 0.0;
                ty = 0.0;
                tz = 2.0 * k_;
                break;
            default:
                assert(false);
            }
        }

        void SetupReflect(const Vec3& n_) {
            assert(std::fabs(n_.Dot(n_) - 1.0) < 0.01);

            const double ax{-2.0 * n_.x};
            const double ay{-2.0 * n_.y};
            const double az{-2.0 * n_.z};

            m11 = 1.0 + ax * n_.x;
            m22 = 1.0 + ay * n_.y;
            m33 = 1.0 + az * n_.z;

            m12 = ax * n_.y;
            m21 = ax * n_.y;
            m13 = ax * n_.z;
            m31 = ax * n_.z;
            m23 = ay * n_.z;
            m32 = ay * n_.z;

            tx = 0.0;
            ty = 0.0;
            tz = 0.0;
        }

        constexpr Mat4x3 operator*(const Mat4x3& m_) const {
            return {
                m11 * m_.m11 + m12 * m_.m21 + m13 * m_.m31,       // m11
                m11 * m_.m12 + m12 * m_.m22 + m13 * m_.m32,       // m12
                m11 * m_.m13 + m12 * m_.m23 + m13 * m_.m33,       // m13
                m21 * m_.m11 + m22 * m_.m21 + m23 * m_.m31,       // m21
                m21 * m_.m12 + m22 * m_.m22 + m23 * m_.m32,       // m22
                m21 * m_.m13 + m22 * m_.m23 + m23 * m_.m33,       // m23
                m31 * m_.m11 + m32 * m_.m21 + m33 * m_.m31,       // m31
                m31 * m_.m12 + m32 * m_.m22 + m33 * m_.m32,       // m32
                m31 * m_.m13 + m32 * m_.m23 + m33 * m_.m33,       // m33
                tx * m_.m11 + ty * m_.m21 + tz * m_.m31 + m_.tx,  // tx
                tx * m_.m12 + ty * m_.m22 + tz * m_.m32 + m_.ty,  // ty
                tx * m_.m13 + ty * m_.m23 + tz * m_.m33 + m_.tz   // tz
            };
        }

        constexpr Mat4x3& operator*=(const Mat4x3& m_) {
            *this = *this * m_;
            return *this;
        }

        constexpr double Determinant() const {
            return (m11 * (m22 * m33 - m23 * m32) + m12 * (m23 * m31 - m21 * m33) + m13 * (m21 * m32 - m22 * m31));
        }

        Mat4x3 Inverse() const {
            const double det{Determinant()};
            assert(std::fabs(det) > 0.000001);

            const double oneOverDet{1.0 / det};

            Mat4x3 r{
                (m22 * m33 - m23 * m32) * oneOverDet,  // m11
                (m13 * m32 - m12 * m33) * oneOverDet,  // m12
                (m12 * m23 - m13 * m22) * oneOverDet,  // m13
                (m23 * m31 - m21 * m33) * oneOverDet,  // m21
                (m11 * m33 - m13 * m31) * oneOverDet,  // m22
                (m13 * m21 - m11 * m23) * oneOverDet,  // m23
                (m21 * m32 - m22 * m31) * oneOverDet,  // m31
                (m12 * m31 - m11 * m32) * oneOverDet,  // m32
                (m11 * m22 - m12 * m21) * oneOverDet,  // m33
                0.0,                                   // tx
                0.0,                                   // ty
                0.0                                    // tz
            };

            r.tx = -(tx * r.m11 + ty * r.m21 + tz * r.m31);
            r.ty = -(tx * r.m12 + ty * r.m22 + tz * r.m32);
            r.tz = -(tx * r.m13 + ty * r.m23 + tz * r.m33);

            return r;
        }

        constexpr Vec3 GetTranslation() const {
            return {tx, ty, tz};
        }

        constexpr Vec3 GetPositionFromParentToLocalMatrix() const {
            return {
                -(tx * m11 + ty * m12 + tz * m13),  // x
                -(tx * m21 + ty * m22 + tz * m23),  // y
                -(tx * m31 + ty * m32 + tz * m33)   // z
            };
        }

        constexpr Vec3 GetPositionFromLocalToParentMatrix() const {
            return {tx, ty, tz};
        }
    };

    inline constexpr Vec3 operator*(const Vec3& v_, const Mat4x3& m_) {
        return {v_.x * m_.m11 + v_.y * m_.m21 + v_.z * m_.m31 + m_.tx, v_.x * m_.m12 + v_.y * m_.m22 + v_.z * m_.m32 + m_.ty,
                v_.x * m_.m13 + v_.y * m_.m23 + v_.z * m_.m33 + m_.tz};
    }

    inline constexpr Vec3& operator*=(Vec3& v_, const Mat4x3& m_) {
        v_ = v_ * m_;
        return v_;
    }
}  // namespace mathUtl

/// <summary>
/// インライン関数で実装できない関数の実装
/// </summary>
namespace mathUtl {
    /// <summary>
    /// 四元数をオイラー角形式に変換する
    /// 入力される四元数は、その名前が示すようにオブジェクト空間から慣性空間、
    /// または、慣性空間からオブジェクト空間への回転を実行するものとする
    /// </summary>
    void EulerAngles::FromObjectToInertialQuaternion(const Quaternion& q_) {
        const double sp{-2.0 * (q_.y * q_.z - q_.w * q_.x)};

        if (std::fabs(sp) > 0.9999) {
            pitch = PiOverTwo * sp;
            heading = std::atan2(-q_.x * q_.z + q_.w * q_.y, 0.5 - q_.y * q_.y - q_.z * q_.z);
            bank = 0.0;
            return;
        }

        pitch = std::asin(sp);
        heading = std::atan2(q_.x * q_.z + q_.w * q_.y, 0.5 - q_.x * q_.x - q_.y * q_.y);
        bank = std::atan2(q_.x * q_.y + q_.w * q_.z, 0.5 - q_.x * q_.x - q_.z * q_.z);
    }
    void EulerAngles::FromInertialToObjectQuaternion(const Quaternion& q_) {
        const double sp{-2.0 * (q_.y * q_.z + q_.w * q_.x)};

        if (std::fabs(sp) > 0.9999) {
            pitch = PiOverTwo * sp;
            heading = std::atan2(-q_.x * q_.z - q_.w * q_.y, 0.5 - q_.y * q_.y - q_.z * q_.z);
            bank = 0.0;
            return;
        }

        pitch = std::asin(sp);
        heading = std::atan2(q_.x * q_.z - q_.w * q_.y, 0.5 - q_.x * q_.x - q_.y * q_.y);
        bank = std::atan2(q_.x * q_.y - q_.w * q_.z, 0.5 - q_.x * q_.x - q_.z * q_.z);
    }

    /// <summary>
    /// 座標変換行列をオイラー角形式に変換する
    /// 入力される行列は、それが示すようにオブジェクト空間からワールド空間、
    /// ワールド空間からオブジェクト空間への座標変換を実行するものとする
    /// この行列は、直交しているものとする
    /// </summary>
    void EulerAngles::FromObjectToWorldMatrix(const Mat4x3& m_) {
        const double sp{-m_.m32};
        if (std::fabs(sp) > 9.99999) {
            pitch = PiOverTwo * sp;
            heading = std::atan2(-m_.m23, m_.m11);
            bank = 0.0;
            return;
        }

        heading = std::atan2(m_.m31, m_.m33);
        pitch = std::asin(sp);
        bank = std::atan2(m_.m12, m_.m22);
    }
    void EulerAngles::FromWorldToObjectMatrix(const Mat4x3& m_) {
        const double sp{-m_.m23};
        if (std::fabs(sp) > 9.99999) {
            pitch = PiOverTwo * sp;
            heading = std::atan2(-m_.m31, m_.m11);
            bank = 0.0;
            return;
        }

        heading = std::atan2(m_.m13, m_.m33);
        pitch = std::asin(sp);
        bank = std::atan2(m_.m21, m_.m22);
    }

    /// <summary>
    /// 回転行列をオイラー角形式に変換する
    /// </summary>
    void EulerAngles::FromRotationMatrix(const RotationMatrix& m_) {
        const double sp{-m_.m23};
        if (std::fabs(sp) > 9.99999) {
            pitch = PiOverTwo * sp;
            heading = std::atan2(-m_.m31, m_.m11);
            bank = 0.0;
            return;
        }

        heading = std::atan2(m_.m13, m_.m33);
        pitch = std::asin(sp);
        bank = std::atan2(m_.m21, m_.m22);
    }

    /// <summary>
    /// オブジェクト空間<->慣性空間の回転を実行するようにセットアップする
    /// 方向はオイラー角形式で与えられる
    /// </summary>
    void Quaternion::SetToRotateObjectToInertial(const EulerAngles& orientation_) {
        double sp{}, sb{}, sh{};
        double cp{}, cb{}, ch{};

        SinCos(&sp, &cp, orientation_.pitch * 0.5);
        SinCos(&sb, &cb, orientation_.bank * 0.5);
        SinCos(&sh, &ch, orientation_.heading * 0.5);

        w = ch * cp * cb + sh * sp * sb;
        x = ch * sp * cb + sh * cp * sb;
        y = -ch * sp * sb + sh * cp * cb;
        z = -sh * sp * cb + ch * cp * sb;
    }
    void Quaternion::SetToRotateInertialToObject(const EulerAngles& orientation_) {
        double sp{}, sb{}, sh{};
        double cp{}, cb{}, ch{};

        SinCos(&sp, &cp, orientation_.pitch * 0.5);
        SinCos(&sb, &cb, orientation_.bank * 0.5);
        SinCos(&sh, &ch, orientation_.heading * 0.5);

        w = ch * cp * cb + sh * sp * sb;
        x = -ch * sp * cb - sh * cp * sb;
        y = ch * sp * sb - sh * cb * cp;
        z = sh * sp * cb - ch * cp * sb;
    }

    void RotationMatrix::Setup(const EulerAngles& orientation_) {
        double sh{}, ch{}, sp{}, cp{}, sb{}, cb{};
        SinCos(&sh, &ch, orientation_.heading);
        SinCos(&sp, &cp, orientation_.pitch);
        SinCos(&sb, &cb, orientation_.bank);

        m11 = ch * cb + sh * sp * sb;
        m12 = -ch * sb + sh * sp * cb;
        m13 = sh * cp;
        m21 = sb * cp;
        m22 = cb * cp;
        m23 = -sp;
        m31 = -sh * cb + ch * sp * sb;
        m32 = sb * sh + ch * sp + cb;
        m33 = ch * cp;
    }

    void RotationMatrix::FromInertialToObjectQuaternion(const Quaternion& q_) {
        m11 = 1.0 - 2.0 * (q_.y * q_.y + q_.z * q_.z);
        m12 = 2.0 * (q_.x * q_.y + q_.w * q_.z);
        m13 = 2.0 * (q_.x * q_.z - q_.w * q_.y);
        m21 = 2.0 * (q_.x * q_.y - q_.w * q_.z);
        m22 = 1.0 - 2.0 * (q_.x * q_.x + q_.z * q_.z);
        m23 = 2.0 * (q_.y * q_.z + q_.w * q_.x);
        m31 = 2.0 * (q_.x * q_.z + q_.w * q_.y);
        m32 = 2.0 * (q_.y * q_.z - q_.w * q_.x);
        m33 = 1.0 - 2.0 * (q_.x * q_.x + q_.y * q_.y);
    }
    void RotationMatrix::FromObjectToInertialQuaternion(const Quaternion& q_) {
        m11 = 1.0 - 2.0 * (q_.y * q_.y + q_.z * q_.z);
        m12 = 2.0 * (q_.x * q_.y - q_.w * q_.z);
        m13 = 2.0 * (q_.x * q_.z + q_.w * q_.y);
        m21 = 2.0 * (q_.x * q_.y + q_.w * q_.z);
        m22 = 1.0 - 2.0 * (q_.x * q_.x + q_.z * q_.z);
        m23 = 2.0 * (q_.y * q_.z - q_.w * q_.x);
        m31 = 2.0 * (q_.x * q_.z - q_.w * q_.y);
        m32 = 2.0 * (q_.y * q_.z + q_.w * q_.x);
        m33 = 1.0 - 2.0 * (q_.x * q_.x + q_.y * q_.y);
    }

    void Mat4x3::SetupLocalToParent(const Vec3& pos_, const EulerAngles& orient_) {
        RotationMatrix orientMatrix;
        orientMatrix.Setup(orient_);

        SetupLocalToParent(pos_, orientMatrix);
    }
    void Mat4x3::SetupLocalToParent(const Vec3& pos_, const RotationMatrix& orient_) {
        m11 = orient_.m11;
        m12 = orient_.m21;
        m13 = orient_.m31;
        m21 = orient_.m12;
        m22 = orient_.m22;
        m23 = orient_.m32;
        m31 = orient_.m13;
        m32 = orient_.m23;
        m33 = orient_.m33;

        tx = pos_.x;
        ty = pos_.y;
        tz = pos_.z;
    }
    void Mat4x3::SetupParentToLocal(const Vec3& pos_, const EulerAngles& orient_) {
        RotationMatrix orientMatrix;
        orientMatrix.Setup(orient_);

        SetupParentToLocal(pos_, orientMatrix);
    }
    void Mat4x3::SetupParentToLocal(const Vec3& pos_, const RotationMatrix& orient_) {
        m11 = orient_.m11;
        m12 = orient_.m12;
        m13 = orient_.m13;
        m21 = orient_.m21;
        m22 = orient_.m22;
        m23 = orient_.m23;
        m31 = orient_.m31;
        m32 = orient_.m32;
        m33 = orient_.m33;

        tx = -(pos_.x * m11 + pos_.y * m21 + pos_.z * m31);
        ty = -(pos_.x * m12 + pos_.y * m22 + pos_.z * m32);
        tz = -(pos_.x * m13 + pos_.y * m23 + pos_.z * m33);
    }

    void Mat4x3::FromQuaternion(const Quaternion& q_) {
        const double ww{2.0 * q_.w};
        const double xx{2.0 * q_.x};
        const double yy{2.0 * q_.y};
        const double zz{2.0 * q_.z};

        m11 = 1.0 - yy * q_.y - zz * q_.z;
        m12 = xx * q_.y + ww * q_.z;
        m13 = xx * q_.z - ww * q_.y;
        m21 = xx * q_.y - ww * q_.z;
        m22 = 1.0 - xx * q_.x - zz * q_.z;
        m23 = yy * q_.z + ww * q_.x;
        m31 = xx * q_.z + ww * q_.y;
        m32 = yy * q_.z - ww * q_.x;
        m33 = 1.0 - xx * q_.x - yy * q_.y;

        tx = 0.0;
        ty = 0.0;
        tz = 0.0;
    }
}  // namespace mathUtl
