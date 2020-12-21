#pragma once
#include <cstdint>
#include <cmath>

namespace mathUtl {
    using int8 = std::int8_t;
    using int16 = std::int16_t;
    using int32 = std::int32_t;
    using int64 = std::int64_t;
    using uint8 = std::uint8_t;
    using uint16 = std::uint16_t;
    using uint32 = std::uint32_t;
    using uint64 = std::uint64_t;

    constexpr double Pi = 3.14159265;
    constexpr double TwoPi = Pi * 2.0;
    constexpr double PiOverTwo = Pi / 2.0;
    constexpr double OneOverPi = 1.0 / Pi;
    constexpr double OneOverTwoPi = 1.0 / TwoPi;

    /// <summary>
    /// 適切に2πの倍数を加える事で角度を-π...πの範囲にラップする
    /// </summary>
    inline double WrapPi(double theta_) {
        theta_ += Pi;
        theta_ -= std::floor(theta_ * OneOverTwoPi) * TwoPi;
        theta_ -= Pi;
        return theta_;
    }

    /// <summary>
    /// 安全な逆三角関数
    /// </summary>
    inline constexpr double SafeAcos(const double x_) {
        // 制限条件をチェックする
        if (x_ <= -1.0) {
            return Pi;
        }

        if (x_ >= 1.0) {
            return 0.0;
        }

        return std::acos(x_);
    }

    inline void SinCos(double* sin_, double* cos_, const double theta_) {
        *sin_ = std::sin(theta_);
        *cos_ = std::cos(theta_);
    }
}