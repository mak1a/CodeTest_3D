#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <string>
#include <iomanip>
// 3D座標系はこのファイル
#include "include/Transformer.hpp"

class Astro {
private:
    mathUtl::Vec3 m_shipPos;
    mathUtl::Vec3 m_garbageLocalPos;

    mathUtl::RotationMatrix m_shipRotateMatrix;

    mathUtl::Vec3 m_answerPos;

public:
    Astro() : m_answerPos(0.0, 0.0, 0.0) {
        std::cin >> m_garbageLocalPos;
        std::cin >> m_shipPos;

        mathUtl::Vec3 vx;
        mathUtl::Vec3 vy;
        mathUtl::Vec3 vz;

        std::cin >> vx >> vy >> vz;

        m_shipRotateMatrix = mathUtl::RotationMatrix(vx.x, vx.y, vx.z, vy.x, vy.y, vy.z, vz.x, vz.y, vz.z);
    }

    Astro& Calculate() {
        mathUtl::EulerAngles e;
        e.FromRotationMatrix(m_shipRotateMatrix);
        mathUtl::Quaternion q;
        q.SetToRotateInertialToObject(e);
        mathUtl::Mat4x3 m;
        m.FromQuaternion(q);
        m_answerPos.Set(m_garbageLocalPos * m + m_shipPos);

        return *this;
    }

    void Answer() const {
        std::cout << std::fixed << std::setprecision(10) << m_answerPos.x << " " << m_answerPos.y << " " << m_answerPos.z << std::endl;
    }
};

int main() {
    Astro astro;
    astro.Calculate().Answer();

    return 0;
}

// 入力例
// 20.0 20.0 100.0
// 100.0 200.0 300.0
// 0.70711 0.70711 0.0 - 0.70711 0.70711 0.0 0.0 0.0 1.0
