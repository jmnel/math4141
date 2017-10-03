// -- Calculations for assignment 1 --
// -- Chapter 3: Exercise 1 ----------

#include <cstdlib>
#include <iostream>
#include <vector>

#include <Vec.hpp>

#include <MinNelderMead.hpp>

using std::cout;
using std::endl;
using std::string;
using std::vector;

using namespace arc;

int main(int argc, char* argv[]) {
    const double alpha = 1.0;
    const double beta = 2.0;
    const double gamma = 0.5;
    const double sigma = 0.5;

    Simplex2d s0;
    s0[0] = Vec2(0.1, 0.0);
    s0[1] = Vec2(0.0, 0.1);
    s0[2] = Vec2(0.0, 0.0);

    const double tol = 1e-24;

    const size_t n = 2;

    NelderMeadAlgoInfo info;

    auto f = [](Vec2d x) -> double {
        return 2.0 * x[0] * x[0] - x[0] - 2.0 * x[1] + x[1] * x[1] + 4.0;
    };

    auto r = minimizeNelderMead(f, s0, tol, n, alpha, beta, gamma, sigma, &info,
                                true);

    cout << "r=" << r << endl;

    for (size_t i = 0; i < info.ops.size(); i++) {
        auto x0 = info.simplices[i][0];
        auto x1 = info.simplices[i][1];
        auto x2 = info.simplices[i][2];
        auto f0 = info.fValues[i][0];
        auto f1 = info.fValues[i][1];
        auto f2 = info.fValues[i][2];
        //cout << i << " - " << (int)info.ops[i] << endl;
        //cout << "  x0=" << x0 << ", x1=" << x1 << ", " << x2 << endl;
        //cout << "  f0=" << f0 << ", f1=" << f1 << ", f2=" << f2 << endl;
    }
}
