#include "MinNelderMead.hpp"

using std::stable_sort;

namespace arc {

    using OpType = NelderMeadAlgoInfo::OpType;

    // -- minimizeNelderMead function --
    Vec2d minimizeNelderMead(const function<double(Vec2d)> &obj,
                             Simplex2d initialSimplex,
                             double tol,
                             size_t maxIterations,
                             double alpha,
                             double beta,
                             double gamma,
                             double sigma,
                             NelderMeadAlgoInfo *info,
                             bool verbose) {
        assert(tol > 0.0);
        assert(maxIterations > 0);

        // cout << "Running Nelder-Mead minimzation ..." << endl;
        // cout << "-----------------------------------" << endl;
        // cout << "α=" << alpha << endl;
        // cout << "β=" << beta << endl;
        // cout << "γ=" << gamma << endl;
        // cout << "with max iterations N=" << maxIterations << endl;
        // cout << "to tolerance TOL ε=" << tol << endl;
        // cout << "-----------------------------------" << endl;
        // cout << endl;

        // struct ObjAtX {
        // Vec2d x;
        // double f;
        // ObjAtX() = default;
        // ObjAtX(Vec2d const &x, double f) {
        // this->x = x;
        // this->f = f;
        //}
        //};

        auto n = maxIterations;

        // Parameter 1: Reflection coefficient
        // double alpha = 1.0;
        // Parameter 2: Expansion coefficient
        // double beta = 2.0;
        // Parameter 3: Contraction coefficient
        // double gamma = 0.5;
        // auto x = initialSimplex;

        using Pair = pair<Vec2d, double>;

        array<Pair, 3> pairs = {{make_pair(initialSimplex[0], 0.0),
                                 make_pair(initialSimplex[1], 0.0),
                                 make_pair(initialSimplex[2], 0.0)}};

        auto &x0 = pairs[0].first;
        auto &x1 = pairs[1].first;
        auto &x2 = pairs[2].first;

        auto &f0 = pairs[0].second;
        auto &f1 = pairs[1].second;
        auto &f2 = pairs[2].second;

        f0 = obj(x0);
        f1 = obj(x1);
        f2 = obj(x2);

        // cout << "Before:" << endl;
        // cout << x0 << "=" << f0 << endl;
        // cout << x1 << "=" << f1 << endl;
        // cout << x2 << "=" << f2 << endl;
        // cout << endl;

        // Sort initial simplex before starting.
        stable_sort(pairs.begin(), pairs.end(),
                    [](const Pair &p0, const Pair &p1) -> bool {
                        return p0.second <= p1.second;
                    });

        // -- Main algorithm iteration loop --
        for (int k = 0; k < (int)n; k++) {
            // Output internal state info
            info->simplices.emplace_back();
            auto &s = info->simplices.back();
            s[0] = x0;
            s[1] = x1;
            s[2] = x2;

            info->fValues.emplace_back();
            info->fValues.back()[0] = f0;
            info->fValues.back()[1] = f1;
            info->fValues.back()[2] = f2;

            // Calculate midpoint of best side
            auto xC = 0.5 * (x0 + x1);

            // Calculate reflection over best side
            auto xR = xC + alpha * (xC - x2);
            auto fR = obj(xR);

            // Case 1: f_1 <= f^r < f_n
            if (f0 <= fR && fR < f1) {
                x2 = xR;
                f2 = fR;

                info->opCountReflect++;
                info->ops.push_back(OpType::Reflect);
            }
            // Case 2
            else if (fR < f0) {
                // Calculate extrapolated point x^e
                auto xE = xC + beta * (xR - xC);
                auto fE = obj(xE);

                if (fE < fR) {
                    x2 = xE;
                    f2 = fE;
                } else {
                    x2 = xR;
                    f2 = fR;
                }
                info->ops.push_back(OpType::Expand);
            }
            // Case 3: f^r > f_n : the simplex seems too big
            else if (fR >= f1) {
                info->ops.push_back(OpType::Contract);
                if (fR >= f2) {
                    auto xK = xC + gamma * (x2 - xC);
                    auto fK = obj(xK);
                    if (fK < f2) {
                        x2 = xK;
                        f2 = fK;
                    } else {
                        x1 = 0.5 * (x1 + x0);
                        x2 = 0.5 * (x2 + x0);
                        f1 = obj(x1);
                        f2 = obj(x2);
                        assert(x0 == 0.5 * (x0 + x0));
                        // This is not necessary skip it
                        // x0 = 0.5 * (x0 + x0);
                    }
                }
                //
                else if (fR < f2) {
                    auto xK = xC + gamma * (xR - xC);
                    auto fK = obj(xK);
                    if (fK <= fR) {
                        x2 = xK;
                        f2 = fK;
                    } else {
                        x1 = 0.5 * (x1 + x0);
                        x2 = 0.5 * (x2 + x0);
                        f1 = obj(x1);
                        f2 = obj(x2);
                    }
                }
            } else {
                assert(false);
            }

            // Sort simplex vertices.
            stable_sort(pairs.begin(), pairs.end(),
                        [](const Pair &p0, const Pair &p1) -> bool {
                            return p0.second <= p1.second;
                        });

            // Check stop condition.

            // Recalculate midpoint on best side
            xC = 0.5 * (x0 + x1);
            auto fC = obj(xC);

            auto diff0 = std::abs(f0 - fC);
            auto diff1 = std::abs(f1 - fC);
            auto diff2 = std::abs(f2 - fC);
            auto sum = diff0 * diff0 + diff1 * diff1 + diff2 * diff2;
            auto measureSq = sum / 3.0;

            if (measureSq < tol * tol) {
                cout << "done" << endl;
                auto solution = (x0 + x1 + x2) / 3.0;
                return solution;
            }

            if (k + 1 >= (int)n) {
                cout << "Warning: Failed to reach TOL in n=" << n;
                cout << " iterations!" << endl;
                auto solution = (x0 + x1 + x2) / 3.0;
                return solution;
            }
        }


        return Vec2d::zero();
    }
}  // namespace arc
