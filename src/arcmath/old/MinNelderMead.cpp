#include "MinNelderMead.hpp"

using std::stable_sort;

namespace arc {

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

        cout << "Running Nelder-Mead minimzation ..." << endl;
        cout << "-----------------------------------" << endl;
        cout << "α=" << alpha << endl;
        cout << "β=" << beta << endl;
        cout << "γ=" << gamma << endl;
        cout << "with max iterations N=" << maxIterations << endl;
        cout << "to tolerance TOL ε=" << tol << endl;
        cout << "-----------------------------------" << endl;
        cout << endl;

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

        array<pair<Vec2d, double>, 3> pairs = {
            {make_pair(initialSimplex[0], 0.0),
             make_pair(initialSimplex[1], 0.0),
             make_pair(initialSimplex[2], 0.0)}};

        auto &x0 = pairs[0].first;
        auto &x1 = pairs[1].first;
        auto &x2 = pairs[2].first;

        // Sort polytope
        stable_sort(objAtX.begin(), objAtX.end(),
                    [](const ObjAtX &a, const ObjAtX &b) -> bool {
                        return a.f <= b.f;
                    });

        // Main algorithm iteration loop
        for (int k = 0; k < (int)n; k++) {
            // Set some reference variables to simply notation
            auto &x0 = objAtX[0].x;
            auto &x1 = objAtX[1].x;
            auto &x2 = objAtX[2].x;
            auto &f0 = objAtX[0].f;
            auto &f1 = objAtX[1].f;
            auto &f2 = objAtX[2].f;

            polytopes.push_back(x0);
            polytopes.push_back(x1);
            polytopes.push_back(x2);

            // Calculate midpoint of best simplex side
            auto xC = 0.5 * (x0 + x1);

            centers.push_back(xC);

            // Calculate reflection across best side
            auto xR = xC + alpha * (xC - x2);
            auto fR = obj(xR);

            // Case 1: f_1 <= f^r < f_n
            if (f0 <= fR && fR < f1) {
                // cout << "case 1" << endl;
                x2 = xR;
                f2 = fR;
            }
            // Case 2: f^r < f_1
            else if (fR < f0) {
                // cout << "Case 2" << endl;
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
            }
            // Case 3: f^r > f_n : the polytope seems too big
            else if (fR >= f1) {
                // cout << "Case 3: Polytope too big." << endl;
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
                        // This doesn't do anything skip it
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
                // cout << "Error: No case" << endl;
            }

            // Sort polytop
            stable_sort(objAtX.begin(), objAtX.end(),
                        [](const ObjAtX &a, const ObjAtX &b) -> bool {
                            return a.f <= b.f;
                        });

            // Check termination condition
            double sum = 0.0;
            // Recalculate midpoint on best side
            xC = 0.5 * (objAtX[0].x + objAtX[1].x);
            auto fC = obj(xC);

            for (auto objX : objAtX) {
                auto diff = abs(obj(objX.x) - fC);
                diff = diff * diff;
                sum += diff;
            }

            auto measure = sum * (1.0 / 3.0);

            if (measure < tol * tol) {
                cout << "Tolerance reached after k=" << k + 1 << " iterations."
                     << endl;
                auto solution = (objAtX[0].x + objAtX[1].x + objAtX[2].x) / 3.0;
                cout << "Solution x˟=" << solution << "ᵀ" << endl;
                cout << "f(x˟)=" << obj(solution);
                cout << "ᵀ" << endl;

                return solution;
            }
        }

        cout << "Error: Failed to find solution in n=" << n;
        cout << " iterations" << endl;
        return Vec2d(0.0);
    }
}  // namespace arc
