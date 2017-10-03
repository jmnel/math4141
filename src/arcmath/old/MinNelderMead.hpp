#pragma once

#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <utility>
#include <vector>

#include <CoreMath.hpp>
#include <Vec.hpp>

using std::array;
using std::cout;
using std::endl;
using std::function;
using std::make_pair;
using std::pair;
using std::vector;

namespace arc {

    typedef array<Vec2d, 3> Simplex2d;

    /**
     *
     * @brief Storage object for Nelder-Mead algorithm internal state.
     *
     **/
    struct NelderMeadAlgoInfo {
        Simplex2d simplex0;
        size_t maxN;
        double tol, alpha, beta, gamma, sigma;
        Vec2d result;
        vector<Simplex2d> simplices;
        vector<Vec2d> baryCenters;

        size_t opCountReflect = 0;
        size_t opcountExpand = 0;
        size_t opCountContract = 0;
        size_t opCountShrink = 0;

        enum class OpType { Reflect, Expand, Contract, Shrink };
        vector<OpType> ops;
    };

    /**
     *
     * @brief Minimizes function using Nelder-Mead polytope algorithm.
     * @brief Iterates f
     * @param obj Objective function to minimize
     * @param initialSimplex Starting polytop
     * @param tol Desired tolerance for solution
     * @param maxIterations Maximum number of iterations
     * @param alpha Reflection coefficient $\alpha$
     * @param beta Reflection coefficient $\beta$
     * @param gamma Contraction coefficient $\gamma$
     * @param verbose Enable verbose debug info
     * @return Solution to minimization of objective function
     *
     **/
    Vec2d minimizeNelderMead(const function<double(Vec2d)> &obj,
                             Simplex2d initialSimplex,
                             double tol,
                             size_t maxIterations,
                             double aplha,
                             double beta,
                             double gamma,
                             double sigma,
                             NelderMeadAlgoInfo *info,
                             bool verbose = false);
}  // namespace arc
