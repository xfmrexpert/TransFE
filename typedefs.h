/***************************************************************************
 *   Copyright (C) 2005-2025 by T. C. Raymond                              *
 *   tcraymond@inductivereasoning.com                                      *
 *                                                                         *
 *   Use of this source code is governed by an MIT-style                   *
 *   license that can be found in the LICENSE.txt file or at               *
 *   https://opensource.org/licenses/MIT.                                  *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                  *
 *                                                                         *
 ***************************************************************************/

#pragma once

#include <vector>
#include "./Eigen/Dense"
#include "./Eigen/Sparse"
#include "./Eigen/IterativeLinearSolvers"
#include <complex>
#include <memory>
#include <iostream>

#ifdef NDEBUG
    #define TFEM_ASSERT(condition, message) ((void)0)
#else
    #define TFEM_ASSERT(condition, message) \
        do { \
            if (!(condition)) { \
                std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                          << " line " << __LINE__ << ": " << message << std::endl; \
                std::terminate(); \
            } \
        } while (false)
#endif

using point = Eigen::Matrix<double, 3, 1>;

template <typename T>
using BigMatrix =  Eigen::SparseMatrix<T>;

template <typename T>
using BigVector = Eigen::VectorX<T>;

constexpr double PI = 3.141592653589793238512808959406186204433;

template <typename T>
using Vector = Eigen::VectorX<T>;

template <typename T>
using Matrix = Eigen::MatrixX<T>;