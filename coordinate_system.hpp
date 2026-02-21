// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include "constants.hpp"
#include "Mesh/point.h"

namespace TFEM
{
    class CoordinateSystem {
    public:
        virtual ~CoordinateSystem() = default;
        virtual double measure(double detJ, const Point& p) const = 0;
    };

    class Cartesian : public CoordinateSystem {
        double measure(double detJ, const Point&) const override { return detJ; }
    };

    class Axisymmetric : public CoordinateSystem {
        double measure(double detJ, const Point& p) const override { 
            return detJ * Constants::TWO_PI * p.x; 
        }
    };
}