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

#include <type_traits>
#include <concepts>
#include "Mesh/point.h"
#include "typedefs.h"
#include "shape_function.hpp"
#include "element_transform.hpp"

namespace TFEM
{
    class FiniteElementBase {
    public:
        virtual ~FiniteElementBase() = default;

        // Pure virtual functions to be implemented by derived classes
        virtual int referenceDimensions() const = 0;
        virtual int numLocalDOFs() const = 0;

        // Access to shape functions via base class pointers
        virtual const ShapeFunction* ShapeFunction() const = 0;
    };


    template <typename SF, typename ET>
        requires std::derived_from<SF, ShapeFunction> && std::derived_from<ET, ElementTransform>
    class FiniteElement : public FiniteElementBase {

    public:
        FiniteElement(int spatial_dim, int ref_dim, int order = 1)
            : shape_function(ref_dim, order), transform(std::make_unique<ET>(spatial_dim, ref_dim, order)) { };

        virtual ~FiniteElement() = default;

        int referenceDimensions() const override {
            return transform->referenceDimensions();
        }

        int numLocalDOFs() const override {
            return shape_function.numShapeFunctions();
        }

        const SF* ShapeFunction() const override {
            return &shape_function;
        }

    protected:
        SF shape_function;

    };
}

