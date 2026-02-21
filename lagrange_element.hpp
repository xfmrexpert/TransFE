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
#include "finite_element.hpp"
#include "element_transform.hpp"

namespace TFEM
{
    /// Lagrange finite element (scalar, H^1 conforming).
    class LagrangeElement : public FiniteElement<LagrangeShapeFunction, ElementTransform>
    {
    public:
        LagrangeElement(size_t dim, int order = 1)
            : FiniteElement<LagrangeShapeFunction, ElementTransform>(dim, dim, order)
        {
        }

        virtual ~LagrangeElement() = default;

        int numLocalDOFs() const override {
            return this->shape_function.N(Point()).size();
        }

    private:
    
    };
}
