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

#include "typedefs.h"
#include "Mesh/point.h"
#include "Mesh/mesh.h"
#include "coordinate_system.hpp"
#include "fe_space.hpp"
#include "element_transform.hpp"
#include "integration_rule.hpp"
#include "assembler.h"

namespace TFEM
{
    template<typename Scalar = double>
    class LinearFormIntegrator
    {
    public:
        using VectorType = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

        LinearFormIntegrator(std::shared_ptr<CoordinateSystem> cs, FiniteElementSpace* feSpace) : coords(cs), feSpace(feSpace) { };
        
        ~LinearFormIntegrator() = default;
        
        void evaluate(TFEM::CellView cell, VectorType& f) const
        {
            const auto* fe = feSpace->GetFiniteElement(cell);
            size_t nDofs = fe->numLocalDOFs();

            f.setZero(nDofs);

            const auto* integration = fe.getIntegrationRule();

            const auto& quadPointsRef = integration->IntPts();
            const auto& quadWeights = integration->Weights();

            // Retrieve geometric transform from FE Space
            const auto* transform = feSpace->GetElementTransform(cell);

            // For each quadrature point in the reference domain
            for (size_t q = 0; q < integration->numIntPts(); q++)
            {
                IntegrationPointData int_data;
                const Point& ptRef = quadPointsRef[q];
                double wq = quadWeights[q];

                // Transform: Reference Point -> Physical Point
                int_data.ptPhys = transform->mapReferencePointToPhysical(ptRef, cell); // x(xi)

                // Calculate Jacobian and its Determinant
                // J = dx/dxi
                auto J = transform->Jacobian(ptRef, cell);
                double detJ = J.determinant();

                // Calculate Physical Gradients
                // dN/dx = dN/dxi * J^{-1}
                // Transpose logic: dN/dx (row vectors) = (J^{-T} * (dN/dxi)^T)^T = dN/dxi * J^{-1}
                // Eigen representation: shapeFunction->Gradients returns [dNi/dxi]
                auto Jinv = J.inverse();
                auto dNdxi = fe.ShapeFunction()->Gradients(ptRef);

                // Assuming Gradients returns (nShape x nDim), we post-multiply by Jinv
                int_data.dNdx = dNdxi * Jinv;

                // Shape Function Values
                int_data.N = fe.ShapeFunction()->Values(ptRef);

                // Compute Total Measure (Weight * DetJacobian * CoordinateSystems)
                int_data.measure = wq * coords->measure(detJ, int_data.ptPhys);

                f += evaluatePt(int_data);

            }
            //std::cout << "Adding force contributor:\n";
            //std::cout << f;
        };

    protected:
		std::shared_ptr<CoordinateSystem> coords;
        FiniteElementSpace* feSpace;

        virtual Vector<double> evaluatePt(const IntegrationPointData& quadData) const = 0;
    };
}