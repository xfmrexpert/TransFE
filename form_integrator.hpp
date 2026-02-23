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

#include <iostream>
#include "typedefs.h"
#include "Mesh/point.h"
#include "Mesh/mesh.h"
#include "fe_space.hpp"
#include "element_transform.hpp"
#include "integration_rule.hpp"

namespace TFEM
{
    struct IntegrationPointData
    {
        Point ptPhys;                        // Physical coordinate (x,y,z)
        double measure;                      // Differential measure (detJ * weight * geometry_factor)
        const Matrix<double>& dNdx;          // Physical Gradients (dN/dx)
        const Vector<double>& N;             // Shape Function values
    };

    template<typename Scalar = double>
    class FormIntegrator
    {
    public:
        using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        using VectorType = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

        FormIntegrator(std::shared_ptr<CoordinateSystem> cs, FiniteElementSpace* feSpace) : coords(cs), feSpace(feSpace) {};

        virtual ~FormIntegrator() = default;

        void evaluate(TFEM::CellView cell, MatrixType& Ke)
        {
            const auto* fe = feSpace->GetFiniteElement(cell);
            int nDofs = fe->numLocalDOFs();

            // Resize and zero the output matrix
            Ke.setZero(nDofs, nDofs);

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

                Ke += evaluatePt(int_data);
            }
        };


    protected:
        std::shared_ptr<CoordinateSystem> coords;
        FiniteElementSpace* feSpace;

        virtual MatrixType evaluatePt(const IntegrationPointData& data) const = 0;
    };

    using BilinearFormIntegrator = FormIntegrator<double>;
    using SesquilinearFormIntegrator = FormIntegrator<std::complex<double>>;
} // namespace TFEM