// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

#include "typedefs.h"
#include "coordinate_system.hpp"
#include "form_integrator.hpp"
#include "constants.hpp"
#include "Mesh/point.h"
#include "finite_element.hpp"
#include "coefficient.hpp"

namespace TFEM
{
    // -----------------------------------------------------------------------------
    // 1. Stiffness Integrator: -div( eps * grad(u) )
    // -----------------------------------------------------------------------------
    // Solves: Integral( eps * grad(u) . grad(v) * 2*pi*r * dr * dz )
    class DiffusionIntegrator : public BilinearFormIntegrator
    {
    private:
        Coefficient *Q;

    public:
        DiffusionIntegrator(std::shared_ptr<CoordinateSystem> cs, FiniteElementSpace* feSpace, Coefficient* q) 
        : BilinearFormIntegrator(cs, feSpace), Q(q)
        {
            TFEM_ASSERT(Q != nullptr, "Coefficient cannot be null");
        };

        MatrixType evaluatePt(const IntegrationPointData& data) const override
        {
            // Get dimensions from the data provided by the FE
            // data.dNdx is (nDofs x SDim)
            const int nDofs = data.dNdx.rows();

            double coeff = Q->EvaluateAtPt(data.ptPhys);

            MatrixType Ke = (coeff * data.measure) * (data.dNdx * data.dNdx.transpose());
            //std::cout << "Ke: " << std::endl;
            //std::cout << Ke << std::endl;
            return Ke;
        };

    };
}