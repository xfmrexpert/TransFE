// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once

namespace TFEM
{
    template <typename T>
    class GridFunction
    {
    public:
        GridFunction(FiniteElementSpace* fes)
            : fespace(fes)
        {
            if (!fespace) {
                throw std::runtime_error("GridFunction initialized with null FiniteElementSpace");
            }
            // Resize data to match the number of DOFs in the space
            // TODO: This won't get set properly until fespace->SetupGlobalDofs gets called
            data.resize(fespace->GetNDofs());
            data.setZero();
        }

        // Allow access to the underlying data vector for solvers to write into
        Vector<T>& GetData() { return data; }
        const Vector<T>& GetData() const { return data; }

        /// Evaluates the field at a physical point.
        /// Returns std::nullopt if the point is outside the mesh.
        std::optional<T> EvaluateAtPt(const Point& pt) const
        {
            // Find which element contains this point
            auto cell_opt = fespace->GetMesh()->FindElement(pt); 
            
            if (!cell_opt) {
                return std::nullopt; // Point not in mesh
            }

            const auto& cell = *cell_opt;

            // Map physical point to reference coordinates (xi, eta, zeta)
            // This requires the inverse Jacobian / Newton-Raphson for non-linear elements
            const ElementTransform& trans = fespace->GetElementTransform(cell.index);
            Point ref_pt = trans.MapPhysicalToReference(pt, cell);

            // Get Shape functions values at reference point
            // N = [N1, N2, N3, ...]
            const Vector<double>& shape_vals = fespace->GetShapeFunctions(cell.index, ref_pt);

            // Get the DOFs indices for this element
            // dofs = [idx1, idx2, idx3, ...]
            const std::vector<int>& dof_indices = fespace->GetElementDofIndices(cell.index);

            // Interpolate: Sum( N_i * u_i )
            T result;
            if constexpr (std::is_arithmetic_v<T>) {
               result = 0; // Scalar zero
            } else {
               result.setZero(); // Vector/Matrix zero
            }

            for (size_t i = 0; i < dof_indices.size(); ++i) {
                result += shape_vals(i) * data(dof_indices[i]);
            }

            return result;
        }

    private:
        FiniteElementSpace* fespace; //Non-owning pointer
        Vector<T> data;
    };
}