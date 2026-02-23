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

#include "Mesh/mesh.h"
#include "finite_element.hpp"
#include "element_transform.hpp"

namespace TFEM
{
    class FiniteElementSpace {
    public:
        // FESpace takes ownership of the finite element definition
        FiniteElementSpace(const Mesh* mesh, std::unique_ptr<FiniteElementBase> fe)
            : mesh(mesh), fe(std::move(fe)) 
        {
            Setup();
        }

        // Disables copy to prevent accidental heavy mesh/dof copying
        FiniteElementSpace(const FiniteElementSpace&) = delete;

        virtual ~FiniteElementSpace() = default;

        // Access to the underlying finite element
        const FiniteElementBase* GetFiniteElement(CellView cell) const {
            ElementType type = cell.Type();
        
            auto it = elements.find(type);
            if (it != elements.end()) {
                return it->second.get();
            }
            return nullptr;
        }

        // Retrieve or create the geometric transform for the given cell
        const ElementTransform* GetElementTransform(const CellView& cell) const {
            ElementType type = cell.Type();
            size_t numNodes = cell.Nodes().size();

            // Determine geometric order based on node count (heuristic)
            int order = 1;
            if ((type == ElementType::Triangle && numNodes > 3) ||
                (type == ElementType::Quad && numNodes > 4) ||
                (type == ElementType::Tetrahedron && numNodes > 4) ||
                (type == ElementType::Hexahedron && numNodes > 8) ||
                (type == ElementType::Segment && numNodes > 2)) {
                order = 2;
            }

            // Check cache
            std::pair<ElementType, int> key = { type, order };
            auto it = transforms.find(key);
            if (it != transforms.end()) {
                return it->second.get();
            }

            // Create new transform
            size_t nsd = mesh->SpaceDimension();
            size_t npd = mesh->Dimension(); // Default to mesh dimension

            // Refine parametric dimension based on element type
            if (type == ElementType::Segment) npd = 1;
            else if (type == ElementType::Point) npd = 0;
            else if (type == ElementType::Triangle || type == ElementType::Quad) npd = 2;
            else if (type == ElementType::Tetrahedron || type == ElementType::Hexahedron) npd = 3;

            transforms[key] = std::make_unique<ElementTransform>(nsd, npd, order);
            return transforms[key].get();
        }

        const Mesh* GetMesh() const { return mesh; }

        // The critical method for assembly
        // Given an element index, return the global equation numbers
        const std::vector<int>& GetElementDofs(int elemIdx) const
        {

        }

    protected:
        const Mesh* mesh;
        const std::unique_ptr<FiniteElementBase> fe;

        // Cache for finite elements: <ElementType, Order> -> Transform
        mutable std::map<ElementType, std::unique_ptr<FiniteElementBase>> elements;

        // Cache for element transforms: <ElementType, Order> -> Transform
        mutable std::map<std::pair<ElementType, int>, std::unique_ptr<ElementTransform>> transforms;

        // Finalize the numbering (call after mesh is ready)
        void Setup()
        {

        }
    };
}

