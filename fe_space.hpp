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
            TFEM::GeometryType type = cell.type(); 
        
            auto it = elements.find(type);
            if (it != elements.end()) {
                return it->second.get();
            }
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

        // Finalize the numbering (call after mesh is ready)
        void Setup()
        {

        }
    };
}

