// Copyright (c) 2026 T. C. Raymond
// SPDX-License-Identifier: MIT

#pragma once
#include "json.hpp"
#include "Mesh/mesh.h"

using json = nlohmann::json;

/**
 * @brief Base class for physics solvers using MFEM
 *
 * @warning The mesh and config references must outlive this solver instance.
 *          Do not destroy the mesh or config objects before the solver is done.
 */
class PhysicsSolver {
protected:
    TFEM::Mesh &mesh;
    const json &config;

public:
    PhysicsSolver(TFEM::Mesh &m, const json &c) : mesh(m), config(c) {}
    
    // Virtual destructor is essential for unique_ptr polymorphism
    virtual ~PhysicsSolver() = default;

    virtual void Setup() = 0;
    virtual void Run() = 0;
    virtual void Save() = 0;
};