# Project Guidelines

## Code Style
- **Standard**: C++20.
- **Paradigm**: Prefer Object-Oriented Programming (OOP) with sensible performance considerations.
    - Encapsulate logic with data (Classes > Structs).
    - Use smart pointers (`std::unique_ptr`, `std::shared_ptr`) for ownership.
    - Start with `virtual` methods for polymorphism; optimize only when proven hot.
- **Naming**:
    - Classes: `PascalCase`.
    - Methods: `camelCase`.
    - Member variables: `PascalCase` or `camelCase` (match local file context).
- **Formatting**:
    - `#pragma once` for header guards.

## Architecture
- **Mesh Systems** (Distinct roles):
    - **`MeshDB/` (Core)**: Original OOP mesh database code that is heavy and doesn't handle orientations well.  To be deprecated.**
    - **`Mesh/` (`TFEM` namespace)**: A lightweight, flat-array mesh in development. Optimized for speed but with high code readability and ease of understanding.
- **Components**:
    - `FEProg`: Main solver driver.
    - `Eigen/`: Use this vendored library for all linear algebra (vectors, matrices).
    - `MeshViewer/`: Raylib-based visualization tool.
- **Excluded**:
    - `Old/`: Deprecated code. Do not use or reference.

## Build and Test
- **Build System**: CMake (min 3.12) with vcpkg integration.
- **Build**:
    ```bash
    cmake -B build -S . -DCMAKE_TOOLCHAIN_FILE=vcpkg/scripts/buildsystems/vcpkg.cmake
    cmake --build build
    ```
- **Targets**:
    - `TransFE`: Main solver executable.
    - `MeshViewer`: Visualization tool (requires Raylib).

## Integration Points
- **Eigen**: Header-only linear algebra.
- **Raylib**: Visualization (managed by vcpkg).
- **Gmsh**: Primary mesh input format (parsed/supported in `Mesh/`).

## Project Conventions
- **"New Age" check**: Avoid excessive template metaprogramming or pure Data-Oriented Design (DOD) separation unless strictly required for IO/GPU interop. 
- **Error Handling**: Exceptions (`std::runtime_error`) for setup/IO failures. Assertions for logic errors in debug.
