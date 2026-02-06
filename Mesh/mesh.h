#pragma once
#include <vector>
#include <span>
#include <array>
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <string>
#include <fstream>

namespace TFEM {

    // Forward Declarations
    class CellView;
    class CellIterator;
    struct CellAccessor;

    // =========================================================================
    // 1. TOPOLOGY DEFINITIONS
    // =========================================================================
    
    enum class ElementType : uint8_t {
        Point = 0,
        Segment = 1,      // 2 Nodes
        Triangle = 2,     // 3 Nodes, 3 Edges, 1 Face
        Quad = 3,         // 4 Nodes, 4 Edges, 1 Face
        Tetrahedron = 4,  // 4 Nodes, 6 Edges, 4 Faces
        Hexahedron = 5,   // 8 Nodes, 12 Edges, 6 Faces
        Prism = 6,        // 6 Nodes, 9 Edges, 5 Faces
        Pyramid = 7
    };

    // =========================================================================
    // 2. MESH DATABASE (The Storage Layer)
    // =========================================================================
    // Core Mesh class using Structure of Arrays (SoA) layout.
    class Mesh {
        // Friends for Iteration
        friend class CellView;
        friend class CellIterator;

    public:
        Mesh() = default;

        // --- Core API ---
        void Clear();
        
        // Populate the mesh from a GMsh file (format 4.1)
        void ReadGmsh(const std::string& filename);

        // Analyzes the cell-to-node connectivity to build:
        // 1. Global unique Edge list
        // 2. Global unique Face list
        // 3. Cell-to-Edge and Cell-to-Face connectivity
        void BuildTopology();

        // --- Iterators ---
        CellAccessor Cells() const;

        // --- Queries ---
        int NumNodes() const { return coordinates.size() / 3; }
        int NumCells() const { return cell_types.size(); }
        int NumEdges() const { return edges.size() / 2; }
        int NumFaces() const { return face_node_offsets.empty() ? 0 : face_node_offsets.size() - 1; }

        // Returns {minX, minY, minZ, maxX, maxY, maxZ}
        std::array<double, 6> BoundingBox() const {
            if (coordinates.empty()) return {0,0,0,0,0,0};
            double minX = coordinates[0], minY = coordinates[1], minZ = coordinates[2];
            double maxX = minX, maxY = minY, maxZ = minZ;
            int n = NumNodes();
            for (int i = 1; i < n; ++i) {
                double x = coordinates[i*3], y = coordinates[i*3+1], z = coordinates[i*3+2];
                if (x < minX) minX = x; if (x > maxX) maxX = x;
                if (y < minY) minY = y; if (y > maxY) maxY = y;
                if (z < minZ) minZ = z; if (z > maxZ) maxZ = z;
            }
            return {minX, minY, minZ, maxX, maxY, maxZ};
        }

    private:
        // --- DATA STORAGE (Structure of Arrays) ---
        
        // Geometry
        // Stored as interleaved X, Y, Z. 
        std::vector<double> coordinates; 

        // Global Topology
        // Edges: Always 2 nodes per edge. Flattened [n1, n2, n1, n2...]
        std::vector<int> edges; 
        
        // Faces: Variable node count (Tri vs Quad vs Poly). CSR format.
        std::vector<int> face_to_node;
        std::vector<int> face_node_offsets; 

        // Cell Connectivity
        std::vector<int> cell_to_node;
        std::vector<int> cell_node_offsets;

        std::vector<int> cell_to_edge;
        std::vector<int> cell_edge_offsets;

        std::vector<int> cell_to_face;
        std::vector<int> cell_face_offsets;

        // Orientations
        // Edge: Does local edge direction match global edge direction? (+1/-1)
        std::vector<int8_t> cell_edge_orientations; 
        // Face: Does local face normal match global face normal? (+1/-1)
        std::vector<int8_t> cell_face_orientations;

        // Metadata
        std::vector<ElementType> cell_types;
        // Physical tags (e.g. 1="Coil", 2="Core"). 
        // Note: Currently assumes one tag per element.
        std::vector<int> cell_attributes;

        // --- Internal Parsing Helpers ---
        void ReadGmsh2(std::ifstream& file);
        void ReadGmsh4(std::ifstream& file);
    };

    // =========================================================================
    // 3. ERGONOMIC VIEWS (The API Layer)
    // =========================================================================

    // A lightweight wrapper for a point in 3D space.
    struct NodeView {
        std::span<const double, 3> data;

        double X() const { return data[0]; }
        double Y() const { return data[1]; }
        double Z() const { return data[2]; }
        
        // Vector arithmetic helpers can be added here
    };

    // The primary interface for your Solvers and Formulations.
    // Created on the stack, discarded immediately. Zero allocation.
    class CellView {
        const Mesh& mesh;
        int index;

    public:
        CellView(const Mesh& m, int i) : mesh(m), index(i) {}

        int ID() const { return index; }
        ElementType Type() const { return mesh.cell_types[index]; }
        int Attribute() const { return mesh.cell_attributes[index]; }

        std::span<const int> Nodes() const {
            assert(index >= 0 && index < mesh.NumCells());
            int start = mesh.cell_node_offsets[index];
            int count = mesh.cell_node_offsets[index + 1] - start;
            assert(start + count <= static_cast<int>(mesh.cell_to_node.size()));
            return { &mesh.cell_to_node[start], static_cast<size_t>(count) };
        }

        NodeView Node(int localIndex) const {
            assert(index >= 0 && index < mesh.NumCells());
            assert((mesh.cell_node_offsets[index] + localIndex) < static_cast<int>(mesh.cell_to_node.size()));

            int globalID = mesh.cell_to_node[mesh.cell_node_offsets[index] + localIndex];
            assert(globalID * 3 + 2 < static_cast<int>(mesh.coordinates.size()));

            return NodeView{ std::span<const double, 3>(mesh.coordinates.data() + globalID * 3, 3) };
        }

        std::span<const int> Edges() const {
            assert(index >= 0 && index < mesh.NumCells());
            int start = mesh.cell_edge_offsets[index];
            if (start >= static_cast<int>(mesh.cell_to_edge.size())) return {}; 
            
            int count = mesh.cell_edge_offsets[index + 1] - start;
            return { &mesh.cell_to_edge[start], static_cast<size_t>(count) };
        }

        std::span<const int8_t> EdgeOrientations() const {
            assert(index >= 0 && index < mesh.NumCells());
            if (mesh.cell_edge_orientations.empty()) return {};

            int start = mesh.cell_edge_offsets[index];
            int count = mesh.cell_edge_offsets[index + 1] - start;
            return { &mesh.cell_edge_orientations[start], static_cast<size_t>(count) };
        }

        std::span<const int> Faces() const {
            assert(index >= 0 && index < mesh.NumCells());
            int start = mesh.cell_face_offsets[index];
            if (start >= static_cast<int>(mesh.cell_to_face.size())) return {};

            int count = mesh.cell_face_offsets[index + 1] - start;
            return { &mesh.cell_to_face[start], static_cast<size_t>(count) };
        }

        std::span<const int8_t> FaceOrientations() const {
            assert(index >= 0 && index < mesh.NumCells());
            if (mesh.cell_face_orientations.empty()) return {};
            
            int start = mesh.cell_face_offsets[index];
            int count = mesh.cell_face_offsets[index + 1] - start;
            return { &mesh.cell_face_orientations[start], static_cast<size_t>(count) };
        }
        
        std::array<double, 3> Centroid() const {
            std::array<double, 3> c = {0,0,0};
            auto nodes = Nodes();
            if (nodes.empty()) return c;
            
            for(int id : nodes) {
                c[0] += mesh.coordinates[id*3];
                c[1] += mesh.coordinates[id*3+1];
                c[2] += mesh.coordinates[id*3+2];
            }
            double scale = 1.0 / nodes.size();
            return {c[0]*scale, c[1]*scale, c[2]*scale};
        }
    };

    // =========================================================================
    // 4. ITERATORS (Syntactic Sugar)
    // =========================================================================
    
    // Enables: for (auto cell : mesh.Cells())
    class CellIterator {
        const Mesh& mesh;
        int idx;
    public:
        CellIterator(const Mesh& m, int i) : mesh(m), idx(i) {}
        CellView operator*() const { return CellView{mesh, idx}; }
        CellIterator& operator++() { ++idx; return *this; }
        bool operator!=(const CellIterator& other) const { return idx != other.idx; }
    };

    struct CellAccessor {
        const Mesh& mesh;
        CellIterator begin() const { return CellIterator(mesh, 0); }
        CellIterator end() const { return CellIterator(mesh, mesh.NumCells()); }
    };

    inline CellAccessor Mesh::Cells() const {
        return CellAccessor{*this};
    }
}