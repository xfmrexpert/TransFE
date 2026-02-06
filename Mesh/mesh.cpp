#include "mesh.h"
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <span>

namespace TFEM {

    void Mesh::Clear() {
        coordinates.clear();
        edges.clear();
        face_to_node.clear();
        face_node_offsets.clear();
        cell_to_node.clear();
        cell_node_offsets.clear();
        cell_to_edge.clear();
        cell_edge_offsets.clear();
        cell_to_face.clear();
        cell_face_offsets.clear();
        cell_edge_orientations.clear();
        cell_face_orientations.clear();
        cell_types.clear();
        cell_attributes.clear();
    }

    // --- Topology Building Helpers ---

    struct EdgeKey {
        int n1, n2;
        bool operator<(const EdgeKey& other) const {
            if (n1 != other.n1) return n1 < other.n1;
            return n2 < other.n2;
        }
    };

    struct FaceKey {
        std::vector<int> nodes;
        bool operator<(const FaceKey& other) const {
            return nodes < other.nodes;
        }
    };

    static const std::vector<std::pair<int, int>> TetEdges = {
        {0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}
    };
    static const std::vector<std::vector<int>> TetFaces = {
        {0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}
    };
    
    static const std::vector<std::pair<int, int>> HexEdges = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}, // Base
        {4, 5}, {5, 6}, {6, 7}, {7, 4}, // Top
        {0, 4}, {1, 5}, {2, 6}, {3, 7}  // Pillars
    };
    static const std::vector<std::vector<int>> HexFaces = {
        {0, 3, 2, 1}, // Bottom
        {4, 5, 6, 7}, // Top
        {0, 1, 5, 4}, // Front
        {1, 2, 6, 5}, // Right
        {2, 3, 7, 6}, // Back
        {3, 0, 4, 7}  // Left
    };

    static void GetLocalEdges(ElementType type, std::vector<std::pair<int, int>>& edges) {
        edges.clear();
        switch (type) {
            case ElementType::Triangle:
                edges = {{0, 1}, {1, 2}, {2, 0}};
                break;
            case ElementType::Quad:
                edges = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
                break;
            case ElementType::Tetrahedron:
                edges = TetEdges;
                break;
            case ElementType::Hexahedron:
                edges = HexEdges;
                break;
            default:
                break;
        }
    }

    static void GetLocalFaces(ElementType type, std::vector<std::vector<int>>& faces) {
        faces.clear();
        switch (type) {
            case ElementType::Tetrahedron:
                faces = TetFaces;
                break;
            case ElementType::Hexahedron:
                faces = HexFaces;
                break;
            default:
                break;
        }
    }

    // Helper to determine relative face orientation
    // Returns 1 if 'current' matches 'canonical' (rotated)
    // Returns -1 if 'current' matches 'canonical' reversed
    // Returns 0 if no match found (topology error)
    static int8_t DetermineFaceOrientation(const std::vector<int>& current, const std::vector<int>& canonical) {
        if (current.size() != canonical.size()) return 0;
        size_t n = current.size();
        if (n == 0) return 0;

        // Find match for the first node
        auto it = std::find(canonical.begin(), canonical.end(), current[0]);
        if (it == canonical.end()) return 0; // Node mismatch

        size_t startIdx = std::distance(canonical.begin(), it);

        // Check forward (Orientation +1)
        bool matchForward = true;
        for (size_t i = 1; i < n; ++i) {
            if (current[i] != canonical[(startIdx + i) % n]) {
                matchForward = false;
                break;
            }
        }
        if (matchForward) return 1;

        // Check backward (Orientation -1)
        // If current is reversed, current[1] should match canonical[(startIdx - 1 + n) % n]
        bool matchBackward = true;
        for (size_t i = 1; i < n; ++i) {
            if (current[i] != canonical[(startIdx - i + n) % n]) {
                matchBackward = false;
                break;
            }
        }
        if (matchBackward) return -1;

        return 0; // Nodes match set, but order/permutation is invalid
    }

    void Mesh::BuildTopology() {
        // 1. Reset Global Topology
        edges.clear();
        face_to_node.clear();
        face_node_offsets.clear();
        face_node_offsets.push_back(0);

        cell_to_edge.clear();
        cell_edge_offsets.clear();
        cell_edge_offsets.push_back(0);

        cell_to_face.clear();
        cell_face_offsets.clear();
        cell_face_offsets.push_back(0);

        cell_edge_orientations.clear();
        cell_face_orientations.clear();

        // 2. Maps for Uniqueness
        // Using std::map for simplicity. 
        // Logic: FaceKey stores sorted nodes to identify the manifold geometry.
        // We recover the canonical ordering from the stored global definition.
        std::map<EdgeKey, int> edgeMap;
        std::map<FaceKey, int> faceMap;

        // 3. Iterate Cells
        for (auto cell : Cells()) {
            auto nodes = cell.Nodes();
            ElementType type = cell.Type();

            // --- Handle Edges ---
            std::vector<std::pair<int, int>> localEdges;
            GetLocalEdges(type, localEdges);

            for (const auto& localEdge : localEdges) {
                int n1 = nodes[localEdge.first];
                int n2 = nodes[localEdge.second];

                EdgeKey key = (n1 < n2) ? EdgeKey{n1, n2} : EdgeKey{n2, n1};
                
                int globalID;
                if (edgeMap.find(key) == edgeMap.end()) {
                    globalID = (int)edgeMap.size();
                    edgeMap[key] = globalID;
                    // Store in global edge definitions
                    edges.push_back(key.n1); // Always stored min->max technically
                    edges.push_back(key.n2);
                } else {
                    globalID = edgeMap[key];
                }

                cell_to_edge.push_back(globalID);
                
                // If local n1 == key.n1 (min), then local is min->max. Global is min->max. Match (+1).
                int8_t orientation = (n1 == key.n1) ? 1 : -1;
                cell_edge_orientations.push_back(orientation);
            }
            cell_edge_offsets.push_back((int)cell_to_edge.size());


            // --- Handle Faces ---
            std::vector<std::vector<int>> localFaces;
            GetLocalFaces(type, localFaces);

            for (const auto& localFace : localFaces) {
                // Construct the local face global node indices
                std::vector<int> currentFaceNodes;
                currentFaceNodes.reserve(localFace.size());
                for (int idx : localFace) {
                    currentFaceNodes.push_back(nodes[idx]);
                }

                // Create Key for Uniqueness (Sorted Nodes)
                FaceKey key;
                key.nodes = currentFaceNodes;
                std::sort(key.nodes.begin(), key.nodes.end());

                int globalID;
                std::vector<int> canonicalNodes;

                if (faceMap.find(key) == faceMap.end()) {
                    // New Face Found
                    globalID = (int)faceMap.size();
                    faceMap[key] = globalID;

                    // Store UN-SORTED canonical definition from this first encounter
                    // This defines the "+1" orientation for this face.
                    for (int n : currentFaceNodes) face_to_node.push_back(n);
                    face_node_offsets.push_back((int)face_to_node.size());
                    
                    canonicalNodes = currentFaceNodes;
                } else {
                    // Existing Face
                    globalID = faceMap[key];
                    
                    // Retrieve Canonical Definition
                    int start = face_node_offsets[globalID];
                    int end = face_node_offsets[globalID+1];
                    canonicalNodes.assign(face_to_node.begin() + start, face_to_node.begin() + end);
                }

                cell_to_face.push_back(globalID);

                // Determine Orientation
                int8_t orientation = DetermineFaceOrientation(currentFaceNodes, canonicalNodes);
                cell_face_orientations.push_back(orientation);
            }
            cell_face_offsets.push_back((int)cell_to_face.size());
        }
    }
}
