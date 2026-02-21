#include "mesh.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <stdexcept>

namespace TFEM {

    // Maps a Gmsh element type integer to our ElementType enum
    static ElementType MapGmshType(int t) {
        switch(t) {
            case 15: return ElementType::Point;
            case 1:  return ElementType::Segment;
            case 8:  return ElementType::Segment;       // 3-node second-order line
            case 2:  return ElementType::Triangle;
            case 9:  return ElementType::Triangle;       // 6-node second-order triangle
            case 3:  return ElementType::Quad;
            case 10: return ElementType::Quad;           // 9-node second-order quad
            case 4:  return ElementType::Tetrahedron;
            case 11: return ElementType::Tetrahedron;    // 10-node second-order tet
            case 5:  return ElementType::Hexahedron;
            case 12: return ElementType::Hexahedron;     // 27-node second-order hex
            case 6:  return ElementType::Prism;
            case 7:  return ElementType::Pyramid;
            default: throw std::runtime_error("Unsupported Gmsh element type: " + std::to_string(t));
        }
    }

    // Total number of nodes for a given Gmsh element type (including mid-edge nodes)
    static int GetGmshNodeCount(int t) {
        switch(t) {
            case 15: return 1;   // Point
            case 1:  return 2;   // 2-node line
            case 8:  return 3;   // 3-node second-order line
            case 2:  return 3;   // 3-node triangle
            case 9:  return 6;   // 6-node second-order triangle
            case 3:  return 4;   // 4-node quad
            case 10: return 9;   // 9-node second-order quad
            case 4:  return 4;   // 4-node tetrahedron
            case 11: return 10;  // 10-node second-order tet
            case 5:  return 8;   // 8-node hexahedron
            case 12: return 27;  // 27-node second-order hex
            case 6:  return 6;   // 6-node prism
            case 7:  return 5;   // 5-node pyramid
            default: return 0;
        }       
    }

    // Number of corner (vertex) nodes for a given Gmsh element type.
    // For first-order types this equals GetGmshNodeCount.
    // For second-order types this is the number of corner nodes only.
    static int GetGmshCornerNodeCount(int t) {
        switch(t) {
            case 15: return 1;
            case 1:  return 2;
            case 8:  return 2;   // 3-node line -> 2 corners
            case 2:  return 3;
            case 9:  return 3;   // 6-node tri -> 3 corners
            case 3:  return 4;
            case 10: return 4;   // 9-node quad -> 4 corners
            case 4:  return 4;
            case 11: return 4;   // 10-node tet -> 4 corners
            case 5:  return 8;
            case 12: return 8;   // 27-node hex -> 8 corners
            case 6:  return 6;
            case 7:  return 5;
            default: return 0;
        }
    }

    // -------------------------------------------------------------------------
    // MSH 2.x Reader
    // -------------------------------------------------------------------------
    // Format reference: https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
    //
    // $Nodes
    //   number-of-nodes
    //   node-number  x  y  z
    //   ...
    // $EndNodes
    //
    // $Elements
    //   number-of-elements
    //   elm-number  elm-type  number-of-tags  <tags...>  node-numbers...
    //   ...
    // $EndElements

    void Mesh::ReadGmsh2(std::ifstream& file) {
        std::string line;

        // Tag-to-physical map from $PhysicalNames (physical_tag -> name)
        // In MSH 2.x, the first tag of an element is the physical group tag.
        std::unordered_map<int, int> entityToPhysical; // entityTag -> physicalTag

        while (std::getline(file, line)) {
            if (line == "$PhysicalNames") {
                int numNames;
                file >> numNames;
                // We read but don't strictly need the names for now.
                // Format: dimension  physical-tag  "name"
                for (int i = 0; i < numNames; ++i) {
                    int dim, tag;
                    std::string name;
                    file >> dim >> tag >> name;
                    // name includes quotes, strip them if needed
                }
                std::getline(file, line); // consume newline
                std::getline(file, line); // $EndPhysicalNames
            }
            else if (line == "$Nodes") {
                int numNodes;
                file >> numNodes;

                coordinates.resize(numNodes * 3);

                // In MSH 2.x, node tags are 1-based and typically contiguous.
                // We build a tag-to-index map for safety.
                std::vector<int> nodeTagToIdx;

                for (int i = 0; i < numNodes; ++i) {
                    int tag;
                    double x, y, z;
                    file >> tag >> x >> y >> z;

                    int idx = i; // 0-based sequential index
                    coordinates[idx * 3 + 0] = x;
                    coordinates[idx * 3 + 1] = y;
                    coordinates[idx * 3 + 2] = z;

                    if (static_cast<size_t>(tag) >= nodeTagToIdx.size()) {
                        nodeTagToIdx.resize(tag + 1, -1);
                    }
                    nodeTagToIdx[tag] = idx;
                }
                std::getline(file, line); // consume newline
                std::getline(file, line); // $EndNodes

                // Now read $Elements (must come after $Nodes in our parsing)
                // Continue to next section
                while (std::getline(file, line)) {
                    if (line == "$Elements") {
                        int numElements;
                        file >> numElements;

                        cell_node_offsets.reserve(numElements + 1);
                        cell_node_offsets.push_back(0);
                        cell_types.reserve(numElements);
                        cell_attributes.reserve(numElements);

                        for (int i = 0; i < numElements; ++i) {
                            int elmNumber, elmType, numTags;
                            file >> elmNumber >> elmType >> numTags;

                            // Read tags. First tag = physical entity, second = elementary entity
                            int physicalTag = 0;
                            for (int t = 0; t < numTags; ++t) {
                                int tag;
                                file >> tag;
                                if (t == 0) physicalTag = tag;
                            }

                            int totalNodes = GetGmshNodeCount(elmType);
                            int cornerNodes = GetGmshCornerNodeCount(elmType);

                            // Read all node tags
                            std::vector<int> allNodes(totalNodes);
                            for (int n = 0; n < totalNodes; ++n) {
                                file >> allNodes[n];
                            }

                            // Store only corner nodes (first cornerNodes entries)
                            for (int n = 0; n < cornerNodes; ++n) {
                                cell_to_node.push_back(nodeTagToIdx[allNodes[n]]);
                            }

                            cell_node_offsets.push_back(static_cast<int>(cell_to_node.size()));
                            cell_attributes.push_back(physicalTag);
                            cell_types.push_back(MapGmshType(elmType));
                        }

                        // Initialize empty topology offsets
                        cell_edge_offsets.assign(NumCells() + 1, 0);
                        cell_face_offsets.assign(NumCells() + 1, 0);

                        std::getline(file, line); // consume newline
                        std::getline(file, line); // $EndElements
                        break; // Done
                    }
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // MSH 4.x Reader
    // -------------------------------------------------------------------------

    void Mesh::ReadGmsh4(std::ifstream& file) {
        std::string line;

        // Map from GMsh node tag to 0-based index in coordinates
        std::vector<int> nodeTagToIdx;

        while (std::getline(file, line)) {
            if (line == "$Nodes") {
                int numEntityBlocks, numNodes, minNodeTag, maxNodeTag;
                file >> numEntityBlocks >> numNodes >> minNodeTag >> maxNodeTag;
                
                coordinates.resize(numNodes * 3);
                if (static_cast<size_t>(maxNodeTag) >= nodeTagToIdx.size()) {
                    nodeTagToIdx.resize(maxNodeTag + 1, -1);
                }

                int currentNodeIdx = 0;

                for (int i = 0; i < numEntityBlocks; ++i) {
                    int entityDim, entityTag, parametric, numNodesInBlock;
                    file >> entityDim >> entityTag >> parametric >> numNodesInBlock;

                    // Read Tags
                    std::vector<int> blockTags(numNodesInBlock);
                    for (int n = 0; n < numNodesInBlock; ++n) {
                        file >> blockTags[n];
                    }

                    // Read Coordinates
                    for (int n = 0; n < numNodesInBlock; ++n) {
                        double x, y, z;
                        file >> x >> y >> z;

                        int idx = currentNodeIdx++;
                        
                        coordinates[idx * 3 + 0] = x;
                        coordinates[idx * 3 + 1] = y;
                        coordinates[idx * 3 + 2] = z;

                        if (static_cast<size_t>(blockTags[n]) >= nodeTagToIdx.size()) {
                            nodeTagToIdx.resize(blockTags[n] + 1, -1);
                        }
                        nodeTagToIdx[blockTags[n]] = idx;
                    }
                }
                std::getline(file, line); // consume newline
                std::getline(file, line); // $EndNodes
            }
            else if (line == "$Elements") {
                int numEntityBlocks, numElements, minEleTag, maxEleTag;
                file >> numEntityBlocks >> numElements >> minEleTag >> maxEleTag;
                
                cell_node_offsets.reserve(numElements + 1);
                cell_node_offsets.push_back(0);
                
                cell_types.reserve(numElements);
                cell_attributes.reserve(numElements);

                for (int i = 0; i < numEntityBlocks; ++i) {
                    int entityDim, entityTag, elementType, numElementsInBlock;
                    file >> entityDim >> entityTag >> elementType >> numElementsInBlock;
                    
                    int nNodes = GetGmshNodeCount(elementType);
                    int cornerNodes = GetGmshCornerNodeCount(elementType);
                    ElementType tfemType = MapGmshType(elementType);

                    for (int k = 0; k < numElementsInBlock; ++k) {
                         size_t eleTag;
                         file >> eleTag;
                         
                         std::vector<int> allNodes(nNodes);
                         for (int n = 0; n < nNodes; ++n) {
                             file >> allNodes[n];
                         }

                         // Store only corner nodes
                         for (int n = 0; n < cornerNodes; ++n) {
                             cell_to_node.push_back(nodeTagToIdx[allNodes[n]]);
                         }
                         
                         cell_node_offsets.push_back(static_cast<int>(cell_to_node.size()));
                         cell_attributes.push_back(entityTag); 
                         cell_types.push_back(tfemType);
                    }
                }

                // Initialize empty topology offsets
                cell_edge_offsets.assign(NumCells() + 1, 0);
                cell_face_offsets.assign(NumCells() + 1, 0);
                
                std::getline(file, line); // consume rest of line
                std::getline(file, line); // $EndElements
            }
        }
    }

    // -------------------------------------------------------------------------
    // Public Entry Point
    // -------------------------------------------------------------------------

    void Mesh::ReadGmsh(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Could not open file: " + filename);
        }

        Clear();

        std::string line;
        double version = 0.0;
        int fileType = 0;
        int dataSize = 0;

        // Read format header
        std::getline(file, line);
        if (line != "$MeshFormat") {
            throw std::runtime_error("Invalid Gmsh file: missing $MeshFormat header.");
        }
        file >> version >> fileType >> dataSize;
        std::getline(file, line); // consume newline
        std::getline(file, line); // $EndMeshFormat

        if (fileType != 0) {
            throw std::runtime_error("Only ASCII Gmsh files are supported (got file-type " + std::to_string(fileType) + ").");
        }

        if (version >= 4.0) {
            ReadGmsh4(file);
        } else if (version >= 2.0) {
            ReadGmsh2(file);
        } else {
            throw std::runtime_error("Unsupported Gmsh format version: " + std::to_string(version));
        }

        // Determine mesh dimension based on element types
        int max_dim = 0;
        for (auto type : cell_types) {
            int d = 0;
            switch(type) {
                case ElementType::Point: d = 0; break;
                case ElementType::Segment: d = 1; break;
                case ElementType::Triangle: 
                case ElementType::Quad: d = 2; break;
                case ElementType::Tetrahedron:
                case ElementType::Hexahedron:
                case ElementType::Prism:
                case ElementType::Pyramid: d = 3; break;
            }
            if (d > max_dim) max_dim = d;
        }
        dim = max_dim; 

        // Determine spatial dimension based on nodal coordinates
        // Identify if Y or Z are ever non-zero
        bool has_y = false;
        bool has_z = false;
        
        // Coordinates are stored as x0, y0, z0, x1, y1, z1 ...
        for (size_t i = 0; i < coordinates.size(); i += 3) {
            if (std::abs(coordinates[i+1]) > 1e-12) has_y = true;
            if (std::abs(coordinates[i+2]) > 1e-12) has_z = true;
            if (has_y && has_z) break; 
        }

        if (has_z) space_dim = 3;
        else if (has_y) space_dim = 2;
        else space_dim = 1;

        // Ensure space_dim is at least dim 
        // (e.g., a planar mesh in X-Y plane is dim=2, space_dim=2.
        //  a surface mesh in 3D might be dim=2, space_dim=3)
        if (space_dim < dim) space_dim = dim;
    }
}
