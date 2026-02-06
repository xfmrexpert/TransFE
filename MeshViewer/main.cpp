#include "raylib.h"
#include "raymath.h"
#include "../Mesh/mesh.h"
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace TFEM;

// A simple palette of distinguishable colors for region coloring
static Color AttributeColor(int attribute) {
    static const Color palette[] = {
        { 31, 119, 180, 255},  // blue
        {255, 127,  14, 255},  // orange
        { 44, 160,  44, 255},  // green
        {214,  39,  40, 255},  // red
        {148, 103, 189, 255},  // purple
        {140,  86,  75, 255},  // brown
        {227, 119, 194, 255},  // pink
        {127, 127, 127, 255},  // gray
        {188, 189,  34, 255},  // olive
        { 23, 190, 207, 255},  // cyan
        {174, 199, 232, 255},  // light blue
        {255, 187, 120, 255},  // light orange
        {152, 223, 138, 255},  // light green
        {255, 152, 150, 255},  // light red
        {197, 176, 213, 255},  // light purple
        {196, 156, 148, 255},  // light brown
    };
    static const int paletteSize = sizeof(palette) / sizeof(palette[0]);
    int idx = attribute % paletteSize;
    if (idx < 0) idx += paletteSize;
    return palette[idx];
}

void DrawElementWireframe(const CellView& cell, Color color) {
    auto nodes = cell.Nodes();
    auto V = [&](int i) -> Vector3 {
        auto n = cell.Node(i);
        return Vector3{ (float)n.X(), (float)n.Y(), (float)n.Z() };
    };

    switch (cell.Type()) {
        case ElementType::Segment:
            DrawLine3D(V(0), V(1), color);
            break;
        case ElementType::Triangle:
            DrawLine3D(V(0), V(1), color);
            DrawLine3D(V(1), V(2), color);
            DrawLine3D(V(2), V(0), color);
            break;
        case ElementType::Quad:
            DrawLine3D(V(0), V(1), color);
            DrawLine3D(V(1), V(2), color);
            DrawLine3D(V(2), V(3), color);
            DrawLine3D(V(3), V(0), color);
            break;
        case ElementType::Tetrahedron:
            DrawLine3D(V(0), V(1), color);
            DrawLine3D(V(1), V(2), color);
            DrawLine3D(V(2), V(0), color);
            DrawLine3D(V(0), V(3), color);
            DrawLine3D(V(1), V(3), color);
            DrawLine3D(V(2), V(3), color);
            break;
        case ElementType::Hexahedron:
            DrawLine3D(V(0), V(1), color);
            DrawLine3D(V(1), V(2), color);
            DrawLine3D(V(2), V(3), color);
            DrawLine3D(V(3), V(0), color);
            DrawLine3D(V(4), V(5), color);
            DrawLine3D(V(5), V(6), color);
            DrawLine3D(V(6), V(7), color);
            DrawLine3D(V(7), V(4), color);
            DrawLine3D(V(0), V(4), color);
            DrawLine3D(V(1), V(5), color);
            DrawLine3D(V(2), V(6), color);
            DrawLine3D(V(3), V(7), color);
            break;
        default:
            break;
    }
}

int main(int argc, char* argv[]) {
    // Initialization
    const int screenWidth = 1280;
    const int screenHeight = 720;

    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(screenWidth, screenHeight, "TFEM Mesh Viewer");

    // Load Mesh
    TFEM::Mesh mesh;
    if (argc > 1) {
        try {
            std::cout << "Loading mesh: " << argv[1] << "..." << std::endl;
            mesh.ReadGmsh(argv[1]);
            mesh.BuildTopology();
            std::cout << "Loaded " << mesh.NumNodes() << " nodes and " << mesh.NumCells() << " elements." << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Error loading mesh: " << e.what() << std::endl;
        }
    } else {
        std::cout << "Usage: MeshViewer <mesh.msh>" << std::endl;
    }

    // Compute bounding box and auto-fit camera
    auto bb = mesh.BoundingBox();
    Vector3 center = {
        (float)(bb[0] + bb[3]) * 0.5f,
        (float)(bb[1] + bb[4]) * 0.5f,
        (float)(bb[2] + bb[5]) * 0.5f
    };
    float extentX = (float)(bb[3] - bb[0]);
    float extentY = (float)(bb[4] - bb[1]);
    float extentZ = (float)(bb[5] - bb[2]);
    float maxExtent = std::max({extentX, extentY, extentZ, 1.0f});

    std::cout << "BBox: (" << bb[0] << ", " << bb[1] << ", " << bb[2]
              << ") -> (" << bb[3] << ", " << bb[4] << ", " << bb[5] << ")" << std::endl;

    Camera3D camera = { 0 };
    camera.target = center;
    // Position camera looking straight down at the XY plane
    camera.position = Vector3{ center.x, center.y, maxExtent * 1.5f };
    camera.up = Vector3{ 0.0f, 1.0f, 0.0f };
    camera.fovy = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;

    float zoomLevel = maxExtent * 1.5f; // Distance from mesh plane

    SetTargetFPS(60);

    while (!WindowShouldClose()) {
        // --- Custom 2D-style camera controls ---

        // Middle-click or left-click drag to pan
        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT) || IsMouseButtonDown(MOUSE_BUTTON_MIDDLE)) {
            Vector2 delta = GetMouseDelta();
            // Scale pan speed relative to zoom level
            float panSpeed = zoomLevel * 0.002f;
            camera.target.x -= delta.x * panSpeed;
            camera.target.y += delta.y * panSpeed;
            camera.position.x -= delta.x * panSpeed;
            camera.position.y += delta.y * panSpeed;
        }

        // Scroll to zoom
        float wheel = GetMouseWheelMove();
        if (wheel != 0.0f) {
            zoomLevel *= (1.0f - wheel * 0.1f);
            if (zoomLevel < 1.0f) zoomLevel = 1.0f;
            camera.position.z = camera.target.z + zoomLevel;
        }

        BeginDrawing();
            ClearBackground(RAYWHITE);

            BeginMode3D(camera);
                // Draw Mesh colored by attribute
                for (auto cell : mesh.Cells()) {
                    Color c = AttributeColor(cell.Attribute());
                    DrawElementWireframe(cell, c);
                }
            EndMode3D();

            DrawText(TextFormat("Cells: %i", mesh.NumCells()), 10, 10, 20, DARKGRAY);
            DrawText("Drag: Pan / Scroll: Zoom", 10, 30, 20, GRAY);
            DrawFPS(10, 50);

        EndDrawing();
    }

    CloseWindow();
    return 0;
}
