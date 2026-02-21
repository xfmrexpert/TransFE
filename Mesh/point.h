/***************************************************************************
 *   Copyright (C) 2005-2024 by T. C. Raymond                              *
 *   tcraymond@inductivereasoning.com                                      *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                  *
 *                                                                         *
 ***************************************************************************/
 
#pragma once

#include <stdexcept>
#include <ostream>
#include <cstdint>
#include <cassert>

namespace TFEM
{
    /// The Point class is a container for 3D (and 2D) spatial coordinates.
    /// Renamed to PascalCase per project guidelines.
    class Point {
    public:
        /// Coordinates are public for direct access (POD-like behavior for math structs is common).
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;

        /// Default constructor initializes to (0.0, 0.0, 0.0).
        constexpr Point() = default;

        /// Constructor for 2D points (z defaults to 0.0).
        constexpr Point(double in_x, double in_y) : x(in_x), y(in_y), z(0.0) {}

        /// Constructor for 3D points.
        constexpr Point(double in_x, double in_y, double in_z) : x(in_x), y(in_y), z(in_z) {}

        /// Access by index (Read-only)
        /// Replaces 'X(n)' with standard operator syntax.
        constexpr double operator[](size_t index) const {
            assert(index < 3 && "Index out of bounds");
            if (index == 0) return x;
            if (index == 1) return y;
            return z;
        }

        /// Access by index (Read/Write)
        /// Allows pt[0] = 5.0;
        constexpr double& operator[](size_t index) {
            assert(index < 3 && "Index out of bounds");
            if (index == 0) return x;
            if (index == 1) return y;
            return z;
        }

        // Keep legacy accessor if strictly needed for backward compatibility during refactor,
        // but prefer operator[].
        double coordinate(uint8_t n) const {
            return (*this)[n];
        }

        /// Equality operator.
        constexpr bool operator==(const Point& other) const = default;

        /// Output stream operator.
        friend std::ostream& operator<<(std::ostream& os, const Point& pt) {
            os << "(" << pt.x << ", " << pt.y << ", " << pt.z << ")";
            return os;
        }

        Point& operator+=(const Point& rhs) {
            x += rhs.x;
            y += rhs.y;
            z += rhs.z;
            return *this;
        }

        Point& operator-=(const Point& rhs) {
            x -= rhs.x;
            y -= rhs.y;
            z -= rhs.z;
            return *this;
        }

        friend Point operator+(Point lhs, const Point& rhs) {
            lhs += rhs;
            return lhs;
        }

        friend Point operator-(Point lhs, const Point& rhs) {
            lhs -= rhs;
            return lhs;
        }

        friend Point operator*(double lhs, const Point& rhs) {
            return Point(lhs * rhs.x, lhs * rhs.y, lhs * rhs.z);
        }
        
        // Allow multiplication on the right side too: pt * 2.0
        friend Point operator*(const Point& lhs, double rhs) {
            return Point(lhs.x * rhs, lhs.y * rhs, lhs.z * rhs);
        }
    };
}
