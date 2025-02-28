#pragma once
#include "declarations.h"

struct Point3d {
    double x, y, z;

    Point3d operator+(const Point3d& other) const { return { x + other.x, y + other.y, z + other.z }; }
    Point3d operator-(const Point3d& other) const { return { x - other.x, y - other.y, z - other.z }; }
    Point3d operator*(double scalar) const { return { x * scalar, y * scalar, z * scalar }; }
    Point3d operator/(double scalar) const { return { x / scalar, y / scalar, z / scalar }; }
};

// Dot product of two vectors
double dot(const Point3d& a, const Point3d& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Cross product of two vectors
Point3d cross(const Point3d& a, const Point3d& b) {
    return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
}

// Magnitude (length) of a vector
double magnitude(const Point3d& a) {
    return std::sqrt(dot(a, a));
}

// Distance between two points
double distance(const Point3d& a, const Point3d& b) {
    return magnitude(a - b);
}

// Normalize a vector (unit vector)
Point3d normalize(const Point3d& a) {
    return a / magnitude(a);
}

// Angle between two vectors
double angle_between(const Point3d& a, const Point3d& b) {
    return std::acos(dot(a, b) / (magnitude(a) * magnitude(b)));
}

// Projection of vector a onto vector b
Point3d projection(const Point3d& a, const Point3d& b) {
    return b * (dot(a, b) / dot(b, b));
}

// Check if three points are collinear
bool are_collinear(const Point3d& a, const Point3d& b, const Point3d& c) {
    return magnitude(cross(b - a, c - a)) < 1e-9;
}

// Check if four points are coplanar
bool are_coplanar(const Point3d& a, const Point3d& b, const Point3d& c, const Point3d& d) {
    return std::abs(dot(cross(b - a, c - a), d - a)) < 1e-9;
}

// Calculate the area of a triangle formed by three points
double triangle_area(const Point3d& a, const Point3d& b, const Point3d& c) {
    return 0.5 * magnitude(cross(b - a, c - a));
}

// Calculate the volume of a tetrahedron formed by four points
double tetrahedron_volume(const Point3d& a, const Point3d& b, const Point3d& c, const Point3d& d) {
    return std::abs(dot(a - d, cross(b - d, c - d))) / 6.0;
}
