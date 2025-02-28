#pragma once
#include "declarations.h"

struct Point {
    int x, y;
    Point(int x = 0, int y = 0) : x(x), y(y) {}
};

int cross_product(Point a, Point b) {
    return a.x * b.y - a.y * b.x;  // 2D cross product of a and b
}

int cross_product(Point o, Point a, Point b) {
    return cross_product(Point(a.x - o.x, a.y - o.y), Point(b.x - o.x, b.y - o.y));
}

bool compare(Point a, Point b) {
    return a.x < b.x || (a.x == b.x && a.y < b.y);
}


// Convex Hull (Graham's Scan) algorithm
// Time complexity: O(n log n)
// Returns the convex hull of a set of points
vector<Point> convex_hull(vector<Point>& points) {
    sort(points.begin(), points.end(), compare);

    vector<Point> hull;

    // Lower hull
    for (const Point& p : points) {
        while (hull.size() >= 2 && cross_product(hull[hull.size() - 2], hull.back(), p) <= 0)
            hull.pop_back();
        hull.push_back(p);
    }

    // Upper hull
    int lower_size = hull.size();
    for (int i = points.size() - 1; i >= 0; i--) {
        while (hull.size() > lower_size && cross_product(hull[hull.size() - 2], hull.back(), points[i]) <= 0)
            hull.pop_back();
        hull.push_back(points[i]);
    }

    hull.pop_back();  // Remove the last point because it's the same as the first one.
    return hull;
}

// Check if a point is inside a polygon
bool point_in_polygon(const Point& pt, const vector<Point>& polygon) {
    int n = polygon.size();
    bool inside = false;
    for (int i = 0, j = n - 1; i < n; j = i++) {
        // Check if pt is inside the polygon using the ray-casting method
        if ((polygon[i].y > pt.y) != (polygon[j].y > pt.y) &&
            pt.x < (polygon[j].x - polygon[i].x) * (pt.y - polygon[i].y) / (polygon[j].y - polygon[i].y) + polygon[i].x) {
            inside = !inside;
        }
    }
    return inside;
}

double distance(Point a, Point b) {
    return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

double square_distance(Point a, Point b) {
	return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
}

// Returns 0 if p, q, r are collinear
int orientation(Point p, Point q, Point r) {
    int val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
    if (val == 0) return 0;  // collinear
    return (val > 0) ? 1 : 2;  // 1 for clockwise, 2 for counterclockwise
}

double area_of_triangle(Point a, Point b, Point c) {
    return abs(a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y)) / 2.0;
}

// Check if point q lies on line segment 'pr'
bool on_segment(Point p, Point q, Point r) {
    return r.x <= max(p.x, q.x) && r.x >= min(p.x, q.x) && r.y <= max(p.y, q.y) && r.y >= min(p.y, q.y);
}

// Check if two line segments intersect
bool do_intersect(Point p1, Point q1, Point p2, Point q2) {
    // Find the four orientations needed for general and special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    // General case
    if (o1 != o2 && o3 != o4)
        return true;

    // Special cases (collinear points)
    if (o1 == 0 && on_segment(p1, q1, p2)) return true;
    if (o2 == 0 && on_segment(p1, q1, q2)) return true;
    if (o3 == 0 && on_segment(p2, q2, p1)) return true;
    if (o4 == 0 && on_segment(p2, q2, q1)) return true;

    return false;  // Doesn't intersect
}



// Check if point q lies on the line segment 'pr' and not on the endpoints
double convex_hull_area(const vector<Point>& hull) {
    int n = hull.size();
    double area = 0;
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        area += hull[i].x * hull[j].y;
        area -= hull[j].x * hull[i].y;
    }
    return abs(area) / 2.0;
}

double point_line_distance(Point p, Point a, Point b) {
    return abs((b.y - a.y) * p.x - (b.x - a.x) * p.y + b.x * a.y - b.y * a.x) /
        sqrt((b.y - a.y) * (b.y - a.y) + (b.x - a.x) * (b.x - a.x));
}

Point project_point_on_line(Point p, Point a, Point b) {
    double r = ((p.x - a.x) * (b.x - a.x) + (p.y - a.y) * (b.y - a.y)) /
        (pow(b.x - a.x, 2) + pow(b.y - a.y, 2));
    return Point(a.x + r * (b.x - a.x), a.y + r * (b.y - a.y));
}

double angle_between_vectors(Point a, Point b) {
    double dot = a.x * b.x + a.y * b.y;
    double len_a = sqrt(a.x * a.x + a.y * a.y);
    double len_b = sqrt(b.x * b.x + b.y * b.y);
    return acos(dot / (len_a * len_b));
}


struct Line {
    long long m, b;
    Line(long long m, long long b) : m(m), b(b) {}
    long long evaluate(long long x) {
        return m * x + b;
    }
};

// Convex Hull Trick: Efficient line maximization / minimization
// Time complexity: O(n) for adding lines, O(log n) for querying
/*
This is a dynamic programming technique used when you have a set of 
lines in the form y = mx + b, and you want to efficiently compute 
the maximum or minimum y-value for a given x. This is typically 
used in problems like optimization with linear functions.
*/
struct ConvexHullTrick {
    vector<Line> lines;
    bool is_bad(const Line& l1, const Line& l2, const Line& l3) {
        return (l3.b - l1.b) * (l1.m - l2.m) < (l2.b - l1.b) * (l1.m - l3.m);
    }
    void add_line(long long m, long long b) {
        Line new_line(m, b);
        while (lines.size() >= 2 && is_bad(lines[lines.size() - 2], lines[lines.size() - 1], new_line))
            lines.pop_back();
        lines.push_back(new_line);
    }
    long long get_max(long long x) {
        while (lines.size() >= 2 && lines[1].evaluate(x) <= lines[0].evaluate(x))
            lines.erase(lines.begin());
        return lines[0].evaluate(x);
    }
};

// Triangulate polygon (useful for polygon area and other operations)
vector<tuple<Point, Point, Point>> triangulate_polygon(const vector<Point>& polygon) {
    vector<tuple<Point, Point, Point>> triangles;
    int n = polygon.size();
    for (int i = 1; i < n - 1; ++i) {
        triangles.push_back({ polygon[0], polygon[i], polygon[i + 1] });
    }
    return triangles;
}

// Check if a point lies on the infinite line passing through two points
bool is_point_on_line(Point p, Point a, Point b) {
    return (b.y - a.y) * (p.x - a.x) == (b.x - a.x) * (p.y - a.y);
}

// Check if point p lies inside the convex polygon
bool point_in_convex_polygon(const Point& p, const vector<Point>& polygon) {
    int n = polygon.size();
    int left = 1, right = n - 1;
    while (right - left > 1) {
        int mid = (left + right) / 2;
        if (cross_product(polygon[0], polygon[mid], p) < 0)
            right = mid;
        else
            left = mid;
    }
    return cross_product(polygon[0], polygon[left], p) >= 0 && cross_product(polygon[left], polygon[right], p) >= 0;
}

// Returns true if point p lies inside the triangle formed by a, b, c
bool point_in_triangle(Point p, Point a, Point b, Point c) {
    double areaABC = area_of_triangle(a, b, c);
    double areaPAB = area_of_triangle(p, a, b);
    double areaPBC = area_of_triangle(p, b, c);
    double areaPCA = area_of_triangle(p, c, a);
    return (areaABC == areaPAB + areaPBC + areaPCA);
}

// Approximate geometric median (Weiszfeld's algorithm)
Point geometric_median(const vector<Point>& points, int max_iterations = 100, double epsilon = 1e-6) {
    double x = 0, y = 0;
    int n = points.size();
    for (auto& pt : points) {
        x += pt.x;
        y += pt.y;
    }
    x /= n;
    y /= n;

    for (int iter = 0; iter < max_iterations; ++iter) {
        double num_x = 0, num_y = 0, denom = 0;
        for (auto& pt : points) {
            double dist = sqrt((pt.x - x) * (pt.x - x) + (pt.y - y) * (pt.y - y));
            if (dist > 0) {
                num_x += pt.x / dist;
                num_y += pt.y / dist;
                denom += 1.0 / dist;
            }
        }
        double new_x = num_x / denom;
        double new_y = num_y / denom;
        if (abs(new_x - x) < epsilon && abs(new_y - y) < epsilon)
            break;
        x = new_x;
        y = new_y;
    }
    return Point(x, y);
}

// Returns the minimum distance between two convex polygons
double min_distance_between_polygons(const vector<Point>& poly1, const vector<Point>& poly2) {
    double min_dist = numeric_limits<double>::infinity();
    int n = poly1.size(), m = poly2.size();

    // For each edge of poly1, find the minimum distance to poly2
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            min_dist = min(min_dist, point_line_distance(poly1[i], poly1[(i + 1) % n], poly2[j]));
        }
    }
    return min_dist;
}

// Rotate point p around point o by angle theta (in radians)
Point rotate_point(Point p, Point o, double theta) {
    double cos_theta = cos(theta);
    double sin_theta = sin(theta);
    double dx = p.x - o.x;
    double dy = p.y - o.y;
    return Point(o.x + dx * cos_theta - dy * sin_theta, o.y + dx * sin_theta + dy * cos_theta);
}

