#include "../includes/edgeflip.h"

// Function to check if a triangle is obtuse
bool is_obtuse(const Point& A, const Point& B, const Point& C) {
    // Compute vectors
    K::Vector_2 vAB = B - A;
    K::Vector_2 vAC = C - A;

    // Compute dot product
    double dot_product = vAB * vAC;

    // If dot product is negative, angle at A is obtuse
    if (dot_product < 0) return true;

    // Repeat for other vertices
    vAB = A - B;
    vAC = C - B;
    dot_product = vAB * vAC;
    if (dot_product < 0) return true;

    vAB = A - C;
    vAC = B - C;
    dot_product = vAB * vAC;
    if (dot_product < 0) return true;

    return false;
}

// Function to eliminate obtuse triangles
void eliminate_obtuse_triangles(CDT& cdt, Polygon_2& steiner_points) {
    bool changes = true;
    int iteration = 0;
    const int max_iterations = 1000; // Prevent infinite loops

    while (changes && iteration < max_iterations) {
        changes = false;
        iteration++;

        // Edge flipping to improve triangle quality
        CGAL::make_conforming_Delaunay_2(cdt);
        CGAL::make_conforming_Gabriel_2(cdt);

        // Check for obtuse triangles
        for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
            Face_handle fh = fit;
            Point p0 = fh->vertex(0)->point();
            Point p1 = fh->vertex(1)->point();
            Point p2 = fh->vertex(2)->point();

            if (is_obtuse(p0, p1, p2)) {
                // Find the vertex opposite the obtuse angle
                Point obtuse_vertex;
                if (angle(p0, p1, p2) > 90) {
                    obtuse_vertex = p0;
                } else if (angle(p1, p2, p0) > 90) {
                    obtuse_vertex = p1;
                } else {
                    obtuse_vertex = p2;
                }

                // Insert Steiner point at the circumcenter of the triangle
                CGAL::Triangle_2<K> triangle(p0, p1, p2);
                Point circumcenter = CGAL::circumcenter(triangle);

                // Check if circumcenter is inside the domain
                if (cdt.is_infinite(cdt.locate(circumcenter))) {
                    // If not, insert midpoint of longest edge
                    double l01 = CGAL::squared_distance(p0, p1);
                    double l12 = CGAL::squared_distance(p1, p2);
                    double l20 = CGAL::squared_distance(p2, p0);

                    double max_length = std::max({l01, l12, l20});
                    Point midpoint;

                    if (max_length == l01) {
                        midpoint = CGAL::midpoint(p0, p1);
                    } else if (max_length == l12) {
                        midpoint = CGAL::midpoint(p1, p2);
                    } else {
                        midpoint = CGAL::midpoint(p2, p0);
                    }

                    cdt.insert(midpoint);
                    steiner_points.push_back(midpoint);
                } else {
                    // Insert circumcenter as Steiner point
                    cdt.insert(circumcenter);
                    steiner_points.push_back(circumcenter);
                }

                changes = true;
                break; // Exit loop to re-triangulate
            }
        }
    }
}