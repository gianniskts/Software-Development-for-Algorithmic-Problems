#include "../includes/edgeflip.h"

// Function to check if a triangle is obtuse
bool is_obtuse(const Point& A, const Point& B, const Point& C) {
    // Define the angles of the triangle (middle is the vertex)
    double a1 = angle(A, B, C);
    double a2 = angle(B, A, C);
    double a3 = angle(A, C, B);
    // Check if any angle is obtuse
    if (a1 == CGAL::OBTUSE || a2 == CGAL::OBTUSE || a3 == CGAL::OBTUSE) {
        return true;
    }
    return false;
}

int count_obtuse_triangles(CDT& cdt) {
    int obtuse_triangle_count = 0;

    for (Face_iterator fi = cdt.finite_faces_begin(); fi != cdt.finite_faces_end(); ++fi) {
        Triangle_2 triangle = cdt.triangle(fi);

        if (is_obtuse(triangle[0], triangle[1], triangle[2])) {
            obtuse_triangle_count++;
        }
    }

    return obtuse_triangle_count;
}

// Function to eliminate obtuse triangles
void eliminate_obtuse_triangles(CDT& cdt, Polygon_2& steiner_points) {

    int obtuse_triangles_begin = count_obtuse_triangles(cdt);
    std::cout<<obtuse_triangles_begin<<std::endl;
    int obtuse_triangles_end = 0;

    while (obtuse_triangles_begin > obtuse_triangles_end) {
        
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
                Triangle_2 triangle(p0, p1, p2);
                Point circumcenter = CGAL::circumcenter(triangle);

                // Check if circumcenter is inside the domain
                if (cdt.is_infinite(cdt.locate(circumcenter))) {
                    // If not, insert midpoint of longest edge
                    K::FT l01 = CGAL::squared_distance(p0, p1);
                    K::FT l12 = CGAL::squared_distance(p1, p2);
                    K::FT l20 = CGAL::squared_distance(p2, p0);

                    K::FT max_length = std::max({l01, l12, l20});
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

                // Re-calculate obtuse triangles
                obtuse_triangles_end = count_obtuse_triangles(cdt);
                std::cout<<obtuse_triangles_end<<std::endl;
                break; // Exit loop to re-triangulate
            }
        }
    }
}