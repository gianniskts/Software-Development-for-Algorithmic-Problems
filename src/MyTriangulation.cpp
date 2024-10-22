#include "../includes/MyTriangulation.h"

// Forward declaration
bool is_obtuse(const Point& A, const Point& B, const Point& C);

// Function to mark triangulation domain
void Triangulation::mark_domain() {
    CGAL::mark_domain_in_triangulation(cdt, in_domain);
}

// Function to check if triangle is inside the boundary
bool Triangulation::is_face_in_domain(CDT::Face_handle face) const {
    return get(in_domain, face);
}

// Function to count triangulation's obtuse triangles
int Triangulation::count_obtuse_triangles() {
    int obtuse_triangle_count = 0;

    // Iterate over finite faces
    for (auto fi = cdt.finite_faces_begin(); fi != cdt.finite_faces_end(); ++fi) {
        // Check if the face is inside the polygon (domain)
        if (is_face_in_domain(fi)) {
            // Check if the triangle is obtuse
            if (is_obtuse(fi->vertex(0)->point(), fi->vertex(1)->point(), fi->vertex(2)->point())) {
                obtuse_triangle_count++;
            }
        }
    }
    
    return obtuse_triangle_count;
}