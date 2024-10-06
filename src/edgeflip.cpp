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

// Function to eliminate obtuse triangles (if possible)
void eliminate_obtuse_triangles(CDT& cdt) {
    bool flipped = true;

    // Loop until no flips are performed
    while (flipped) {
        flipped = false;

        // Iterate through all faces (triangles)
        for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
            Face_handle fh = fit;
            
            // Check if the triangle is obtuse
            if (is_obtuse(fh->vertex(0)->point(), fh->vertex(1)->point(), fh->vertex(2)->point())) {
                // Try to flip the edge opposite to the obtuse angle
                for (int i = 0; i < 3; ++i) {
                    if (cdt.is_flipable(fh, i)) {
                        cdt.flip(fh, i);
                        flipped = true; // Set flag to continue checking
                        break;
                    }
                }
            }
        }
    }
}