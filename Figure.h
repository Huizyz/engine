#ifndef ENGINE_FIGURE_H
#define ENGINE_FIGURE_H

#include <list>

#include "vector3d.h"
#include "Face.h"
#include "Color.h"
#include "Line2D.h"

class Figure {
public:
    Figure();

    std::vector<Vector3D> points;
    std::vector<Face> faces;
    Color color;


};
typedef std::list<Figure> Figures3D;

Matrix scalefigure(const double scale);
Matrix rotateX(const double angle);
Matrix rotateY(const double angle);
Matrix rotateZ(const double angle);
Matrix translate(const Vector3D &vector);

Matrix Transformation(double scale, double ro_X, double ro_Y, double ro_Z, const Vector3D& V3D);
void applyTransformation(Figure &fig, const Matrix &m);

Matrix eyePointTrans(Vector3D &eyepoint);

std::tuple<double, double, double>toPolar(Vector3D &point);

Point2D doProjection(Vector3D &point, const double d);
Lines2D doProjection(const Figures3D &figure3D);

//platonic figures
Figure createCube();
Figure createTetrahedron();
Figure createOctahedron();
Figure createIcosahedron();
Figure createDodecahedron();

Figure createSphere(const double radius, const int n);
Figure createCone(const int n, const double h);
Vector3D calculateMidpoint(const Vector3D& p1, const Vector3D& p2);
#endif //ENGINE_FIGURE_H
