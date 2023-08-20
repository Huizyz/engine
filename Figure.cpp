#define _USE_MATH_DEFINES
#include "Figure.h"
#include "Line2D.h"
#include <cmath>
#include <tuple>




Matrix scalefigure(const double scale){
Matrix scaled;
    scaled(1,1) = scale;
    scaled(1,2) = 0;
    scaled(1,3) = 0;
    scaled(1,4) = 0;

    scaled(2, 1) = 0;
    scaled(2, 2) = scale;
    scaled(2, 3) = 0;
    scaled(2, 4) = 0;

    scaled(3, 1) = 0;
    scaled(3, 2) = 0;
    scaled(3, 3) = scale;
    scaled(3, 4) = 0;

    scaled(4, 1) = 0;
    scaled(4, 2) = 0;
    scaled(4, 3) = 0;
    scaled(4, 4) = 1;
    return scaled;
}

Matrix rotateX(const double angle){
Matrix rotationX;
    rotationX(1,1) = 1;
    rotationX(1,2) = 0;
    rotationX(1,3) = 0;
    rotationX(1,4) = 0;

    rotationX(2,1) = 0;
    rotationX(2,2) = cos(angle * (M_PI / 180));
    rotationX(2,3) = sin(angle * (M_PI / 180));
    rotationX(2,4) = 0;

    rotationX(3,1) = 0;
    rotationX(3,2) = -sin(angle * (M_PI / 180));
    rotationX(3,3) = cos(angle * (M_PI / 180));
    rotationX(3,4) = 0;

    rotationX(4,1) = 0;
    rotationX(4,2) = 0;
    rotationX(4,3) = 0;
    rotationX(4,4) = 1;
    return rotationX;
}

Matrix rotateY(const double angle){
Matrix rotationY;
    rotationY(1,1) = cos(angle * (M_PI / 180));
    rotationY(1,2) = 0;
    rotationY(1,3) = -sin(angle * (M_PI / 180));
    rotationY(1,4) = 0;

    rotationY(2,1) = 0;
    rotationY(2,2) = 1;
    rotationY(2,3) = 0;
    rotationY(2,4) = 0;

    rotationY(3,1) = sin(angle * (M_PI / 180));
    rotationY(3,2) = 0;
    rotationY(3,3) = cos(angle * (M_PI / 180));
    rotationY(3,4) = 0;

    rotationY(4,1) = 0;
    rotationY(4,2) = 0;
    rotationY(4,3) = 0;
    rotationY(4,4) = 1;
    return rotationY;
}

Matrix rotateZ(const double angle){
Matrix rotationZ;
    rotationZ(1,1) = cos(angle * (M_PI / 180));
    rotationZ(1,2) = sin(angle * (M_PI / 180));
    rotationZ(1,3) = 0;
    rotationZ(1,4) = 0;

    rotationZ(2,1) = -sin(angle * (M_PI / 180));
    rotationZ(2,2) = cos(angle * (M_PI / 180));
    rotationZ(2,3) = 0;
    rotationZ(2,4) = 0;

    rotationZ(3,1) = 0;
    rotationZ(3,2) = 0;
    rotationZ(3,3) = 1;
    rotationZ(3,4) = 0;

    rotationZ(4,1) = 0;
    rotationZ(4,2) = 0;
    rotationZ(4,3) = 0;
    rotationZ(4,4) = 1;
    return rotationZ;
}

Matrix translate(const Vector3D &vector) {
    Matrix trans;
    trans(1,1) = 1;
    trans(1,2) = 0;
    trans(1,3) = 0;
    trans(1,4) = 0;

    trans(2, 1) = 0;
    trans(2, 2) = 1;
    trans(2, 3) = 0;
    trans(2, 4) = 0;

    trans(3, 1) = 0;
    trans(3, 2) = 0;
    trans(3, 3) = 1;
    trans(3, 4) = 0;

    trans(4, 1) = vector.x;
    trans(4, 2) = vector.y;
    trans(4, 3) = vector.z;
    trans(4, 4) = 1;
    return trans;
}

Matrix Transformation(double scale, double ro_X, double ro_Y, double ro_Z, const Vector3D& V3D){
    Matrix transformation;
    transformation = scalefigure(scale) * rotateX(ro_X) * rotateY(ro_Y) * rotateZ(ro_Z) * translate(V3D);
    return transformation;
}

void applyTransformation(Figure &fig, const Matrix &m){
    for (auto &elements : fig.points) {
        elements *= m;
    }
}

std::tuple<double, double, double>toPolar(Vector3D &point){
    double r = sqrt(pow(point.x, 2) + pow(point.y, 2) + pow(point.z, 2));
    double theta = atan2(point.y, point.x);
    double phi = acos(point.z / r);

    return {theta, phi, r};
}

Matrix eyePointTrans(Vector3D &eyepoint){
    std::tuple<double, double, double> Bolcoordinaat = toPolar(eyepoint);

    Matrix eyePointMatrix;
    eyePointMatrix(1, 1) = -sin(std::get<0>(Bolcoordinaat));
    eyePointMatrix(1, 2) = -(cos(std::get<0>(Bolcoordinaat)) * cos(std::get<1>(Bolcoordinaat)));
    eyePointMatrix(1, 3) = cos(std::get<0>(Bolcoordinaat)) * sin(std::get<1>(Bolcoordinaat));
    eyePointMatrix(1,4) = 0;

    eyePointMatrix(2, 1) = cos(std::get<0>(Bolcoordinaat));
    eyePointMatrix(2, 2) = -(sin(std::get<0>(Bolcoordinaat)) * cos(std::get<1>(Bolcoordinaat)));
    eyePointMatrix(2, 3) = sin(std::get<0>(Bolcoordinaat)) * sin(std::get<1>(Bolcoordinaat));
    eyePointMatrix(2,4) = 0;

    eyePointMatrix(3,1) = 0;
    eyePointMatrix(3, 2) = sin(std::get<1>(Bolcoordinaat));
    eyePointMatrix(3, 3) = cos(std::get<1>(Bolcoordinaat));
    eyePointMatrix(3,4) = 0;

    eyePointMatrix(4,1) = 0;
    eyePointMatrix(4,2) = 0;
    eyePointMatrix(4, 3) = -std::get<2>(Bolcoordinaat);
    eyePointMatrix(4, 4) = 1;
    return eyePointMatrix;
}

Point2D doProjection(Vector3D &point, const double d){
    double xACC = (d * point.x) / (-point.z);
    double yACC = (d * point.y) / (-point.z);
    return Point2D(xACC,yACC);
}

Lines2D doProjection(const Figures3D &figures){
    Lines2D lines;
    for (auto figure:figures) {
        for (auto face:figure.faces) {
            for (int i = 0; i < face.point_indexes.size(); ++i) {
                Point2D p1 = doProjection(figure.points[face.point_indexes[i]],1);
                Point2D p2 = doProjection(figure.points[face.point_indexes[(i+1)%face.point_indexes.size()]],1);
                lines.emplace_back(Line2D(p1,p2,figure.color));

            }
        }
    }
    return lines;
}

Figure createCube(){
    Figure cube;
    //points
    Vector3D p1 = Vector3D::point(1,-1,-1);
    Vector3D p2 = Vector3D::point(-1,1,-1);
    Vector3D p3 = Vector3D::point(1,1,1);
    Vector3D p4 = Vector3D::point(-1,-1,1);
    Vector3D p5 = Vector3D::point(1,1,-1);
    Vector3D p6 = Vector3D::point(-1,-1,-1);
    Vector3D p7 = Vector3D::point(1,-1,1);
    Vector3D p8 = Vector3D::point(-1,1,1);
    std::vector<Vector3D> points = {p1, p2, p3, p4, p5, p6, p7, p8};
    for (const auto& point : points) {
        cube.points.push_back(point);
    }

    //faces
    vector<vector<int>> faceIndices = {
            {1, 5, 3, 7}, {5, 2, 8, 3}, {2, 6, 4, 8},
            {6, 1, 7, 4}, {7, 3, 8, 4}, {1, 6, 2, 5}
    };

    for (const auto& indices : faceIndices) {
        cube.faces.emplace_back(indices);
    }


    return cube;
}

Figure createTetrahedron(){
    Figure tetrahedron;
    Vector3D p1 = Vector3D::point(1,-1,-1);
    Vector3D p2 = Vector3D::point(-1,1,-1);
    Vector3D p3 = Vector3D::point(1,1,1);
    Vector3D p4 = Vector3D::point(-1,-1,1);
    std::vector<Vector3D> points = {p1, p2, p3, p4};
    for (const auto& point : points) {
        tetrahedron.points.push_back(point);
    }

    vector<vector<int>> faceIndices = {
            {1, 2, 3}, {2, 4, 3}, {1, 4, 2}, {1, 3, 4}
    };

    for (const auto& indices : faceIndices) {
        tetrahedron.faces.emplace_back(indices);
    }


    return tetrahedron;
}
Figure createOctahedron(){
    Figure octahedron;
    Vector3D p1 = Vector3D::point(1, 0, 0);
    Vector3D p2 = Vector3D::point(0,1,0);
    Vector3D p3 = Vector3D::point(-1,0,0);
    Vector3D p4 = Vector3D::point(0,-1,0);
    Vector3D p5 = Vector3D::point(0,0,-1);
    Vector3D p6 = Vector3D::point(0,0,1);
    vector<Vector3D> points = {p1, p2, p3, p4, p5, p6};
    for (const auto& point : points) {
        octahedron.points.push_back(point);
    }

    vector<vector<int>> faceIndices = {
            {1, 2, 6}, {2, 3, 6}, {3, 4, 6}, {4, 1, 6},
            {2, 1, 5}, {3, 2, 5}, {4, 3, 5}, {1, 4, 5}
    };

    for (const auto& indices : faceIndices) {
        octahedron.faces.emplace_back(indices);
    }


    return octahedron;
}
Figure createIcosahedron(){
    Figure icosahedron;

    //points
    Vector3D p1 = Vector3D::point(0, 0, (sqrt(5)/2));
    Vector3D p2 = Vector3D::point(cos((2 - 2) * 2 * M_PI / 5), sin((2 - 2) * 2 * M_PI / 5), 0.5);
    Vector3D p3 = Vector3D::point(cos((3 - 2) * 2 * M_PI / 5), sin((3 - 2) * 2 * M_PI / 5), 0.5);
    Vector3D p4 = Vector3D::point(cos((4 - 2) * 2 * M_PI / 5), sin((4 - 2) * 2 * M_PI / 5), 0.5);
    Vector3D p5 = Vector3D::point(cos((5 - 2) * 2 * M_PI / 5), sin((5 - 2) * 2 * M_PI / 5), 0.5);
    Vector3D p6 = Vector3D::point(cos((6 - 2) * 2 * M_PI / 5), sin((6 - 2) * 2 * M_PI / 5), 0.5);
    Vector3D p7 = Vector3D::point(cos(M_PI / 5 + (7 - 7) * 2 * M_PI / 5), sin(M_PI / 5 + (7 - 7) * 2 * M_PI / 5),-0.5);
    Vector3D p8 = Vector3D::point(cos(M_PI / 5 + (8 - 7) * 2 * M_PI / 5), sin(M_PI / 5 + (8 - 7) * 2 * M_PI / 5),-0.5);
    Vector3D p9 = Vector3D::point(cos(M_PI / 5 + (9 - 7) * 2 * M_PI / 5), sin(M_PI / 5 + (9 - 7) * 2 * M_PI / 5),-0.5);
    Vector3D p10 = Vector3D::point(cos(M_PI / 5 + (10 - 7) * 2 * M_PI / 5), sin(M_PI / 5 + (10 - 7) * 2 * M_PI / 5),-0.5);
    Vector3D p11 = Vector3D::point(cos(M_PI / 5 + (11 - 7) * 2 * M_PI / 5), sin(M_PI / 5 + (11 - 7) * 2 * M_PI / 5),-0.5);
    Vector3D p12 = Vector3D::point(0,0,-(sqrt(5)/2));
    vector<Vector3D> points = {p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12};
    for (const auto& point : points) {
        icosahedron.points.push_back(point);
    }

    //faces
    vector<vector<int>> faceIndices = {
            {1, 2, 3}, {1, 3, 4}, {1, 4, 5}, {1, 5, 6}, {1, 6, 2},
            {2, 7, 3}, {3, 7, 8}, {3, 8, 4}, {4, 8, 9}, {4, 9, 5},
            {5, 9, 10}, {5, 10, 6}, {6, 10, 11}, {6, 11, 2}, {2, 11, 7},
            {12, 8, 7}, {12, 9, 8}, {12, 10, 9}, {12, 11, 10}, {12, 7, 11}
    };

    for (const auto& indices : faceIndices) {
        icosahedron.faces.emplace_back(indices);
    }

    return icosahedron;
}

Figure createDodecahedron(){
    Figure dodecahedron;
    Figure icosahedron = createIcosahedron();

    vector<Vector3D> points;
    for (const auto& icoface : icosahedron.faces) {
        double sumX = 0, sumY = 0, sumZ = 0;
        for (const int pointindex : icoface.point_indexes) {
            const Vector3D& icopoint = icosahedron.points[pointindex];
            sumX += icopoint.x;
            sumY += icopoint.y;
            sumZ += icopoint.z;
        }
        double avgX = sumX / 3;
        double avgY = sumY / 3;
        double avgZ = sumZ / 3;

        Vector3D point = Vector3D::point(avgX, avgY, avgZ);
        dodecahedron.points.push_back(point);
    }



    vector<vector<int>> faceIndices = {
            {1, 2, 3, 4, 5}, {1, 6, 7, 8, 2}, {2, 8, 9, 10, 3},
            {3, 10, 11, 12, 4}, {4, 12, 13, 14, 5}, {5, 14, 15, 6, 1},
            {20, 19, 18, 17, 16}, {20, 15, 14, 13, 19}, {19, 13, 12, 11, 18},
            {18, 11, 10, 9, 17}, {17, 9, 8, 7, 16}, {16, 7, 6, 15, 20}
    };

    for (const auto& indices : faceIndices) {
        dodecahedron.faces.emplace_back(indices);
    }


    return dodecahedron;
}

void TriangleDivision(Figure& icosahedron){
    for (Face& face: icosahedron.faces) {
        //points ABC
        Vector3D fp1, fp2, fp3;
        int p1, p2, p3;
        for(int i = 0; i < face.point_indexes.size(); i++) {
            if(i==0) {
                p1 = face.point_indexes[i];
                fp1 = icosahedron.points[i];
            }
            else if(i==1) {
                p2 = face.point_indexes[i];
                fp2 = icosahedron.points[i];
            }
            else if(i==2){
                p3 = face.point_indexes[i];
                fp3 = icosahedron.points[i];
            }
        }
        //points DEF
        Vector3D mp1 = calculateMidpoint(fp1, fp2);
        Vector3D mp2 = calculateMidpoint(fp2, fp3);
        Vector3D mp3 = calculateMidpoint(fp3, fp1);
        int mpi1 = icosahedron.points.size() + 1;//D
        int mpi2 = icosahedron.points.size() + 2;//E
        int mpi3 = icosahedron.points.size() + 3;//F

        vector<Vector3D> midpoints = {mp1, mp2, mp3};
        for (const auto& midpoint : midpoints) {
            icosahedron.points.push_back(midpoint);
        }

        //faces
        vector<int> smallface1 = {p1, mpi1, mpi2};
        Face sf1(smallface1);
        vector<int> smallface2 = {p2, mpi3, mpi1};
        Face sf2(smallface2);
        vector<int> smallface3 = {p3, mpi2, mpi3};
        Face sf3(smallface3);
        vector<int> smallface4 = {mpi1, mpi3, mpi2};
        Face sf4(smallface4);

        vector<Face> faces = {sf1, sf2, sf3, sf4};
        for(const auto& face: faces){
            icosahedron.faces.push_back(face);
        }
    }
}

Figure createSphere(const double r, const int n){
    Figure icosahedron = createIcosahedron();
    for(int i = 0; i < n; i++){
        TriangleDivision(icosahedron);
    }
    for(auto& point: icosahedron.points){
        double r = sqrt((point.x*point.x)+(point.y*point.y)+(point.z*point.z));
        point.x /= r;
        point.y /= r;
        point.z /= r;
    }
    return icosahedron;
}

Vector3D calculateMidpoint(const Vector3D& p1, const Vector3D& p2) {
    double X = (p1.x + p2.x) / 2.0;
    double Y = (p1.y + p2.y) / 2.0;
    double Z = (p1.z + p2.z) / 2.0;

    return Vector3D::point(X, Y, Z);
}

Figure createCone(const int n, const double h){
    Figure cone;
    vector<Vector3D> points;
    for(int i = 0; i < n; i++) {
        Vector3D point = Vector3D::point(cos((2 * i * M_PI) / n), sin((2 * i * M_PI) / n), 0);
        points.push_back(point);


    }
    Vector3D pn = Vector3D::point(0,0,h);
    points.emplace_back(pn);

    for(int i = 0; i < n; i++){
        vector<int> sideface = {i, (i+1) % n, n};
        Face sf(sideface);
        cone.faces.emplace_back(sf);
    }
    // Ground face: {p(n-1), p(n-2), ..., p(0)}
    vector<int> fn;

    for(int m = n - 1; m >= 0; m--){
        fn.emplace_back(m);
    }
    Face GF(fn);
    cone.faces.emplace_back(GF);

    return cone;
}

Figure::Figure() : color(1.0, 1.0, 1.0) {

}
