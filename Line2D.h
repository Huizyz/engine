#ifndef ENGINE_LINE2D_H
#define ENGINE_LINE2D_H
#include "Point2D.h"
#include "Color.h"
#include "list"
#include "easy_image.h"
using namespace std;

class Line2D {
public:
    Point2D p1;
    Point2D p2;
    Color color;

    //deze velden moeten worden ingevuld tijdens het projecteren
    double z1;
    double z2;

    Line2D(const Point2D &p1, const Point2D &p2, const Color &color);
};


using Lines2D = std::list<Line2D>;
img::EasyImage draw2DLines(Lines2D& lines, int size, Color &color);
#endif //ENGINE_LINE2D_H
