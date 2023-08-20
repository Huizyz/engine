#include <limits>
#include "Line2D.h"
#include "easy_image.h"
#include "cmath"

img::EasyImage draw2DLines(Lines2D& lines, const int size, Color &color) {

    //background colour
    img::Color backGroundColor;
    backGroundColor.blue = lround(color.blue * 255);
    backGroundColor.green = lround(color.green * 255);
    backGroundColor.red = lround(color.red * 255);

    // Determine the range of x and y coordinates
    double xmin = numeric_limits<double>::infinity();
    double xmax = xmin * -1;
    double ymin = numeric_limits<double>::infinity();
    double ymax = ymin * -1;

    // iterate over the list of lines to update the variables
    for (Line2D& line : lines) {
        if (line.p1.x < xmin) {
            xmin = line.p1.x;
        }
        if (line.p1.x > xmax) {
            xmax = line.p1.x;
        }
        if (line.p2.x < xmin) {
            xmin = line.p2.x;
        }
        if (line.p2.x > xmax) {
            xmax = line.p2.x;
        }
        if (line.p1.y < ymin) {
            ymin = line.p1.x;
        }
        if (line.p1.y > ymax) {
            ymax = line.p1.y;
        }
        if (line.p2.y < ymin) {
            ymin = line.p2.y;
        }
        if (line.p2.y > ymax) {
            ymax = line.p2.y;
        }
    }

    // calculate the size of the image
    double x_range = xmax - xmin;
    double y_range = ymax - ymin;
    double max_range = std::max(x_range, y_range);
    double image_x = size*(x_range / max_range);
    double image_y = size*(y_range / max_range);

    img::EasyImage image(lround(image_x) + 1,lround(image_y) + 1, backGroundColor);

    // calculate the scale factor d
    double d = 0.95 *(image_x/x_range);

    // multiply the coordinates of all points by d
    Lines2D copy_line;
    for (Line2D& line : lines) {
        line.p1.x *= d;
        line.p1.y *= d;
        line.p2.x *= d;
        line.p2.y *= d;
        copy_line.emplace_back(line);
    }

    // Move the line drawing
    double DC_x = d*((xmin+xmax)/2);
    double DC_y = d*((ymin+ymax)/2);
    double dx = (image_x/2.0) - DC_x;
    double dy = (image_y/2.0) - DC_y;

    // Sum of Linepoints with (dx,dy) + rounding of points
    for (Line2D& line : copy_line) {
        image.draw_line(lround(line.p1.x + dx),lround(line.p1.y + dy),lround(line.p2.x + dx),lround(line.p2.y + dy),img::Color(::lround(line.color.red * 255), ::lround(line.color.green * 255),::lround(line.color.blue * 255)));
    }
    return image;
}

Line2D::Line2D(const Point2D &p1, const Point2D &p2, const Color &color) : p1(p1), p2(p2), color(color) {}
