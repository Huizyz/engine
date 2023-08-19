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
    double xmin = lines.front().p1.x;
    double xmax = lines.front().p1.x;
    double ymin = lines.front().p1.y;
    double ymax = lines.front().p1.y;

    // iterate over the list of lines to update the variables
    for (Line2D& line : lines) {
        xmin = std::min(xmin, std::min(line.p1.x, line.p2.x));
        xmax = std::max(xmax, std::max(line.p1.x, line.p2.x));
        ymin = std::min(ymin, std::min(line.p1.y, line.p2.y));
        ymax = std::max(ymax, std::max(line.p1.y, line.p2.y));
    }

    // calculate the size of the image
    double x_range = xmax - xmin;
    double y_range = ymax - ymin;
    double max_range = std::max(x_range, y_range);
    double image_x = size*(x_range / max_range);
    double image_y = size*(y_range / max_range);

    img::EasyImage image(lround(image_x),lround(image_y), backGroundColor);

    // calculate the scale factor d
    double d = 0.95 *(image_x/x_range);

    // multiply the coordinates of all points by d
    Lines2D copy_line;
    for (Line2D& line : lines) {
        line.p1.x *= d;
        line.p1.y *= d;
        line.p2.x *= d;
        line.p2.y *= d;
        copy_line.push_back(line);
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
