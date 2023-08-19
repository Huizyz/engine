#define _USE_MATH_DEFINES
#include "easy_image.h"
#include "ini_configuration.h"
#include "l_parser.h"
#include "Line2D.h"
#include "vector3d.h"
#include "Color.h"
#include "Figure.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>



Lines2D drawLSystem(std::string filenaam, const Color &color) {
    //aan te passen
    LParser::LSystem2D l_system;
    std::ifstream in(filenaam);
    in >> l_system;
    in.close();

    set<char> alphabet = l_system.get_alphabet();
    double angle = l_system.get_angle();
    double current_angle = l_system.get_starting_angle();
    string initiator = l_system.get_initiator();
    unsigned int nr_iterations = l_system.get_nr_iterations();

    string draw;
    Lines2D  lines;

    for (int i=0; i < nr_iterations; i++) {
        for(char c: initiator){
            if(alphabet.find(c) != alphabet.end()){
                draw += l_system.get_replacement(c);
            }
            else if(c == '+' || c == '-' || c == '(' || c == ')'){
                draw += c;
            }
            else {
                cerr << "wrong character" << endl;
            }
        }
        initiator = draw;
        draw = "";
    }

    double x1 = 0.0;
    double y1 = 0.0;
    double x2;
    double y2;

    vector<pair<pair<double,double>,double>> punten;
    for(char c: initiator) {
        if(alphabet.find(c) != alphabet.end()){
            x2 = x1 + cos(current_angle * (M_PI / 180));
            y2 = y1 + sin(current_angle * (M_PI / 180));
            if(l_system.draw(c)){
                lines.emplace_back(Point2D(x1,y1), Point2D(x2,y2), color);
            }
            x1 = x2;
            y1 = y2;
        }
        else if(c == '+'){
            current_angle += angle;
        }
        else if(c == '-'){
            current_angle -= angle;
        }
        else if(c == '('){
            punten.emplace_back(make_pair(x1,y1),current_angle);
        }
        if(c == ')'){
            x1 = punten.back().first.first;
            y1 = punten.back().first.second;
            current_angle = punten.back().second;
            punten.pop_back();
        }

    }
    return lines;
}

img::EasyImage generate_image(const ini::Configuration &configuration)
{
    // img::EasyImage image;
    string type = configuration["General"]["type"].as_string_or_die();
    int size = configuration["General"]["size"].as_int_or_die();
    std::vector<double> backGroundInfo = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
    Color Color(backGroundInfo[0], backGroundInfo[1], backGroundInfo[2]);

    if (type == "2DLSystem") {
        // Read L-system configuration
        std::string inputFile = configuration["2DLSystem"]["inputfile"].as_string_or_die();
        std::vector<double> colorInfo = configuration["2DLSystem"]["color"].as_double_tuple_or_die();
        class Color lineColor(colorInfo[0], colorInfo[1], colorInfo[2]);
        // Generate+draw lines from L-system
        Lines2D lines = drawLSystem(inputFile, lineColor);
        return draw2DLines(lines, size, Color);
    }

    if (type == "Wireframe") {
        // Read wireframe configuration
        int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
        std::vector<double> eyeInfo = configuration["General"]["eye"].as_double_tuple_or_die();
        Vector3D eyeLocation = Vector3D::point(eyeInfo[0], eyeInfo[1], eyeInfo[2]);

        Figures3D figures;
        Lines2D lines;

        for (int l = 0; l < nrFigures; l++) {
            Figure figure;

            string figureName = "Figure" + to_string(l);
            std::string figureType = configuration[figureName]["type"].as_string_or_die();

            // get all the inputs for the matrix
            double rotationX = configuration[figureName]["rotateX"].as_double_or_die();
            double rotationY = configuration[figureName]["rotateY"].as_double_or_die();
            double rotationZ = configuration[figureName]["rotateZ"].as_double_or_die();
            double scale = configuration[figureName]["scale"].as_double_or_die();

            // get the location of the center in the figures for the matrix transformations
            std::vector<double> centerInfo = configuration[figureName]["center"].as_double_tuple_or_die();
            Vector3D centerLocation = Vector3D::vector(centerInfo[0], centerInfo[1], centerInfo[2]);

            // transformations
            Matrix transformation = Transformation(scale, rotationX, rotationY, rotationZ,centerLocation);

            if (figureType == "LineDrawing") {
                // make all the figure points + faces (lines)
                int nrPoints = configuration[figureName]["nrPoints"].as_int_or_die();
                int nrLines = configuration[figureName]["nrLines"].as_int_or_die();
                for (int j = 0; j < nrPoints; j++) {
                    string pointName = "point" + to_string(j);
                    std::vector<double> pointsInfo = configuration[figureName][pointName].as_double_tuple_or_die();
                    Vector3D point = Vector3D::point(pointsInfo[0], pointsInfo[1], pointsInfo[2]);
                    figure.points.push_back(point);
                }
                //make all the figure lines
                for (int h = 0; h < nrLines; h++) {
                    string lineName = "line" + to_string(h);
                    std::vector<int> linesInfo = configuration[figureName][lineName].as_int_tuple_or_die();
                    Face pointIndex = Face(linesInfo);
                    figure.faces.push_back(pointIndex);
                }
                // transform figure + push in list of all figures
            }
            else if (figureType == "Cube") {
                figure = createCube();
                applyTransformation(figure,transformation);
            }
            else if (figureType == "Tetrahedron") {
                figure = createTetrahedron();
                applyTransformation(figure,transformation);
            }
            else if (figureType == "Icosahedron") {
                figure = createIcosahedron();
                applyTransformation(figure,transformation);
            }
            else if (figureType == "Octahedron") {
                figure = createOctahedron();
                applyTransformation(figure,transformation);
            }
            else if (figureType == "Dodecahedron") {
                figure = createDodecahedron();
                applyTransformation(figure,transformation);
            }
            else if (figureType == "Cone") {
                int n = configuration[figureName]["n"].as_int_or_die();
                double height = configuration[figureName]["height"].as_double_or_die();
                figure = createCone(n,height);
                applyTransformation(figure,transformation);
            }
            else if (figureType == "Sphere") {
                //bekijk dit na!!
                int n = configuration[figureName]["n"].as_int_or_die();
                if (n == 0){
                    figure = createIcosahedron();
                    applyTransformation(figure,transformation);

                }
                else{
                    figure = createSphere(1, n);
                    applyTransformation(figure, transformation);

                }
            }
            //get the color information
            std::vector<double> colorInfo = configuration[figureName]["color"].as_double_tuple_or_die();
            class Color linesColor(colorInfo[0], colorInfo[1], colorInfo[2]);
            figure.color = linesColor;
            figures.emplace_back(figure);
            // apply the eyepoint transformation on all the figures and project everything to 2D to draw.
            Matrix eyepointMatrix = eyePointTrans(eyeLocation);
            lines = doProjection(figures);
            return draw2DLines(lines, size, Color);
        }
    }

}

int main(int argc, char const* argv[])
{
        int retVal = 0;
        try
        {
                std::vector<std::string> args = std::vector<std::string>(argv+1, argv+argc);
                if (args.empty()) {
                        std::ifstream fileIn("filelist");
                        std::string filelistName;
                        while (std::getline(fileIn, filelistName)) {
                                args.push_back(filelistName);
                        }
                }
                for(std::string fileName : args)
                {
                        ini::Configuration conf;
                        try
                        {
                                std::ifstream fin(fileName);
                                if (fin.peek() == std::istream::traits_type::eof()) {
                                    std::cout << "Ini file appears empty. Does '" <<
                                    fileName << "' exist?" << std::endl;
                                    continue;
                                }
                                fin >> conf;
                                fin.close();
                        }
                        catch(ini::ParseException& ex)
                        {
                                std::cerr << "Error parsing file: " << fileName << ": " << ex.what() << std::endl;
                                retVal = 1;
                                continue;
                        }

                        img::EasyImage image = generate_image(conf);
                        if(image.get_height() > 0 && image.get_width() > 0)
                        {
                                std::string::size_type pos = fileName.rfind('.');
                                if(pos == std::string::npos)
                                {
                                        //filename does not contain a '.' --> append a '.bmp' suffix
                                        fileName += ".bmp";
                                }
                                else
                                {
                                        fileName = fileName.substr(0,pos) + ".bmp";
                                }
                                try
                                {
                                        std::ofstream f_out(fileName.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
                                        f_out << image;

                                }
                                catch(std::exception& ex)
                                {
                                        std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                                        retVal = 1;
                                }
                        }
                        else
                        {
                                std::cout << "Could not generate image for " << fileName << std::endl;
                        }
                }
        }
        catch(const std::bad_alloc &exception)
        {
    		//When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
    		//Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
    		//(Unless of course you are already consuming the maximum allowed amount of memory)
    		//If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
		//mark the test as failed while in reality it just needed a bit more memory
                std::cerr << "Error: insufficient memory" << std::endl;
                retVal = 100;
        }
        return retVal;
}
