#ifndef ENGINE_ZBUFFER_H
#define ENGINE_ZBUFFER_H


#include <vector>
#include "easy_image.h"
#include "Color.h"

class ZBuffer: public std::vector<std::vector<double> >
{
public:
    //Constructor: maakt een Z-Buffer van de correcte
    //grootte aan en initialiseert alle velden op +inf
    ZBuffer(const int width, const int height);

};

#endif //ENGINE_ZBUFFER_H
