#ifndef COMPONENT_H
#define COMPONENT_H

#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

using namespace cv;

class Component
{
public:
    Component(std::vector<Point2d> data);
    ~Component();
    std::vector<Point2d> points;
    float mean, variance, median;
    int minx, miny, maxx, maxy;
    bool stats;
    void getStats(Mat *SWTimage);
    int getWidth();
    int getHeigth();
    int isValid();
    Point2d getCenter();
    std::vector<Point2d> rotate_component(float angle);
    int isSimilar(Component similar);
    float angle2Components(Component *compon);
};



#endif // COMPONENT_H
