#include "component.h"

#define PI 3.14159265

Component::Component(std::vector<Point2d> data)
{
    points = data;
    stats = false;
}


Component::~Component()
{
}

void Component::getStats(Mat * SWTimage){
    std::vector<float> temp;
    temp.reserve(this->points.size());
    this->variance = 0;
    this->mean = 0;
    this->minx = 1000000;
    this->miny = 1000000;
    this->maxx = 0;
    this->maxy = 0;
    for (std::vector<Point2d>::iterator it = this->points.begin(); it != this->points.end(); it++){
        float stroke = SWTimage->at<float>(it->y, it->x);
        this->mean += stroke;
        temp.push_back(stroke);
        this->minx = min(minx, (int)it->x);
        this->miny = min(miny, (int)it->y);
        this->maxx = max(maxx, (int)it->x);
        this->maxy = max(maxy, (int)it->y);
    }
    mean = mean / ((float)points.size());
    for (std::vector<float>::const_iterator it = temp.begin(); it != temp.end(); it++) {
        variance += (*it - mean) * (*it - mean);
    }
    variance = variance / ((float)points.size());
    std::sort(temp.begin(), temp.end());
    median = temp[temp.size() / 2];
    stats = true;
}

int Component::getWidth(){
    assert(stats);
    return (maxx - minx);
}

int Component::getHeigth(){
    assert(stats);
    return(maxy - miny);
}

Point2d Component::getCenter(){
    assert(stats);
    Point2d center;
    center.x = (int)(maxx - getWidth()/2. + 1);
    center.y = (int)(maxy - getHeigth() / 2. + 1);
    return center;
}

std::vector<Point2d> Component::rotate_component(float rads){


    Point2d center = getCenter();
    std::vector<float> angles;
    std::vector<float> hipotenusas;
    for (std::vector<Point2d>::iterator it = points.begin(); it != points.end(); it++){
        float h = sqrt((it->x - center.x)*(it->x - center.x) + (it->y - center.y)*(it->y - center.y));
        hipotenusas.push_back(h);
        float cateto_opuesto = it->y - center.y;
        float cateto_contiguo = it->x - center.x;
        //float angle = asin((h / cateto_opuesto) + 2*PI);
        float angle = atan2(cateto_opuesto, -cateto_contiguo) + PI; //angle 0 to 360

        angles.push_back(angle);
    }

    std::vector<Point2d> rotated_component;
    float tminx = this->minx;
    float tmaxx = this->maxx;
    float tminy = this->miny;
    float tmaxy = this->maxy;
    for (int i = 0; i < points.size(); i++){
        Point2d rotated_point;
        float rojo = angles[i];
        float verde = hipotenusas[i];
        float new_y = -(hipotenusas[i] * sin(rads + angles[i])) + center.y;
        float new_x = hipotenusas[i] * cos(rads + angles[i]) + center.x;
        rotated_point.x = (int)(new_x);
        rotated_point.y = (int)(new_y);
        rotated_component.push_back(rotated_point);
    }

    return rotated_component;
}

int Component::isValid(){
    assert(stats);
    /*int width = (int)(maxx - minx + 1);
    int height = (int)(maxy - miny + 1);*/
    if (variance > 0.6 * mean){
        return 0;
    }
    float width = (float)(maxx - minx);
    float height = (float)(maxy - miny);

    // check font height
    if (height > 200) {
        return 0;
    }

    if (width > 200) {
        return 0;
    }

    if (height < 8) {
        return 0;
    }
    if (points.size() < 10) {
        return 0;
    }

    int min_width_rotated = 10000000, min_height_rotated = 10000000;
    int height_ofminwidth, width_ofminheight;
    float increment = 1. / 36.;
    for (float theta = increment * PI; theta < PI / 2.0; theta += increment * PI){
        std::vector<Point2d> rotated_points = rotate_component(theta);
        int minx_rotated = 1000000;
        int miny_rotated = 1000000;
        int maxx_rotated = 0;
        int maxy_rotated = 0;
        for (std::vector<Point2d>::iterator it = rotated_points.begin(); it != rotated_points.end(); it++){
            minx_rotated = min(minx_rotated, (int)it->x);
            miny_rotated = min(miny_rotated, (int)it->y);
            maxx_rotated = max(maxx_rotated, (int)it->x);
            maxy_rotated = max(maxy_rotated, (int)it->y);
        }
        int width_rotated = maxx_rotated - minx_rotated;
        int height_rotated = maxy_rotated - miny_rotated;
        if (min_width_rotated != min(min_width_rotated, width_rotated)){
            min_width_rotated = min(min_width_rotated, width_rotated);
            height_ofminwidth = maxy_rotated - miny_rotated;
        }
        if (min_height_rotated != min(min_height_rotated, height_rotated)){
            min_height_rotated = min(min_height_rotated, height_rotated);
            width_ofminheight = maxx - minx;
        }

    }
    if (min_width_rotated == 0 || min_height_rotated == 0 || height_ofminwidth == 0 || width_ofminheight == 0){
        return 0;
    }
    if (min_width_rotated < 2 || min_height_rotated < 2){
        return 0;
    }
    float ratio = ((float)min_width_rotated / (float)height_ofminwidth);
    float ratio2 = min_height_rotated / width_ofminheight;
    if ((float)min_width_rotated / (float)height_ofminwidth < 1. / 10. || (float)min_height_rotated / (float)width_ofminheight < 1. /10.) {
        return 0;
    }

    /*float area = length * width;
    float rminx = (float)minx;
    float rmaxx = (float)maxx;
    float rminy = (float)miny;
    float rmaxy = (float)maxy;
     compute the rotated bounding box
    float increment = 1. / 36.;
    for (float theta = increment * PI; theta<PI / 2.0; theta += increment * PI) {
        float xmin, xmax, ymin, ymax, xtemp, ytemp, ltemp, wtemp;
        xmin = 1000000;
        ymin = 1000000;
        xmax = 0;
        ymax = 0;
        for (unsigned int i = 0; i < points.size(); i++) {
            xtemp = points[i].x * cos(theta) + points[i].y * -sin(theta);
            ytemp = points[i].x * sin(theta) + points[i].y * cos(theta);
            xmin = std::min(xtemp, xmin);
            xmax = std::max(xtemp, xmax);
            ymin = std::min(ytemp, ymin);
            ymax = std::max(ytemp, ymax);
        }
        ltemp = xmax - xmin + 1;
        wtemp = ymax - ymin + 1;
        if (ltemp*wtemp < area) {
            area = ltemp*wtemp;
            length = ltemp;
            width = wtemp;
        }
    }
     //check if the aspect ratio is between 1/10 and 10
    if (length / width < 1. / 10. || length / width > 10.) {
        return 0;
    }*/

    return 1;
}

int Component::isSimilar(Component similar){
    assert(similar.stats);
    assert(stats);
    //similar median
    //if ((float)(similar.median) > (float)(median*1.6) || (float)(similar.median) < (float)(median*0.5))
    if (float(max(similar.median, this->median) / min(similar.median, this->median)) > 2.0)
        return 0;

    //similar height
    //if ((float)(similar.getHeigth()) > (float)(this->getHeigth()*1.6) || (float)(similar.getHeigth()) < (float)(this->getHeigth()*0.5))
    if (float(max(similar.getHeigth(), this->getHeigth()) / min(similar.getHeigth(), this->getHeigth())) > 2.0)
        return 0;

    //similar widht
    //if ((float)(similar.getWidth()) > (float)(this->getWidth()*2.0) || (float)(similar.getWidth()) < (float)(this->getWidth()*0.5))
    if (float(max(similar.getWidth(), this->getWidth()) / min(similar.getWidth(), this->getWidth())) > 2.0)
        return 0;

    //Center distances are higher than 3 times candidate width
    int distance_horiz = abs(similar.getCenter().x - this->getCenter().x);
    int max_width = max(this->getWidth(), similar.getWidth());
    if ((float)distance_horiz > (float)(3 * max_width))
        return 0;

    int distance_vert = abs(similar.getCenter().y - this->getCenter().y);
    int min_heigth = max(this->getHeigth(), similar.getHeigth());
    if ((float)distance_vert > (float)(1 * min_heigth))
        return 0;

    return 1;
}

float Component::angle2Components(Component *compon) {
    Point2d origin = this->getCenter();
    Point2d destiny = compon->getCenter();
    float cateto_opuesto = destiny.y - origin.y;
    float cateto_contiguo = destiny.x - origin.x;
    float angle = atan2(cateto_opuesto, -cateto_contiguo) + PI; //angle 0 to 360
    return angle;
}
