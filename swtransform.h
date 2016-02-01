#ifndef SWTRANSFORM_H
#define SWTRANSFORM_H
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
using namespace cv;

int strokeMedianFilter(Mat *SWTimage, std::vector< std::vector < Point2d> > *saved_rays);
int strokeWidth(int i, int j, short g_x, short g_y, Mat *SWTimage, std::vector< std::vector< Point2d> > *saved_rays, Mat *edge_image,
    Mat* gradient_x, Mat* gradient_y);
int getSteps(int g_x, int g_y, float& step_x, float& step_y);
Mat SWTransform(Mat* edge_image, Mat* gradient_x, Mat* gradient_y);

#endif // SWTRANSFORM_H
