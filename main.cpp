#include "component.h"
#include "disjointsets.h"
#include "swtransform.h"

#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <unordered_map>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <boost/graph/connected_components.hpp>
#include <boost/unordered_map.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#define PI 3.14159265
using namespace cv;

std::vector< Component > legallyConnectedComponents(Mat* SWTimage);
std::vector< Component > filterComponents(Mat * SWTimage, std::vector< Component > components);
std::vector< std::vector< Component > > connectedLetters(std::vector< Component > valid_components);
void saveWords(std::vector < std::vector< Component > > words, Mat components_image);
void paintComponents(Mat *image, std::vector< Component > components);
void paintComponent(Mat *image, Component component);
void renderWordBoxes(std::vector<std::vector< Component > > words, Mat src);

/** @function main */
int main(int argc, char** argv)
{
    Mat src, src_gray;
    /// Load an image
    src = imread("road_signs.png");
    //src = imread("sign.jpg");

    if (!src.data)
    {
        std::cout << "Couldn't read the image" << std::endl;
        return -1;
    }

    /// Convert the image to grayscale
    cvtColor(src, src_gray, CV_BGR2GRAY);
    Mat filtered;
    // Reduce noise with a 3x3 gaussian kernel
    blur(src_gray, filtered, Size(3, 3));
    /// Canny detector
    Mat detected_edges;
    Canny(filtered, detected_edges, 124, 204, 3);
    //gradients
    Mat gradient_x, gradient_y;
    Sobel(filtered, gradient_x, CV_16S, 1, 0, 3, 1, 0, BORDER_DEFAULT);
    Sobel(filtered, gradient_y, CV_16S, 0, 1, 3, 1, 0, BORDER_DEFAULT);

    ///Get Stroke Width Transform
    Mat SWTimage = SWTransform(&detected_edges, &gradient_x, &gradient_y);

    ///Segement components by similar stroke width (Blob labeling)
    std::vector < Component > components;
    components = legallyConnectedComponents(&SWTimage);

    ///Filter components
    std::vector < Component > valid_components;
    valid_components = filterComponents(&SWTimage, components);

    ///Create words joining components
    std::vector< std::vector< Component > > words;
    words = connectedLetters(valid_components);
    //paint boxes on detected text
    renderWordBoxes(words, src);
    //Create windo canvas to show img
    namedWindow("detectedtext", CV_WINDOW_AUTOSIZE);
    //Show image in the name of the window
    imshow("detectedtext", src);


//////////////// UNCOMMENT IF YOU WANT TO PRINT INTERMEDIATE STAGES ///////////////////////////
//    Mat segmented_image = cv::Mat(src.rows, src.cols, CV_8UC3, cv::Scalar(0, 0, 0));
//    paintComponents(&segmented_image, components);
//    Mat filtered_components = cv::Mat(src.rows, src.cols, CV_8UC3, cv::Scalar(0, 0, 0));
//    paintComponents(&filtered_components, valid_components);
//    Mat words_image = cv::Mat(src.rows, src.cols, CV_8UC3, cv::Scalar(0, 0, 0));
//    for (std::vector<std::vector<Component > >::iterator it = words.begin(); it != words.end(); it++){
//        paintComponents(&words_image, *it);
//    }
//    //////////////////////Prettier SWTimage/////////////////
//    SWTimage.convertTo(SWTimage, CV_8U, 1.0);
//    uchar *ptr4;
//    for (int i = 0; i < SWTimage.rows; i++){
//        ptr4 = SWTimage.ptr<uchar>(i);
//        for (int j = 0; j < SWTimage.cols; j++){
//            if (ptr4[j] == 0){
//                ptr4[j] = 255;
//            }
//            else
//                ptr4[j] = (int)(ptr4[j] * 3 );
//        }
//    }
//    namedWindow("edges", CV_WINDOW_AUTOSIZE);
//    namedWindow("SWT", CV_WINDOW_AUTOSIZE);
//    namedWindow("segemented", CV_WINDOW_AUTOSIZE);
//    namedWindow("components_filtered", CV_WINDOW_AUTOSIZE);
//    namedWindow("words", CV_WINDOW_AUTOSIZE);
//    imshow("edges", detected_edges);
//    imshow("SWT", SWTimage);
//    imshow("segemented", segmented_image);
//    imshow("components_filtered", filtered_components);
//    imshow("words", words_image);

    /// Wait until user exit program by pressing a key
    waitKey(0);

    return 0;
}


void renderWordBoxes(std::vector< std::vector< Component > > words, Mat src){
    ///////////////////Boxes over original image////////////////
    for (std::vector< std::vector< Component > >::iterator it = words.begin(); it != words.end(); it++){
        if (it->size() > 2){
            int maxx = 0; int maxy = 0; int minx = 10000000; int miny = 10000000;
            for (std::vector<Component>::iterator it1 = it->begin(); it1 != it->end(); it1++){
                maxx = max(maxx, it1->maxx);
                maxy = max(maxy, it1->maxy);
                minx = min(minx, it1->minx);
                miny = min(miny, it1->miny);
            }
            Point vertex1;
            vertex1.x = (int)(minx * 0.98);
            vertex1.y = (int)(miny * 0.98);
            Point vertex2;
            vertex2.x = (int)(maxx * 1.02);
            vertex2.y = (int)(maxy * 1.02);
            rectangle(src, vertex1, vertex2, Scalar(0, 0, 255), 1);
        }
    }
}

void paintComponents(Mat *image, std::vector<Component> components){
    for (std::vector<Component>::iterator it = components.begin(); it != components.end(); ++it){
        paintComponent(image, *it);
    }
}

void paintComponent(Mat *image, Component component){
    int color0 = rand() % 256;
    int color1 = rand() % 256;
    int color2 = rand() % 256;
    for (std::vector<Point2d>::iterator it = component.points.begin(); it != component.points.end(); ++it){
        image->at<cv::Vec3b>(it->y, it->x)[0] = color0;
        image->at<cv::Vec3b>(it->y, it->x)[1] = color1;
        image->at<cv::Vec3b>(it->y, it->x)[2] = color2;
    }
}

void saveWords(std::vector < std::vector< Component > > words, Mat src){
    int name = 0;
    for (std::vector < std::vector< Component > >::iterator it = words.begin(); it != words.end(); it++){
        int maxx = 0; int maxy = 0; int minx = 10000000; int miny = 10000000;
        for (std::vector<Component>::iterator it2 = it->begin(); it2 != it->end(); it2++){
            maxx = max(maxx, it2->maxx);
            maxy = max(maxy, it2->maxy);
            minx = min(minx, it2->minx);
            miny = min(miny, it2->miny);
        }
        Point vertex1;
        vertex1.x = minx;
        vertex1.y = miny;
        Point vertex2;
        vertex2.x = maxx;
        vertex2.y = maxy;
        /*Mat dst;
        src(Rect(vertex1, vertex2)).copyTo(dst);*/
        std::string naming = std::to_string(name);
        Mat dst = src(Rect(minx, miny, maxx - minx, maxy - miny)).clone();
        imwrite(naming + ".jpg", dst);
        name++;
    }
}

std::vector<Component> filterComponents(Mat * SWTimage, std::vector<Component> components){
    std::vector<Component> valid_components;
    for (std::vector<Component>::iterator it = components.begin(); it != components.end(); it++){
        it->getStats(SWTimage);
        if (it->isValid())
            valid_components.push_back(*it);
    }
    return valid_components;
}

std::vector< std::vector<Component> > connectedLetters(std::vector<Component> valid_components){

    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
    const int num_vertex = valid_components.size();
    Graph g(num_vertex); //a graph is made with as much nodes as valid components

    int count = 0;
    for (std::vector<Component>::iterator it1 = valid_components.begin(); it1 != valid_components.end(); it1++){
        int count1 = 0;
        for (std::vector<Component>::iterator it2 = valid_components.begin(); it2 != valid_components.end(); it2++){
            if (it1->isSimilar(*it2)){ //if 2 components are similar( distance, size etc) an edge join them
                boost::add_edge(count, count1, g);
            }
            count1++;
        }
        count++;
    }


    //similar components with same directions(lines) are concatenated forming words
    std::vector<int> node_chain_lenght (num_vertex, 0);
    boost::graph_traits<Graph>::vertex_iterator init, end;
    typedef boost::property_map<Graph, boost::vertex_index_t>::type IndexMap;
    IndexMap index = get(boost::vertex_index, g);
    DisjointSets lines(num_vertex);
    for (boost::tie(init, end) = vertices(g); init != end; init++){ //all vertex iterator
        boost::graph_traits<Graph>::adjacency_iterator neig, neig_end;
        //neighbour iterator for each vertex
        for (boost::tie(neig, neig_end) = adjacent_vertices(*init, g); neig != neig_end; neig++){
            std::vector<int> line;
            line.push_back(index[*init]);
            line.push_back(index[*neig]);
            float angle = valid_components[index[*init]].angle2Components(&valid_components[index[*neig]]);

            boost::graph_traits<Graph>::adjacency_iterator neig2, neig2_end;
            boost::tie(neig2, neig2_end) = adjacent_vertices(*neig, g);
            DisjointSets visited_nodes(num_vertex);
            visited_nodes.Union(visited_nodes.FindSet(index[*init]), visited_nodes.FindSet(index[*neig]));

            int old = index[*init];
            int now = index[*neig];
            bool matched = true;
            while (matched == true){ //moving in graph while angles match

                matched = false;
                boost::graph_traits<Graph>::adjacency_iterator match;
                for (boost::tie(neig2, neig2_end); neig2 != neig2_end; neig2++){ //
                    float angle2 = valid_components[index[*neig]].angle2Components(&valid_components[index[*neig2]]);
                    if (std::max(angle, angle2)/std::min(angle, angle2) < 1.05){
                        if (visited_nodes.FindSet(index[*neig2]) != visited_nodes.FindSet(index[*init])){
                            match = neig2;
                            matched = true;
                        }
                    }
                }
                if (matched == true){ //move to matched node
                    visited_nodes.Union(visited_nodes.FindSet(index[*init]), visited_nodes.FindSet(index[*match]));
                    old = now;
                    now = index[*match];
                    boost::tie(neig2, neig2_end) = adjacent_vertices(*match, g);
                    line.push_back(index[*match]);
                }
            }

            int len = line.size();
            if (len > 2){ //if more than 2 nodes in a row it is considered a word
                for (std::vector<int>::iterator it = line.begin(); it != line.end(); it++){
                    lines.Union(lines.FindSet(*it), lines.FindSet(line[0]));
                    node_chain_lenght[*it] = max(len, node_chain_lenght[*it]);
                }
            }
        }
    }

    typedef std::unordered_map<int, std::vector< Component > > componentmap;
    componentmap map;

    //Store each line in a map
    for (int i = 0; i < valid_components.size(); i++){
        map[lines.FindSet(i)].push_back(valid_components[i]);
    }

    //Convert the map of lines in a vector of lines
    std::vector< std::vector< Component > > connected_letters;
    for (componentmap::iterator it = map.begin(); it != map.end(); it++){
        if (it->second.size() > 2)
            connected_letters.push_back(it->second);
    }

    //Write on terminal the length of each word found
    for (std::vector< std::vector< Component > >::iterator it =
        connected_letters.begin(); it != connected_letters.end(); it++){
        std::cout << "Word: ";
        for (std::vector<Component>::iterator it2 = it->begin(); it2 != it->end(); it2++){
            std::cout << "-";
        }
        std::cout << std::endl;
    }

    return connected_letters;
}

std::vector<Component> legallyConnectedComponents(Mat* SWTimage){
    Mat labels = cv::Mat(SWTimage->rows, SWTimage->cols, CV_16U, cv::Scalar(0));
    DisjointSets equivalent_connections(1); //created with one element (background)
    short next_label = 1;

    float* ptr;
    for (int i = 0; i < SWTimage->rows; i++){
        ptr = SWTimage->ptr<float>(i);
        for (int j = 0; j < SWTimage->cols; j++){
            if (ptr[j] > 0){
                std::vector<int> neighbour_labels;
                // check pixel to the west, west-north, north, north-east
                if (j - 1 >= 0) { //if not out of image
                    float west = SWTimage->at<float>(i, j - 1);
                    if (west > 0){
                        //if (ptr[j] / west <= 1.5 || west / ptr[j] <= 1.5){
                        if (max(ptr[j], west) / min(ptr[j], west) <= 2.0){
                            int west_label = labels.at<short>(i, j - 1);
                            if (west_label != 0)
                                neighbour_labels.push_back(west_label);
                        }
                    }
                }
                if (i - 1 >= 0) { //if not out of image
                    if (j - 1 >= 0){ //if not out of image
                        float west_north = SWTimage->at<float>(i - 1, j - 1);
                        if (west_north > 0){
                            //if (ptr[j] / west_north <= 1.5 || west_north / ptr[j] <= 1.5){
                            if (max(ptr[j], west_north) / min(ptr[j], west_north) <= 2.0){
                                int west_north_label = labels.at<short>(i - 1, j - 1);
                                if (west_north_label != 0)
                                    neighbour_labels.push_back(west_north_label);
                            }
                        }
                    }
                    float north = SWTimage->at<float>(i - 1, j);
                    if (north > 0){
                        //if (ptr[j] / north <= 1.5 || north / ptr[j] <= 1.5){
                        if (max(ptr[j], north) / min(ptr[j], north) <= 2.0){
                            int north_label = labels.at<short>(i - 1, j);
                            if (north_label != 0)
                                neighbour_labels.push_back(north_label);
                        }
                    }
                    if (j + 1 < SWTimage->cols){ //if not out of image
                        float north_east = SWTimage->at<float>(i - 1, j + 1);
                        if (north_east > 0){
                            //if (ptr[j] / north_east <= 1.5 || north_east / ptr[j] <= 1.5){
                            if (max(ptr[j], north_east) / min(ptr[j], north_east) <= 2.0){
                                int north_east_label = labels.at<short>(i - 1, j + 1);
                                if (north_east_label != 0)
                                    neighbour_labels.push_back(north_east_label);
                            }
                        }
                    }
                }
                //if neighbours have no label. Set new label(component)
                if (neighbour_labels.empty()){
                    labels.at<short>(i, j) = next_label;
                    next_label++;
                    equivalent_connections.AddElements(1);
                }
                else{//find minium label in connected neighbours. Save the equivalent connections for the second pass
                    int minimum = neighbour_labels[0];
                    if (neighbour_labels.size() > 1){
                        for (std::vector<int>::iterator it = neighbour_labels.begin(); it != neighbour_labels.end(); it++){
                            minimum = std::min(minimum, *it);
                            equivalent_connections.Union(equivalent_connections.FindSet(neighbour_labels[0]), equivalent_connections.FindSet(*it));
                        }
                    }
                    labels.at<short>(i, j) = minimum;
                }

            }
        }
    }
    //second pass. Assign component label from the equivalent connections list
    short * ptr2;
    typedef std::unordered_map< short, std::vector< Point2d > > componentmap;
    componentmap map;

    for (int i = 0; i < labels.rows; i++){
        ptr2 = labels.ptr<short>(i);
        for (int j = 0; j < labels.cols; j++){
            if (ptr2[j] != 0){
                ptr2[j] = equivalent_connections.FindSet(ptr2[j]);
                Point2d point;
                point.x = j;
                point.y = i;
                map[ptr2[j]].push_back(point);
            }
        }
    }

    std::vector<Component> components;
    for (componentmap::iterator it = map.begin(); it != map.end(); it++){
        components.push_back(Component(it->second));
    }


    return components;
}

