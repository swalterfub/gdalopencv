/*
 * gdal_image.cpp -- Load GIS data into OpenCV Containers using the Geospatial Data Abstraction Library
*/
// OpenCV Headers
#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"

// GDAL Dependencies
#include <cpl_conv.h>
#include <gdal_priv.h>

//  Boost Dependencies
#include <boost/filesystem.hpp>

// C++ Standard Libraries
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

// STL Dependencies
#include <string>

#include <sstream>

using namespace cv;
using namespace std;

// define the corner points
//    Note that GDAL library can natively determine this
cv::Point2d tl( -122.441017, 37.815664 );
cv::Point2d tr( -122.370919, 37.815311 );
cv::Point2d bl( -122.441533, 37.747167 );
cv::Point2d br( -122.3715,   37.746814 );
// determine dem corners
cv::Point2d dem_bl( -122.0, 38);
cv::Point2d dem_tr( -123.0, 37);
// range of the heat map colors
std::vector<std::pair<cv::Vec3b,double> > color_range;
// List of all function prototypes
cv::Point2d lerp( const cv::Point2d&, const cv::Point2d&, const double& );
cv::Vec3b get_dem_color( const double& );
cv::Point2d world2dem( const cv::Point2d&, const cv::Size&);
cv::Point2d pixel2world( const int&, const int&, const cv::Size& );
void add_color( cv::Vec3b& pix, const uchar& b, const uchar& g, const uchar& r );
/*
 * Linear Interpolation
 * p1 - Point 1
 * p2 - Point 2
 * t  - Ratio from Point 1 to Point 2
*/
cv::Point2d lerp( cv::Point2d const& p1, cv::Point2d const& p2, const double& t ){
    return cv::Point2d( ((1-t)*p1.x) + (t*p2.x),
                        ((1-t)*p1.y) + (t*p2.y));
}
/*
 * Interpolate Colors
*/
template <typename DATATYPE, int N>
cv::Vec<DATATYPE,N> lerp( cv::Vec<DATATYPE,N> const& minColor,
                          cv::Vec<DATATYPE,N> const& maxColor,
                          double const& t ){
    cv::Vec<DATATYPE,N> output;
    for( int i=0; i<N; i++ ){
        output[i] = (uchar)(((1-t)*minColor[i]) + (t * maxColor[i]));
    }
    return output;
}
/*
 * Compute the dem color
*/
cv::Vec3b get_dem_color( const double& elevation ){
    // if the elevation is below the minimum, return the minimum
    if( elevation < color_range[0].second ){
        return color_range[0].first;
    }
    // if the elevation is above the maximum, return the maximum
    if( elevation > color_range.back().second ){
        return color_range.back().first;
    }
    // otherwise, find the proper starting index
    int idx=0;
    double t = 0;
    for( int x=0; x<(int)(color_range.size()-1); x++ ){
        // if the current elevation is below the next item, then use the current
        // two colors as our range
        if( elevation < color_range[x+1].second ){
            idx=x;
            t = (color_range[x+1].second - elevation)/
                (color_range[x+1].second - color_range[x].second);
            break;
        }
    }
    // interpolate the color
    return lerp( color_range[idx].first, color_range[idx+1].first, t);
}
/*
 * Given a pixel coordinate and the size of the input image, compute the pixel location
 * on the DEM image.
*/
cv::Point2d world2dem( cv::Point2d const& coordinate, const cv::Size& dem_size   ){
    // relate this to the dem points
    // ASSUMING THAT DEM DATA IS ORTHORECTIFIED
    double demRatioX = ((dem_tr.x - coordinate.x)/(dem_tr.x - dem_bl.x));
    double demRatioY = 1-((dem_tr.y - coordinate.y)/(dem_tr.y - dem_bl.y));
    cv::Point2d output;
    output.x = demRatioX * dem_size.width;
    output.y = demRatioY * dem_size.height;
    return output;
}
/*
 * Convert a pixel coordinate to world coordinates
*/
cv::Point2d pixel2world( const int& x, const int& y, const cv::Size& size ){
    // compute the ratio of the pixel location to its dimension
    double rx = (double)x / size.width;
    double ry = (double)y / size.height;
    // compute LERP of each coordinate
    cv::Point2d rightSide = lerp(tr, br, ry);
    cv::Point2d leftSide  = lerp(tl, bl, ry);
    // compute the actual Lat/Lon coordinate of the interpolated coordinate
    return lerp( leftSide, rightSide, rx );
}
/*
 * Add color to a specific pixel color value
*/
void add_color( cv::Vec3b& pix, const uchar& b, const uchar& g, const uchar& r ){
    if( pix[0] + b < 255 && pix[0] + b >= 0 ){ pix[0] += b; }
    if( pix[1] + g < 255 && pix[1] + g >= 0 ){ pix[1] += g; }
    if( pix[2] + r < 255 && pix[2] + r >= 0 ){ pix[2] += r; }
}

/**
 * Compute the scale factor for the conversion between image depths
 *
 * @param[in] gdalDepth   - GDAL Depth Type
 * @param[in] opencvDepth - OpenCV Depth Type
 * @return Scale factor between depths
 */
double gdal2OpenCVScale( const int& gdalDepth, const int& opencvDepth ){

    if( opencvDepth == CV_8U  && gdalDepth == GDT_Byte   ) return 1;
    if( opencvDepth == CV_8U  && gdalDepth == GDT_Int16  ) return 1/16.0;
    if( opencvDepth == CV_8U  && gdalDepth == GDT_UInt16 ) return 1/16.0;
    if( opencvDepth == CV_8U  && gdalDepth == GDT_Float32 ) return 1./255.;
    if( opencvDepth == CV_16U && gdalDepth == GDT_Byte   ) return 16;
    if( opencvDepth == CV_16U && gdalDepth == GDT_Int16  ) return 16;
    if( opencvDepth == CV_16U && gdalDepth == GDT_UInt16 ) return 1;
    if( opencvDepth == CV_16S && gdalDepth == GDT_Byte   ) return 16;
    if( opencvDepth == CV_16S && gdalDepth == GDT_Int16  ) return 16;
    if( opencvDepth == CV_16S && gdalDepth == GDT_UInt16 ) return 1;

    if( opencvDepth == CV_32F && gdalDepth == GDT_Float32 ) return 1.;

    throw string("Error: Unknown OpenCV Type or Unknown GDAL Type");
}

/**
 *  Convert and set the value from the gdal raster to the opencv image
 *
 * @param[in] val Pixel value from GDAL dataset
 * @param[in] image OpenCV Image to load
 * @param[in] point Pixel to set the value to
 * @param[in] scale Scale factor to convert with
*/
void GDAL2OpenCV( double const& val, Mat& image, Point const& point, const double& scale ){

    // get the required depth
    //int cvDepth = image.depth();

    // get the value for the image
    double pixVal = val * scale;

    // place the value
    if( image.depth() == CV_8U ){
        image.at<uchar>(point) = (uchar)pixVal;
    }
    else if( image.depth() == CV_16U ){
        image.at<ushort>(point) = (ushort)pixVal;
    }
    else if( image.depth() == CV_16S ){
        image.at<short>(point) = (short)pixVal;
    }
    else if( image.depth() == CV_32S ){
        image.at<int>(point) = (int)pixVal;
    }
    else if( image.depth() == CV_32F ){
        image.at<float>(point) = (float)pixVal;
    }
    else if( image.depth() == CV_64F ){
        image.at<double>(point) = (double)pixVal;
    }
    else{
        throw string("Error: Unknown Depth");
    }

}

/**
 * Convert Type 2 Depth
 *
 * @param[in] type OpenCV Type to convert
 * @return Associated depth
*/
int cvType2Depth( const int& type ){

    switch( type ){

        case CV_8UC1:
            return CV_8U;
        case CV_8UC2:
            return CV_8U;
        case CV_8UC3:
            return CV_8U;
        case CV_16UC1:
            return CV_16U;
        case CV_16UC2:
            return CV_16U;
        case CV_16UC3:
            return CV_16U;
        case CV_16SC1:
            return CV_16S;
        case CV_16SC2:
            return CV_16S;
        case CV_16SC3:
            return CV_16S;
        case CV_32SC1:
            return CV_32S;
        case CV_32SC2:
            return CV_32S;
        case CV_32SC3:
            return CV_32S;
        case CV_32FC1:
            return CV_32F;
        case CV_32FC2:
            return CV_32F;
        case CV_32FC3:
            return CV_32F;
        case CV_64FC1:
            return CV_64F;
        case CV_64FC2:
            return CV_64F;
        case CV_64FC3:
            return CV_64F;
        default:
            return -1;
    }
}


/**
 * Extract the number of channels from the type
 *
 * @param[in] type
 * @return Number of channels
*/
int cvType2Channels( const int& type ){

    switch( type ){

        case CV_8UC1:
        case CV_16UC1:
        case CV_16SC1:
        case CV_32SC1:
        case CV_32FC1:
        case CV_64FC1:
            return 1;
        case CV_8UC2:
        case CV_16UC2:
        case CV_16SC2:
        case CV_32SC2:
        case CV_32FC2:
        case CV_64FC2:
            return 2;
        case CV_8UC3:
        case CV_16UC3:
        case CV_16SC3:
        case CV_32SC3:
        case CV_32FC3:
        case CV_64FC3:
            return 3;

        default:
            return 0;
    }
}

/**
 * Convert Depth and Channels into Type
 *
 * @param[in] Depth OpenCV Depth value
 * @param[in] Channels Number of channels
 * @return  OpenCV Type
*/
int cvDepthChannel2Type( const int Depth, const int Channels ){

    if( Depth == CV_8U  && Channels == 3 )  return CV_8UC3;
    if( Depth == CV_8U  && Channels == 2 )  return CV_8UC2;
    if( Depth == CV_8U  && Channels == 1 )  return CV_8UC1;
    if( Depth == CV_16U && Channels == 3 ) return CV_16UC3;
    if( Depth == CV_16U && Channels == 2 ) return CV_16UC2;
    if( Depth == CV_16U && Channels == 1 ) return CV_16UC1;
    if( Depth == CV_16S && Channels == 3 ) return CV_16SC3;
    if( Depth == CV_16S && Channels == 2 ) return CV_16SC2;
    if( Depth == CV_16S && Channels == 1 ) return CV_16SC1;
    if( Depth == CV_32S && Channels == 3 ) return CV_32SC3;
    if( Depth == CV_32S && Channels == 2 ) return CV_32SC2;
    if( Depth == CV_32S && Channels == 1 ) return CV_32SC1;
    if( Depth == CV_32F && Channels == 1 ) return CV_32FC1;
    if( Depth == CV_64F && Channels == 1 ) return CV_64FC1;

    throw string("Error: combo not supported");

    return 0;
}

/**
 * Checks if the OpenCV Image Type is supported for this.
 *
 * 2 channel images are not supported only because I am not sure what
 * types of images these are.
*/
bool validOpenCVImageType( const int& imageType ){

    switch( imageType ){

        case CV_8UC1:
        case CV_8UC3:
        case CV_16UC1:
        case CV_16UC2:
        case CV_16UC3:
        case CV_16SC1:
        case CV_16SC2:
        case CV_16SC3:
        case CV_32SC1:
        case CV_32SC2:
        case CV_32SC3:
        case CV_32FC1:
        case CV_32FC2:
        case CV_32FC3:
        case CV_64FC1:
        case CV_64FC2:
        case CV_64FC3:
            return true;
        default:
            return false;
    }
}

/**
 * Read an image format using GDAL as the toolset.
 *
 * @param[in] filename Image filename
 * @param[in] flags Type of image you want to read
 */
Mat  imread_geo( const string& filename, const int& imageType ){

    // parse the flags to determine which color type we want to read
    if( validOpenCVImageType(imageType) == false ) throw string("Invalid Image Type");

    // ensure the file exists
    if( boost::filesystem::exists( filename ) == false ){
        throw string( "Error: File ") + filename + string(" does not exist" );
    }

    // register the gdal driver
    GDALAllRegister();

    // load the dataset
    GDALDataset*  dataset = (GDALDataset*) GDALOpen( filename.c_str(), GA_ReadOnly);

    // if the dataset returns null, then the file does not exist, or there was a read error
    if( dataset == NULL ){ throw string("Error:  GDAL Dataset returned null from read"); }

    // check if pixel data even exists
    if( dataset->GetRasterCount() <= 0 ) throw string("Error: File does not contain pixel data");

    //get the driver infomation
    //GDALDriver*  driver = dataset->GetDriver();

    // get raster image size
    Size imgSize( dataset->GetRasterXSize(), dataset->GetRasterYSize() );

    // create mats for each raster layer
    vector<Mat> layers( dataset->GetRasterCount() );
    for( size_t i=0; i<layers.size(); i++ )
        layers[i] = Mat( imgSize, cvDepthChannel2Type(cvType2Depth( imageType), 1));

    // iterate through each band
    for (int i = 0; i < dataset->GetRasterCount(); i++) {

        //create the band object
        GDALRasterBand *band = dataset->GetRasterBand(i + 1);

        // get the gdal band type
        int datatype = band->GetRasterDataType();

        // compute the scale factor
        double scale = gdal2OpenCVScale( datatype, layers[i].depth());

        //read each image row
        for ( int r = 0; r < layers[i].rows; r++) {

            float* pafScanline;
            pafScanline = (float*) CPLMalloc(sizeof (float) *layers[i].cols);
            CPLErr eErr = band->RasterIO(GF_Read, 0, r, layers[i].cols, 1, pafScanline, layers[i].cols, 1, GDT_Float32, 0, 0);
            if ( eErr != CE_None ) {
                        CPLFree(pafScanline);
                        Mat m;
                        return m;
                    }
            // iterate through each column
            for ( int c = 0; c < layers[i].cols; c++) {
                GDAL2OpenCV( pafScanline[c], layers[i], Point(c,r), scale );
            }
            CPLFree(pafScanline);
        }

    }

    // close our GDAL Dataset and clean up
    GDALClose( dataset );

    //merge channels into single image
    Mat output;
    if( layers.size() > 1 )
        cv::merge( layers, output );
    else
        output = layers[0].clone();

    // do channel conversions if necessary
    if( cvType2Channels(imageType) == output.channels() )
        return output;

    if( cvType2Channels(imageType) == 1 && output.channels() == 3 ){
        cvtColor( output, output, CV_BGR2GRAY );
        return output;
    }
    if( cvType2Channels(imageType) == 3 && output.channels() == 1 ){
        cvtColor( output, output, CV_GRAY2BGR);
        return output;
    }

    // Return the output image
    return output;
}

/*
 * Main Function
*/
int main( int argc, char* argv[] ){
    /*
     * Check input arguments
    */
    /*if( argc < 3 ){
        cout << "usage: " << argv[0] << " <image_name> <dem_model_name>" << endl;
        return -1;
    }
    */
    // load the image (note that we don't have the projection information.  You will
    // need to load that yourself or use the full GDAL driver.  The values are pre-defined
    // at the top of this file
    //cv::Mat image = cv::imread(argv[1], cv::IMREAD_LOAD_GDAL| cv::IMREAD_ANYDEPTH );
    //if (image.empty()){ throw std::runtime_error("ERROR loading image"); }

    //Mat inc_x = imread_geo( argv[1]+string("_0"), CV_32FC1 );
    //Mat inc_y = imread_geo( argv[1]+string("_1"), CV_32FC1 );
    //Mat inc_z = imread_geo( argv[1]+string("_2"), CV_32FC1 );
    //cout << "Using " << inc_z.size() << " channels" << endl;

    std::vector<cv::Mat> inc;
    inc.push_back(imread_geo( argv[1]+string("_0"), CV_32FC1 ));
    inc.push_back(imread_geo( argv[1]+string("_1"), CV_32FC1 ));
    inc.push_back(imread_geo( argv[1]+string("_2"), CV_32FC1 ));

    cout << "Using " << inc.size() << " Size" << endl;
    Mat3f mat_inc;
    cv::merge(inc,mat_inc);

    cout << "Using " << mat_inc.rows << " rows with " << mat_inc.cols << " cols and " << mat_inc.channels() << " channels with " << mat_inc.size() << " size each." << endl;

    //std::vector<Mat> emi;
    //Mat emi_x = imread_geo( argv[1]+string("_3"), CV_32FC1 );
    //Mat emi_y = imread_geo( argv[1]+string("_4"), CV_32FC1 );
    //Mat emi_z = imread_geo( argv[1]+string("_5"), CV_32FC1 );

    std::vector<Mat> norm;
    norm.push_back(imread_geo( argv[1]+string("_6"), CV_32FC1 ));
    norm.push_back(imread_geo( argv[1]+string("_7"), CV_32FC1 ));
    norm.push_back(imread_geo( argv[1]+string("_8"), CV_32FC1 ));
    Mat3f mat_norm;
    cv::merge(norm,mat_norm);

    //cv::Mat dot_array=[inc_norm, mat_norm];

    //Mat dot = Mat::dot(dot_array);
    // load the dem model
    /*cv::Mat dem = cv::imread(argv[2], cv::IMREAD_LOAD_GDAL | cv::IMREAD_ANYDEPTH );
    // create our output products
    cv::Mat output_dem(   image.size(), CV_8UC3 );
    cv::Mat output_dem_flood(   image.size(), CV_8UC3 );
    // for sanity sake, make sure GDAL Loads it as a signed short
    //if( dem.type() != CV_16SC1 ){ throw std::runtime_error("DEM image type must be CV_16SC1"); }
    // define the color range to create our output DEM heat map
    //  Pair format ( Color, elevation );  Push from low to high
    //  Note:  This would be perfect for a configuration file, but is here for a working demo.
    color_range.push_back( std::pair<cv::Vec3b,double>(cv::Vec3b( 188, 154,  46),   -1));
    color_range.push_back( std::pair<cv::Vec3b,double>(cv::Vec3b( 110, 220, 110), 0.25));
    color_range.push_back( std::pair<cv::Vec3b,double>(cv::Vec3b( 150, 250, 230),   20));
    color_range.push_back( std::pair<cv::Vec3b,double>(cv::Vec3b( 160, 220, 200),   75));
    color_range.push_back( std::pair<cv::Vec3b,double>(cv::Vec3b( 220, 190, 170),  100));
    color_range.push_back( std::pair<cv::Vec3b,double>(cv::Vec3b( 250, 180, 140),  200));
    // define a minimum elevation
    double minElevation = -13575;
    // iterate over each pixel in the image, computing the dem point
    for( int y=0; y<image.rows; y++ ){
    for( int x=0; x<image.cols; x++ ){
        // convert the pixel coordinate to lat/lon coordinates
        cv::Point2d coordinate = pixel2world( x, y, image.size() );
        // compute the dem image pixel coordinate from lat/lon
        cv::Point2d dem_coordinate = world2dem( coordinate, dem.size() );
        // extract the elevation
        double dz;
        if( dem_coordinate.x >=    0    && dem_coordinate.y >=    0     &&
            dem_coordinate.x < dem.cols && dem_coordinate.y < dem.rows ){
            dz = dem.at<short>(dem_coordinate);
        }else{
            dz = minElevation;
        }
        // write the pixel value to the file
        output_dem_flood.at<cv::Vec3b>(y,x) = image.at<cv::Vec3b>(y,x);
        // compute the color for the heat map output
        cv::Vec3b actualColor = get_dem_color(dz);
        output_dem.at<cv::Vec3b>(y,x) = actualColor;
        // show effect of a 10 meter increase in ocean levels
        if( dz < 10 ){
            add_color( output_dem_flood.at<cv::Vec3b>(y,x), 90, 0, 0 );
        }
        // show effect of a 50 meter increase in ocean levels
        else if( dz < 50 ){
            add_color( output_dem_flood.at<cv::Vec3b>(y,x), 0, 90, 0 );
        }
        // show effect of a 100 meter increase in ocean levels
        else if( dz < 100 ){
            add_color( output_dem_flood.at<cv::Vec3b>(y,x), 0, 0, 90 );
        }
    }}
    */
    // print our heat map
    cv::imwrite( "heat-map.tif"   ,  mat_inc );
    // print the flooding effect image
    //cv::imwrite( "flooded.jpg",  output_dem_flood);
    return 0;
}
