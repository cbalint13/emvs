/*********************************************************************
 * Software License Agreement (BSD License)
 *
 * Copyright (c) 2013, 2016
 *
 * Authors in alphabetical order:
 *
 * Balint Cristian <cristian dot balint at gmail dot com>
 * Nghia Ho <nghiaho12 at yahoo dot com> http://nghiaho.com
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the Willow Garage nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *********************************************************************/

#ifndef EMVS_H
#define EMVS_H

#include <map>
#include <vector>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/xfeatures2d.hpp>


using namespace std;
using namespace cv;
using namespace cv::xfeatures2d;

struct PointRGB : public Point3f
{
    uchar r;
    uchar g;
    uchar b;
};

typedef pair<PointRGB, uint> PointRGB_hit_type;

namespace emvs
{

// visualise data distribution
template <typename T> Mat plotGraph(vector<T>& vals, int YRange[2], float best, float second_best)
{

    auto it = minmax_element(vals.begin(), vals.end());
    float scale = 1./ceil(*it.second - *it.first); 
    float bias = *it.first;

    int rows = YRange[1] - YRange[0] + 1;
    Mat image = Mat::zeros( rows+20, vals.size(), CV_8UC3 );
    image.setTo(0);
    for (int i = 0; i < (int)vals.size()-1; i++)
    {
        line(image, Point(i, rows - 1 - (vals[i] - bias)*scale*YRange[1]),
             Point(i+1, rows - 1 - (vals[i+1] - bias)*scale*YRange[1]), Scalar(255, 0, 0), 1);

        if (vals[i]==best)
          circle( image, Point(i, rows - 1 - (vals[i] - bias)*scale*YRange[1]), 2, Scalar(0, 0, 255));
        else if (vals[i]==second_best)
          circle( image, Point(i, rows - 1 - (vals[i] - bias)*scale*YRange[1]), 2, Scalar(0, 255, 0));
        else
          circle( image, Point(i, rows - 1 - (vals[i] - bias)*scale*YRange[1]), 1, Scalar(0, 255, 255));
    }
    resize( image, image, Size( image.cols*2, image.rows*2 ) );

    return image;
}


class Camera
{

  public:

    Camera()
    {
      m_daisy_init = false;
    }

    void SetIntrinsic( float focal, int width, int height )
    {
      m_K(0,0) = focal;
      m_K(1,1) = focal;
      m_K(0,2) = width/2;
      m_K(1,2) = height/2;
      m_K(2,2) = 1.0;

      m_inv_K = m_K.inv();
    }

    void SetRotation( const Matx33f &R )
    {
      m_R = R;
      // cache
      m_inv_R = R.t();
      m_C = -m_inv_R * m_t;
    }

    void SetTranslation( const Matx31f &t )
    {
      m_t = t;
      // cache
      m_C = -m_inv_R * m_t;
    }

    Matx31f Position() const
    {
      // http://phototour.cs.washington.edu/bundler/bundler-v0.4-manual.html
      // -R' * t
      return m_C;
    }

    Matx31f ViewDirection() const
    {
      // http://phototour.cs.washington.edu/bundler/bundler-v0.4-manual.html
      // transform (0, 0, -1) from camera to world
      Matx31f d;

      d(0,0) = -m_inv_R(0,2);
      d(1,0) = -m_inv_R(1,2);
      d(2,0) = -m_inv_R(2,2);

      return d;
    }

    Point3f Unproject( float xp, float yp, float depth ) const
    {
      CV_Assert( depth > 0 );

      Matx31f x, ray, ret;
      x(0) = xp; x(1) = yp; x(2) = 1;

      ray = m_inv_K * x;
      ray *= depth;
      // -z into the screen
      ray(2) = -ray(2);

      ret = m_inv_R*( ray - m_t );
      return Point3f( ret(0), ret(1), ret(2) );
    }

    Point2f Project( float x, float y, float z, float *ret_depth = NULL ) const
    {
      Matx31f X, new_pt;

      X(0) = x; X(1) = y; X(2) = z;

      new_pt = m_R*X + m_t;

      if ( ret_depth )
      {
        *ret_depth = -new_pt(2);
      }

      new_pt(0) /= -new_pt(2);
      new_pt(1) /= -new_pt(2);
      new_pt(2) = 1.0f;

      new_pt = m_K*new_pt;

      return Point2f( new_pt(0), new_pt(1) );
    }

    Matx33f GetK() const { return m_K; }
    Matx33f GetR() const { return m_R; }
    Matx31f GetT() const { return m_t; }
    Matx33f GetInvK() const { return m_inv_K; }

    // task specific
    Mat m_image;
    Mat descriptors;
    Ptr<DAISY> m_daisy;
    bool m_daisy_init;
    vector<int> m_index_to_point_cloud;

  private:

    Matx31f m_t, m_C;
    Matx33f m_R, m_inv_R;
    Matx33f m_K, m_inv_K;
};

// IO routines
void LoadBundler( const string &bundle_file, vector<Camera> &cameras );
void LoadVisData( const string &bundle_file, const vector<Camera> &cameras, multimap<uint,uint> &visdata );
void OutputPLY( const string &filename, const vector<PointRGB_hit_type> &point_cloud_hit, uint min_hit );

// Stereo Routines
void UniformSampling( const Camera& camera1, const Camera& camera2, const Point2i& Img2DPoint,
                      float Lmin, float Lmax, float r, vector<Point2f>& epi, vector<Point3f>& Xs,
                      vector<float>& depths );

void DoStereoPair( Camera& camera1, Camera& camera2, vector<PointRGB_hit_type>& point_cloud_hit );

// Miscelaneous
int TermProgress( double dfComplete , int nLastTick = -1 );

} // end namespace emvs

#endif
