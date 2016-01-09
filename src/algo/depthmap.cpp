/*********************************************************************
 * Software License Agreement (BSD License)
 *
 * Copyright (c) 2013, 2015
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

/*

 "Efficient Large Scale Multi-View Stereo for Ultra High Resolution
 Image Sets" E. Tola and C. Strecha and P. Fua, Machine Vision
 and Applications

 */

#include <iostream>


#include "emvs.hpp"

//#define DEBUG


using namespace std;
using namespace emvs;

static inline float FastSumSquaredDifference( const float* vec1, const float* vec2, int N )
{
    // handle multiple word
    CV_Assert(N % 4 == 0);

    const float *ptr1 = vec1;
    const float *ptr2 = vec2;
    float sum = 0.0f;

    for ( int i = 0; i < N; i += 4 )
    {
       __m128 a, b, c, d;
      float res;

      a = _mm_loadu_ps(ptr1);
      b = _mm_loadu_ps(ptr2);
      c = _mm_sub_ps(a, b);

      d = _mm_dp_ps(c, c, 0xff);

      _mm_store_ss(&res, d);

      ptr1 += 4;
      ptr2 += 4;
      sum += res;
    }

    return sum;
}

namespace emvs {

void DoStereoPair( Camera &camera1, Camera &camera2, vector<PointRGB_hit_type>& point_cloud_hit )
{

    const int border = 20;

    const float min_depth = 0.5;
    const float max_depth = 25;

    const float match_ratio_threshold = 0.8;

    const int non_maxima_ratio = 10.0f;
    const float max_dist_sq = 0.025*0.025;


    // -Z reversal identity matrix
    Matx33f Z = Matx33f::eye(); Z(2,2) = -1;

    // reference camera parameters
    Matx31f image_point;
    Matx33f R0t_Z_invK0 = camera1.GetR().t()*Z*camera1.GetInvK();

    CV_Assert( camera1.descriptors.cols == camera2.descriptors.cols );

    // maximum rows*cols amount per image pairs
    if ( camera1.m_index_to_point_cloud.empty() )
    {
        camera1.m_index_to_point_cloud.resize( camera1.m_image.cols * camera1.m_image.rows, -1 );
    }

    if ( camera2.m_index_to_point_cloud.empty() )
    {
        camera2.m_index_to_point_cloud.resize( camera2.m_image.cols * camera2.m_image.rows, -1 );
    }


    int nLastTick = -1;
    #pragma omp parallel for schedule(dynamic,1)
    for ( int y = border; y < camera1.m_image.rows - border; y++ )
    {

      float *desc1 = NULL;
      float *desc2 = NULL;

      // progress bar
      #pragma omp critical
      {
        nLastTick = TermProgress( (double)(y+border) / (double)(camera1.m_image.rows), nLastTick );
      }

      for ( int x = border; x < camera1.m_image.cols - border; x++ )
      {
        vector<Point2f> epi;
        vector<Point3f> Xs;
        vector<float> depths;
        vector<float> dists, supressed_dists;
        Point2i ref;

        ref.x = x;
        ref.y = y;

        /*
         * (1) First SubSampling in integer mode given ref (img1 x,y).
         *  min_rdepth = 0.5, max_rdepth = 25 (sufficient for any SfM ?!)
         *  step in img2 of 8 pixels for epiline scan.
         *  8 pixel should beless then  half of DAISY size so response gaps are fine.
         */
        UniformSampling( camera1, camera2, ref, min_depth, max_depth, 8, epi, Xs, depths );

        // numerical bounds: protect epipolar samples
        if ( epi.empty() ) continue;
        // numerical bounds: too few samples ?
        if ( epi.size() < 18 ) continue;

        desc1 = camera1.descriptors.ptr<float>( y*camera2.m_image.cols + x );

        int best_i = -1;
        float best = FLT_MAX;
        float second_best = FLT_MAX;
        dists.resize(epi.size(), FLT_MAX);

#ifdef DEBUG
        vector<float> debug;
#endif
        int count = 0;
        for ( size_t i = 0; i < epi.size(); i++ )
        {

          int xv = epi[i].x;
          int yv = epi[i].y;

          if ( (xv < border) ||
               (xv >= camera2.m_image.cols - border) ||
               (yv < border) ||
               (yv >= camera2.m_image.rows - border)
             ) continue;


          desc2 = camera2.descriptors.ptr<float>( yv*camera2.m_image.cols + xv) ;
          float dist = FastSumSquaredDifference( desc1, desc2, camera1.descriptors.cols );

#ifdef DEBUG
          debug.push_back( dist );
#endif

          dists[i] = dist;
          if( dist < best )
          {
            best_i = i;
            best = dist;
          }
          count++;
        }

        // too few distances
        if ( count < 8 ) continue;

        // Non-max supression on dists
        int non_maxima_win = dists.size() / non_maxima_ratio;
        for ( int w = -non_maxima_win/2; w <= non_maxima_win/2; w++ )
        {
          if ( w == 0 ) continue;

          int ww = best_i + w;
          if ( ( ww < 0 ) ||
               ( ww >= (int)dists.size() )
             ) continue;

          dists[ww] = FLT_MAX;
        }

        // Find the second best
        second_best = FLT_MAX;
        for( size_t i=0; i < dists.size(); i++ )
        {
          if ( (int)i == best_i ) continue;
          if ( dists[i] < second_best )
            second_best = dists[i];
        }

        float ratio = best / second_best;

#ifdef DEBUG
          printf("epi.size():%lu debug.size():%lu x:%i y:%i best:[%f] best_i:%i, ratio:%.2f < (%.2f)\n", epi.size(), debug.size(), x, y, best, best_i, best / second_best, match_ratio_threshold);

          int range[2] = { 0, 0 };
          range[1] = debug.size();
          cv::Mat lineGraph = plotGraph( debug, range, best, second_best );

          namedWindow( "Display window", WINDOW_AUTOSIZE );
          imshow( "Display window", lineGraph );
          waitKey(0);
#endif


        if ( ratio < match_ratio_threshold )
        {
          // setup
          bool better = false;

          // SUBSAMPLE
          // near minima
          if ( true )
          {
//            printf("DBG d1[%f]->d2[%f] size:%lu [%f]->[%f]<-[%.40f] {%f} x:%i y%i\n",
//                   min_depth, max_depth, epi.size(), depths[best_i-1], depths[best_i], depths[best_i+1], best, (int)epi[best_i].x, (int)epi[best_i].y);

            epi.clear();
            float min_rdepth = depths[best_i-1];
            float max_rdepth = depths[best_i+1];


            float delta_depth = 0.01f;

            // numerical bounds: if goo too deep.
            if (max_rdepth == 0) max_rdepth = depths[best_i] + delta_depth;
            if (min_rdepth == 0) min_rdepth = depths[best_i] - delta_depth;

            // numerical bounds: if go too far.
            if ( depths[best_i] / min_rdepth > 10.0f ) min_rdepth = depths[best_i] - delta_depth;
            if ( max_rdepth / depths[best_i] > 10.0f ) max_rdepth = depths[best_i] + delta_depth;
            if ( depths[best_i] > max_rdepth ) max_rdepth = depths[best_i] + delta_depth;
            if ( depths[best_i] < min_rdepth ) min_rdepth = depths[best_i] - delta_depth;

//             printf("    d1[%f]->d2[%f] size:%lu [%f]->[%f]<-[%f] {%f}\n",
//              min_rdepth, max_rdepth, epi.size(), depths[best_i-1], depths[best_i], depths[best_i+1], best);

            depths.clear();

            /*
             * (2) Second SubSampling given img1 (x,y).
             * min_rdepth = near found minima minus one.
             * max_rdepth = near found minima plus one.
             * TODO: better/relaxed expression of this.
             */
            UniformSampling( camera1, camera2, ref, min_rdepth, max_rdepth, 0.5f, epi, Xs, depths );

            if ( epi.size() > 50000 ) continue;
            if ( epi.empty() ) continue;

            best_i = -1;
            dists.resize( epi.size(), FLT_MAX );
            for ( size_t i=0; i < epi.size(); i++ )
            {
              double xv = epi[i].x;
              double yv = epi[i].y;

              if ( (xv < border) ||
                   (xv >= camera2.m_image.cols - border) ||
                   (yv < border) ||
                   (yv >= camera2.m_image.rows - border)
                 ) continue;

              float desc3[ camera2.descriptors.cols ];
              camera2.m_daisy->GetDescriptor( yv, xv, 0, desc3 );

              // L1 or L2?
              float dist = FastSumSquaredDifference( desc1, desc3, camera1.descriptors.cols );

              dists[i] = dist;
              if( dist < best )
              {
                best_i = i;
                best = dist;
                better = true;
              }
            }

            //if (better)
              //printf("    d1[%f]->d2[%f] size:%lu [%f]->[%f]<-[%f] {%f} x:%f y:%f\n",
              // min_rdepth, max_rdepth, epi.size(), depths[best_i-1], depths[best_i], depths[best_i+1], best, epi[best_i].x, epi[best_i].y);

          } // END second refine

          // Found worse minima
          // when subsampling ?
          // Then skip this point.
          if ( ! better ) continue;


          int x1 = x;
          int y1 = y;
          int idx1 = y1*(camera1.m_image.cols) + x1;

          int x2 = epi[best_i].x;
          int y2 = epi[best_i].y;
          int idx2 = y2*(camera2.m_image.cols) + x2;

          float best_depth = depths[best_i];

          // Compute 3D point position
          Matx31f xx;
          xx(0) = x; xx(1) = y; xx(2) = 1;

          Matx31f X = best_depth * R0t_Z_invK0 * xx + camera1.Position();

          Point3f cur_pt(X(0), X(1), X(2));

          // Add new point
          if ( camera1.m_index_to_point_cloud[idx1] == -1 )
          {
            #pragma omp critical
            {
               PointRGB pt;
               pt.x = cur_pt.x;
               pt.y = cur_pt.y;
               pt.z = cur_pt.z;

               Vec3b bgrPixel = camera1.m_image.at<Vec3b>(y, x);
               // asign colors
               pt.r = bgrPixel.val[2];
               pt.g = bgrPixel.val[1];
               pt.b = bgrPixel.val[0];

               point_cloud_hit.push_back(PointRGB_hit_type(pt, 1));

               camera1.m_index_to_point_cloud[idx1] = point_cloud_hit.size() - 1;
               camera2.m_index_to_point_cloud[idx2] = point_cloud_hit.size() - 1;

             }
          }
          else
          {
            // See if this 3d point is close to the previous match
            int pt_idx = camera1.m_index_to_point_cloud[idx1];

            Point3f &prev_pt = point_cloud_hit[pt_idx].first;

            Point3f diff = prev_pt - cur_pt;
            float dist = diff.x*diff.x + diff.y*diff.y + diff.z*diff.z;

            #pragma omp critical
            {
              if ( dist < max_dist_sq )
              {
                camera2.m_index_to_point_cloud[idx2] = pt_idx;
                point_cloud_hit[pt_idx].second++;
              }
              else
              {
                // No match.
                // Don't know what to do yet.
              }
            }
          }
        }
        else
        {
          // ration threshold failed.
        }
      }
    }
    nLastTick = TermProgress( 1.0f, nLastTick );
}

} // end emvs namespace
