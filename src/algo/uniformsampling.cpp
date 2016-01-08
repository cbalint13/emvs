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

 "A Closed-Form Solution for the Uniform Sampling of the Epipolar Line
 via Non-Uniform Depth Sampling" Tola, Engin EPFL-REPORT-150161

 */

#include <iostream>


#include "emvs.hpp"


using namespace std;
using namespace emvs;

namespace emvs {

void UniformSampling( const Camera& camera1, const Camera& camera2, const Point2i& Img2DPoint,
                      float Lmin, float Lmax, float r, vector<Point2f>& epi,
                      vector<Point3f>& Xs, vector<float>& depths )
{

    // Magic matrix that makes
    // everything work with bundler
    Matx33f Z = Matx33f::eye(); Z(2,2) = -1;

    Matx31f x; x(2) = 1.0f;
    x(0) = Img2DPoint.x; x(1) = Img2DPoint.y;

    Matx33f K1 = camera2.GetK();
    Matx33f inv_K0 = camera1.GetInvK();

    Matx33f R0 = camera1.GetR();
    Matx33f R1 = camera2.GetR();

    Matx31f C0 = camera1.Position();

    /* not used now
      Matx31f C1 = camera2.Position();
    */

    Matx31f t1 = camera2.GetT();

    Matx31f a = K1*Z*R1*R0.t()*Z*inv_K0*x;
    Matx31f b = K1*Z*R1*C0 + K1*Z*t1;

    int width = camera2.m_image.cols;
    int height = camera2.m_image.rows;

    // should be renamed to start, end
    // as umin does not mean it is always
    // smaller than umax

    float umin = (Lmin*a(0) + b(0)) / (Lmin*a(2) + b(2));
    float vmin = (Lmin*a(1) + b(1)) / (Lmin*a(2) + b(2));

    float umax = (Lmax*a(0) + b(0)) / (Lmax*a(2) + b(2));
    float vmax = (Lmax*a(1) + b(1)) / (Lmax*a(2) + b(2));

    // protect against inf values
    if ( Lmax > 1000000000 )
    {
      epi.clear();
      depths.clear();
      Xs.clear();
      return;
    }

    // Re-calc and modify Lmin
    float orig_Lmin = Lmin;
    float orig_Lmax = Lmax;

    /*
     * omptimize
     * quadrants
     */

    float _Lmin = 0, _Lmax = 0;

    if ( umin < 0 ) _Lmin = -b(0) / a(0);
    if ( vmin < 0 ) _Lmin = -b(1) / a(1);
    if ( umin >= width )
      _Lmin = (b(0) - (width -1)*b(2)) / ((width -1)*a(2) - a(0));
    if ( vmin >= height )
      _Lmin = (b(0) - (height-1)*b(2)) / ((height-1)*a(2) - a(0));


    if ( _Lmin >= orig_Lmin && _Lmin <= orig_Lmax )
    {
        Lmin = _Lmin;
        umin = (Lmin*a(0) + b(0)) / (Lmin*a(2) + b(2));
        vmin = (Lmin*a(1) + b(1)) / (Lmin*a(2) + b(2));
    }


    if ( umax < 0 ) _Lmax = -b(0) / a(0);
    if ( vmax < 0 ) _Lmax = -b(1) / a(1);
    if ( umax >= width )
      _Lmax = (b(0) - (width -1)*b(2)) / ((width -1)*a(2) - a(0));
    if ( vmax >= height )
      _Lmax = (b(0) - (height-1)*b(2)) / ((height-1)*a(2) - a(0));

    if (_Lmax >= orig_Lmin && _Lmax <= orig_Lmax)
    {
        Lmax = _Lmax;
        umax = (Lmax*a(0) + b(0)) / (Lmax*a(2) + b(2));
        vmax = (Lmax*a(1) + b(1)) / (Lmax*a(2) + b(2));
    }


    float Du = umax-umin;
    float Dv = vmax-vmin;

    float dl = sqrt(Du*Du + Dv*Dv);
    float du = Du / dl; float dv = Dv / dl;

    float wmin = Lmin*a(2) + b(2);

    float w = wmin; float L = Lmin;
    float u = umin; float v = vmin;

    /* not used now
     * Matx31f R0t_invK0_x = R0.t()*Z*inv_K0*x;
     */

    // reserve storage
    if ( (int)(dl+0.5f) > 0 )
      epi.reserve( (int)(dl+0.5f) );
    else
    {
      epi.clear();
      depths.clear();
      Xs.clear();
      return;
    }

    /* not used now
     * Xs.reserve( (int)dl );
     */

    depths.reserve((int)(dl+0.5f));
    for ( float length=0; length <= dl; length += r )
    {
      float dL;
      float rdu = r*du;
      float rdv = r*dv;

      if ( fabs(du) != 0 )
        dL = w*rdu/(a(0) - a(2)*(u + rdu));
      else
        dL = w*rdv/(a(1) - a(2)*(v + rdv));


      L += dL;
      u += rdu;
      v += rdv;

      w = L*a(2) + b(2);

      // if out of image (should not happen)
      if ( u < 0 || v < 0 || u >= width || v >= height )
        continue;

      /* not used now
       * Matx31f X = L*R0t_invK0_x + C0;
       */

      // accumulate u,v
      epi.push_back( Point2f(u,v) );

      /* not used now
       * Xs.push_back(Point3f(X(0), X(1), X(2)));
       */

      // accumulate depth
      depths.push_back( L );
    }
}

} // end emvs namespace
