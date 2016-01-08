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


#include <vector>
#include <iostream>
#include <fstream>

#include <opencv2/core/core.hpp>

#include "emvs.hpp"


using namespace cv;
using namespace std;

namespace emvs {

void OutputPLY(const string &filename, const vector<PointRGB_hit_type>& point_cloud_hit, uint min_hit)
{
    ofstream output( filename.c_str() );

    CV_Assert( output );

    printf( "Total Points: %lu\n", point_cloud_hit.size() );

    int count = 0;
    for( size_t j = 0; j < point_cloud_hit.size(); j++ )
    {
        uint hit = point_cloud_hit[j].second;
        if ( hit >= min_hit )
          count++;
    }

    printf("Valid Points: %i\n", count);

    output << "ply" << endl;
    output << "format ascii 1.0" << endl;
    output << "element face 0" << endl;
    output << "property list uchar int vertex_indices" << endl;
    output << "element vertex " << count << endl;
    output << "property float x" << endl;
    output << "property float y" << endl;
    output << "property float z" << endl;
    output << "property uchar diffuse_red" << endl;
    output << "property uchar diffuse_green" << endl;
    output << "property uchar diffuse_blue" << endl;
    output << "end_header" << endl;

    for ( size_t j = 0; j < point_cloud_hit.size(); j++ )
    {
      const PointRGB &pt = point_cloud_hit[j].first;
      uint hit = point_cloud_hit[j].second;

      if ( hit >= min_hit )
      {
        output << pt.x << " " << pt.y << " " << pt.z << " " << (int) pt.r << " " << (int) pt.g << " " << (int) pt.b << endl;
      }
    }
}

} // end emvs namespace
