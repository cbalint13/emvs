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

#include <map>
#include <vector>
#include <iostream>
#include <fstream>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "emvs.hpp"

#define MIN_HIT 3


using namespace std;
using namespace emvs;

int main( int argc, char **argv )
{
    // PMVS dir path
    string PMVSDir;
    string OutPLYFilename;

    bool help = false;

    // parse arguments
    for ( int i = 1; i < argc; i++ )
    {
      if ( argv[i][0] != '-' )
      {
         if ( PMVSDir.empty() )
         {
            PMVSDir = argv[i];
            continue;
         }
         if ( OutPLYFilename.empty() )
         {
            OutPLYFilename = argv[i];
            continue;
         }
      }
    }

    if ( ( PMVSDir.empty() ) ||
         ( OutPLYFilename.empty() ) )
      help = true;

    if ( help )
    {
        cout << endl;
        cout << "Usage: emvs path_to_pmvs_dir out_ply_filename" << endl;
        cout << endl;
        exit( 1 );
    }
    cout << "PMVSDir: " << PMVSDir
         << "OutPLYFilename: " << OutPLYFilename << endl;

    double frequency = getTickFrequency();

    // scene cameras
    vector<Camera> cameras;

    // map of camarea pairs
    multimap<uint,uint> visdata;

    // final point cloud, with hit count
    // (how many cameras seen this 3d point)
    vector<PointRGB_hit_type> point_cloud_hit;

    // load bundler, but it doesn't
    // keep info of the image size
    LoadBundler( PMVSDir + "/bundle.rd.out", cameras );

    // load vis data
    LoadVisData( PMVSDir + "/vis.dat", cameras, visdata );

    // sequential processing
    int visnum = visdata.size();
    for ( size_t idx1 = 0; idx1 < cameras.size(); idx1++ )
    {

        pair<multimap<uint,uint>::iterator,
             multimap<uint,uint>::iterator> ret;

        if ( cameras[idx1].m_daisy_init == false )
        {
            char str[1024];

            sprintf(str, "/visualize/%08lu.jpg", idx1);
            cameras[idx1].m_image = imread( PMVSDir + "/" + str );
            Mat img = cameras[idx1].m_image;

            CV_Assert( img.data );

            // set missing width and height instrinsic
            cameras[idx1].SetIntrinsic( cameras[idx1].GetK()(0,0), img.cols, img.rows );
            flip( img, img, 0 );

            cameras[idx1].m_daisy = DAISY::create( 8, 2, 4, 4, DAISY::NRM_FULL, noArray(), true, false );
            //cameras[idx1].m_daisy = DAISY::create( 15, 3, 8, 8, DAISY::NRM_FULL, noArray(), true, false );
            cameras[idx1].m_daisy->compute( img, cameras[idx1].descriptors );

            cameras[idx1].m_daisy_init = true;
        }

        // feed pair cameras for idx1
        ret = visdata.equal_range(idx1);

        multimap<uint,uint>::iterator it;
        for ( it = ret.first; it != ret.second; ++it )
        {
            char str[128];

            int idx2 = it->second;

            cout << "Stereo pair: #" << visnum << " " << idx1 << " <-> " << idx2 << endl;

            if ( cameras[idx2].m_daisy_init == false )
            {
                sprintf(str, "/visualize/%08d.jpg", idx2);
                cameras[idx2].m_image = imread( PMVSDir + "/" + str );
                Mat img = cameras[idx2].m_image;

                CV_Assert( img.data );

                cameras[idx2].SetIntrinsic( cameras[idx2].GetK()(0,0), img.cols, img.rows );

                flip( img, img, 0 );

                cameras[idx2].m_daisy = DAISY::create( 8, 2, 4, 4, DAISY::NRM_FULL, noArray(), true, false );
                //cameras[idx2].m_daisy = DAISY::create( 15, 3, 8, 8, DAISY::NRM_FULL, noArray(), true, false );
                cameras[idx2].m_daisy->compute( img, cameras[idx2].descriptors );

                cameras[idx2].m_daisy_init = true;
            }

            // 3D point cloud from stereo pairs
            int64 pairStartTime = getTickCount();
            DoStereoPair( cameras[idx1], cameras[idx2], point_cloud_hit );
            int64 pairEndTime = getTickCount();
            printf("Compute Time: %.4f sec\n", ( pairEndTime - pairStartTime ) / frequency );

            cout << endl;

            // idx2 done, clean unused memory
            cameras[idx2].m_daisy.release();
            cameras[idx2].descriptors.release();
            cameras[idx2].m_daisy_init = false;

            // periodically write to a .PLY file,
            // so we can see what it's producing
            OutputPLY( OutPLYFilename, point_cloud_hit, MIN_HIT );

            visnum--;
        }
        // idx1 done, clean unused memory
        cameras[idx1].m_daisy.release();
        cameras[idx1].descriptors.release();
        cameras[idx1].m_daisy_init = false;
    }

    return 0;
}
