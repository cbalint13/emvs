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


#include <iostream>
#include <fstream>

#include "emvs.hpp"


using namespace std;
using namespace emvs;

namespace emvs {

void LoadBundler( const string& bundle_file, vector<Camera>& cameras )
{
    ifstream bundle(bundle_file.c_str());

    if (!bundle) {
      cerr << "LoadBundle: Error opening " << bundle_file << " for reading" << endl;
      CV_Assert(bundle);
    }

    int num;
    char line[1024];
    stringstream str;

    // header
    bundle.getline(line, sizeof(line));
    bundle.getline(line, sizeof(line));

    str.str(line);
    str >> num;

    cameras.resize(num);

    for (int i=0; i < num; i++)
    {
      Matx33f R;
      Matx31f t;

      float focal;

      // focal, r1, r2
      bundle.getline(line, sizeof(line));

      str.str(line);
      str.clear();
      str >> focal;

      // image dimensions
      // are missing for now
      // will be set to -1
      cameras[i].SetIntrinsic( focal, -1, -1 );

      // rotation #1
      bundle.getline(line, sizeof(line));
      str.str(line); str.clear();
      str >> R(0,0); str >> R(0,1); str >> R(0,2);
      // rotation #2
      bundle.getline(line, sizeof(line));
      str.str(line); str.clear();
      str >> R(1,0); str >> R(1,1); str >> R(1,2);
      // rotation #3
      bundle.getline(line, sizeof(line));
      str.str(line); str.clear();
      str >> R(2,0); str >> R(2,1); str >> R(2,2);
      // translation
      bundle.getline(line, sizeof(line));
      str.str(line); str.clear();
      str >> t(0); str >> t(1); str >> t(2);

      cameras[i].SetRotation(R);
      cameras[i].SetTranslation(t);
    }
    bundle.close();
}

void LoadVisData( const string& visdata_file, const vector<Camera>& cameras, multimap<uint,uint>& visdata )
{
    ifstream visfile( visdata_file.c_str() );

    if ( !visfile ) {
      cerr << "LoadVisData: Error opening " << visdata_file << " for reading" << endl;
      CV_Assert( visfile );
    }

    int num;
    char line[1024];
    stringstream str;

    // header
    visfile.getline(line, sizeof(line));
    visfile.getline(line, sizeof(line));

    str.str(line);
    str >> num;


    for (int i=0; i < num; i++)
    {

      int idx, count;

      // focal, r1, r2
      visfile.getline(line, sizeof(line));


      str.str(line);
      str.clear();
      str >> idx;
      str >> count;

      cout << "\e[0m" << idx << ":";
      for (int c=0; c<count; c++)
      {
        int idp;
        str >> idp;
        float norm;

        Matx41f ray0, ray1;

        ray0(0) = cameras[idx].GetR()(2,0);
        ray0(1) = cameras[idx].GetR()(2,1);
        ray0(2) = cameras[idx].GetR()(2,2);
        ray0(3) = 0.0f;

        ray1(0) = cameras[idp].GetR()(2,0);
        ray1(1) = cameras[idp].GetR()(2,1);
        ray1(2) = cameras[idp].GetR()(2,2);
        ray1(3) = 0.0f;

        norm = ray0(0)*ray0(0) + ray0(1)*ray0(1)
             + ray0(2)*ray0(2) + ray0(3)*ray0(3);
        ray0 *= 1/sqrt(norm);

        norm = ray1(0)*ray1(0) + ray1(1)*ray1(1)
             + ray1(2)*ray1(2) + ray1(3)*ray1(3);
        ray1 *= 1/sqrt(norm);

        float angle = ray0(0) * ray1(0) + ray0(1) * ray1(1)
                    + ray0(2) * ray1(2) + ray0(3) * ray1(3);
        angle = acos(angle) * 180.0f/CV_PI;

        // 10<->30 angle
        // 0.8<>1.2 scale
        if ( (angle <= 30) &&
             (angle >= 10) )
        {
          // good basleine pair
          cout << " \e[0;32m" << idp;
          visdata.insert( pair<int,int>(idx,idp) );
        } else
        {
          // good basleine pair
          cout << " \e[0;31m" << idp;
        }
      }
      cout << "\e[0m" << endl;
    }

    cout << "Total pairs: " << visdata.size() << endl;

    visfile.close();
}

} // end emvs namespace
