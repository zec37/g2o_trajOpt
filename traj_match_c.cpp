/////////////////////////////////
/////// name: g2o_trajOpt ///////
/////// status: building ////////
/////////////////////////////////

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "/usr/include/boost/algorithm/string.hpp"	// temp usage.

#include <Eigen/Core>
#include <Eigen/Dense>
#include "g2o/types/slam2d/vertex_point_xy.h"
#include "g2o/types/slam2d/edge_pointxy.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
//#include "g2o/core/factory.h"
//#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/solvers/cholmod/linear_solver_cholmod.h"
#include "g2o/solvers/eigen/linear_solver_eigen.h"
#include "g2o/types/slam2d/edge_se2_pointxy_offset.h"

using namespace std;
using namespace Eigen;


/////////////////////////
/// Pre_define param. ///
const int tr_num = 80;
const int loops_each_traj = 6;

int main(int argc, char **argv) {
  

    /////////////////////////////////////
    /// Load all trajectories data. ///
    vector<Vector2d> poses;
    int *pose_sum = new int[tr_num];
    int *traj_size = new int[tr_num];

    int count = 0;
    int diff = 0;
    for(int i = 0; i < tr_num; i++){

        stringstream intostr;
        string traj_name, traj_index;
        intostr << (i + 1);
        intostr >> traj_index;
        traj_name = "../data/traj" + traj_index + ".txt";
        ifstream f(traj_name);

        while(!f.eof()){
            int index = 0;
            float x = 0, y = 0;
            f >> index >> x >> y;
            if(!f.good() || index == 0) break;
            poses.push_back(Vector2d(x, y));
            count ++;
        }
        if(i > 0) diff = count - pose_sum[i - 1];
        else diff = count;
        *(pose_sum + i) = count;
        *(traj_size + i) = diff;

        cout << "Traj pose number: " << diff << endl;
        cout << "Total pose number: " << count << endl;
    }


    ///////////////////////
    /// Load loop data. ///
    string message;
    ifstream loop_file("../data/loops.txt");

      // Stored data. //
    vector<int> loopIndex;      // Store loop points.
    vector<vector<int>> globeIndexList;   // Store loop points' corresponding global ID in each trajectory.
      // Temp data. //
    vector<string> temp, lmGTID, lmIndex;   // temp loop points & global ID container.
    int loop_check[loops_each_traj];	// to find the reference of loop points in loopIndex.
    int tr_n = 0;	// to check current reading trajactory.

    while(getline(loop_file, message)){

        boost::split(temp, message, boost::is_any_of(":"));
        boost::split(lmGTID, temp[2], boost::is_any_of(" "));
        boost::split(lmIndex, temp[3], boost::is_any_of(" "));

        if(lmGTID.size() != loops_each_traj || lmIndex.size() != loops_each_traj)
        cerr << "Insufficient loop points." << endl;

        // Construct loop point list. //
        for(int i = 0; i < loops_each_traj; i++){

	    int int_ID = atoi(lmGTID[i].c_str());
	    
	    // When starting, the loopIndex is blank. //
	    if(loopIndex.empty()){		
	      loopIndex.push_back(int_ID);
              vector<int> temp_vec;	// Need to initialize global index as well.
              globeIndexList.push_back(temp_vec);
	    }
	    // Assign new index if found. //
	    else{
	      vector<int>::iterator iter = find(
		loopIndex.begin(), loopIndex.end(), int_ID);
	      if(iter == loopIndex.end()){
		  loopIndex.push_back(int_ID);
		  vector<int> temp_vec;	// Need to initialize global index as well.
		  globeIndexList.push_back(temp_vec);
	      }
	    }
            loop_check[i] = int_ID;
        }
        //copy(loopIndex.begin(), loopIndex.end(), ostream_iterator<int>(cout, "\n"));	//for debugging.

        // Assign global ID to the corresponding loop point. //
        for(int i = 0; i < loops_each_traj; i++){
            
	    int loop_ref = 0;
	    while(loopIndex[loop_ref] != loop_check[i]) loop_ref++;
	    int int_ID =  atoi(lmIndex[i].c_str());
            int globeID = (*(pose_sum + tr_n)) - (*(traj_size + tr_n)) + int_ID - 1;
	    globeIndexList[loop_ref].push_back(globeID);

            vector<int>::iterator end = globeIndexList[loop_ref].end();
            cerr << "The " << i + 1 << "th loop point for lmGTID for traj " << tr_n + 1 << " is " << loopIndex[loop_ref]  << endl;
            cerr << "The " << i + 1 << "th loop point for lmIndex for traj " << tr_n + 1 << " is " << *(end - 1) << endl;
        }
        tr_n++;
    }
    loop_file.close();


     ////////////////////////////////////
     /// Create Optimization Problem. ///
     g2o::SparseOptimizer optimizer;
     g2o::BlockSolverX::LinearSolverType *linearSolver;
     linearSolver = new g2o::LinearSolverEigen<g2o::BlockSolverX::PoseMatrixType>();
     g2o::BlockSolverX *solver_ptr = new g2o::BlockSolverX( linearSolver );
     g2o::OptimizationAlgorithmLevenberg *solver = new g2o::OptimizationAlgorithmLevenberg( solver_ptr );
     optimizer.setAlgorithm( solver );
     cout << "Solver initiated. " << endl;


    /////////////////////////////////
    /// Add odometry constraints. ///
    cout << "Odometry Optimization... " << endl;
    Matrix2d eye2d; eye2d << 1, 0, 0, 1;

    int vertice_total = 0;
    int edge_total = 0;
    int *gap_ptr = pose_sum;

    for(int i = 0; i < poses.size(); i++){

        /// Add vertex. ///

        g2o::VertexPointXY *vertice = new g2o::VertexPointXY();
        vertice->setId(i);
        vertice->setEstimate(poses[i]);
        optimizer.addVertex(vertice);
        vertice_total++;
    }

        /// Add edges. //

    for(int i = 0; i < poses.size(); i++){
        if(i != (*gap_ptr) - 1){
	  g2o::EdgePointXY *edge = new g2o::EdgePointXY();
	  edge->vertices()[0] = optimizer.vertex(i);
	  edge->vertices()[1] = optimizer.vertex(i + 1);
	  edge->setMeasurement(poses[i + 1] - poses[i]);
	  edge->setInformation(eye2d);
	  optimizer.addEdge(edge);
	  edge_total++;
	}
	else gap_ptr++;
    }
    
    vector<vector<int>>::iterator globe_iter = globeIndexList.begin() - 1;
    while(globe_iter++ != globeIndexList.end()){

       int list_size = (*globe_iter).size();
       for(int i = 0; i < list_size - 1; i++){
           for(int j = i + 1; j < list_size; j++){

               g2o::EdgePointXY *edge = new g2o::EdgePointXY();
               edge->vertices()[0] = optimizer.vertex((*globe_iter)[i]);
               edge->vertices()[1] = optimizer.vertex((*globe_iter)[j]);
               edge->setMeasurement(Vector2d(0, 0));
               edge->setInformation(eye2d);
               optimizer.addEdge(edge);
               edge_total++;
           }
       }
   }


    /////////////////////
    /// Optimization. ///

    cout << vertice_total << " vertices and " << edge_total << " edges" << " has been added to the graph." << endl;

    optimizer.save("loaded.g2o");
    cout << "Initialize optimization..." << endl;
    optimizer.setVerbose(true);
    optimizer.initializeOptimization();
    cout << "Running optimization..." << endl;
    optimizer.optimize(10);

    cout << "Saving optimization results..." << endl;
    optimizer.save("result.g2o");

    delete traj_size, pose_sum, gap_ptr;

    return 0;
}
