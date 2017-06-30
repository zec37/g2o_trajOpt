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

//#include <Eigen/Core>
#include "g2o/types/slam2d/vertex_point_xy.h"
#include "g2o/types/slam2d/edge_pointxy.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
//#include "g2o/core/factory.h"
//#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/solvers/cholmod/linear_solver_cholmod.h"
//#include "g2o/solvers/eigen/linear_solver_eigen.h"

using namespace std;
using namespace Eigen;


/////////////////////////
/// Pre_define param. ///
const int tr_num = 80;
const int loop_point_num = 6;

int main(int argc, char **argv) {


    /////////////////////////////////////
    /// Load all trajectories data. ///
    vector<Vector2d> poses;
    unsigned int *pose_sum = new unsigned int[tr_num];
    unsigned int *traj_size = new unsigned int[tr_num];

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


    /////////////////////////////////////////////
    /// Load loop data as splitted segements. ///
    string message;
    ifstream loop_file("../data/loops.txt");
    vector<string> temp, lmGTID, lmIndex;
    vector<int> *lmGTID_all = new vector<int>[tr_num];
    vector<int> *lmIndex_all = new vector<int>[tr_num];

    int tr_n = 0;
    while(getline(loop_file, message)){

        boost::split(temp, message, boost::is_any_of(":"));
        boost::split(lmGTID, temp[2], boost::is_any_of(" "));
        boost::split(lmIndex, temp[3], boost::is_any_of(" "));

        if(lmGTID.size() != loop_point_num || lmIndex.size() != loop_point_num)
        cerr << "Insufficient loop points." << endl;

        for(int i = 0; i < loop_point_num; i++){

            (*(lmGTID_all + tr_n)).push_back(atoi(lmGTID[i].c_str()));
            (*(lmIndex_all + tr_n)).push_back(atoi(lmIndex[i].c_str()));    
            //cerr << "The " << i << "th loop point for lmGTID for traj " << tr_n + 1 << " is " << (*(lmGTID_all + tr_n))[i] << endl;
            //cerr << "The " << i << "th loop point for lmIndex for traj " << tr_n + 1 << " is " << (*(lmIndex_all + tr_n))[i] << endl;
        }
        tr_n++;
    }
    loop_file.close();


     ////////////////////////////////////
     /// Create Optimization Problem. ///
     g2o::SparseOptimizer optimizer;
     g2o::BlockSolverX::LinearSolverType *linearSolver;
     linearSolver = new g2o::LinearSolverCholmod<g2o::BlockSolverX::PoseMatrixType>();
     g2o::BlockSolverX *solver_ptr = new g2o::BlockSolverX( linearSolver );
     g2o::OptimizationAlgorithmLevenberg *solver = new g2o::OptimizationAlgorithmLevenberg( solver_ptr );
     optimizer.setAlgorithm( solver );
     cout << "Solver initiated. " << endl;


    /////////////////////////////////
    /// Add odometry constraints. ///
    cout << "Odometry Optimization... " << endl;

    unsigned int vertice_total = 0;
    unsigned int edge_total = 0;
    unsigned int *gap_ptr = pose_sum;

    for(int i = 0; i < poses.size(); i++){

        /// Add vertex. ///

        g2o::VertexPointXY *vertice = new g2o::VertexPointXY();
        vertice->setId(i);
        vertice->setEstimate(poses[i]);
        optimizer.addVertex(vertice);
        vertice_total++;

        /// Add edges. //

        if(i == *gap_ptr - 1) { gap_ptr++; continue;}

        g2o::EdgePointXY *edge = new g2o::EdgePointXY();
        edge->vertices()[0] = optimizer.vertex(i);
        edge->vertices()[1] = optimizer.vertex(i + 1);
        edge->setMeasurement(poses[i + 1] - poses[i]);
        //edge->setInformation(Eigen::Matrix2d::Identity());
        optimizer.addEdge(edge);
        edge_total++;
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

    delete traj_size, pose_sum, lmGTID_all, lmIndex_all, gap_ptr;

    return 0;
}
