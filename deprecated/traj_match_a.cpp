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


/////////////////////////////////////////
/// Read trajectory from a data file. ///
vector<Vector2d> traj_read(char* filename){

    vector<Vector2d> traj;

    fstream f(filename);
    if(!f) cerr << "File not found.";
    else cout << "Reading from " << filename << endl;

    while(!f.eof()){

        int index = 0;
        float x = 0, y = 0;
        f >> index >> x >> y;
        if (!f.good() || index == 0)
        break;
        traj.push_back(Vector2d(x, y));
    }
    f.close();
    return traj;
}


// ///////////////////////////////////////////////////////////////////////
// /// Construct possible matching pairs for a single consensus point. ///
// int *pairs[2](int *id_list, int size){
// 
//     int temp[size * (size - 1)][2];
//     for(int i = 0; i < size - 1; i++){
//         for(int j = i + 1; j < size; j++){
//             
//         }
//     }
// }


///////////////////////////////////////////////
// ///// OpenGL Extension | not started ///////
// void init()
// {
//     glClearColor(0.0, 0.0, 0.0, 0.0);
//     glMatrixMode(GL_PROJECTION);
//     glLoadIdentity();
//     gluOrtho2D(-50.0, 50.0, -50.0, 50.0);
// }
// 
// void drawTraj()
// {
//     
//     glColor3d(1.0, 1.0, 1.0);
//     for(int i; i < traj1.size(); i++){
//       glVertex2d(traj1[i](0), traj1[i](1));
//     }
// }


int main(int argc, char **argv) {


    /////////////////////////////////////
    /// Load all trajectories data. ///
    vector<Vector2d> *trajs = new vector<Vector2d>[tr_num];
    unsigned int poses_total = 0;
    unsigned int *traj_size = new unsigned int[tr_num];
    unsigned int *pose_sum = new unsigned int[tr_num];
    unsigned int edge_total = 0;

    for(int i = 0; i < tr_num; i++){

        char *file = new char[20];
        int size = 0;

        sprintf(file, "../data/traj%d.txt", i + 1);
        *(trajs + i) = traj_read(file);
        size = (*(trajs + i)).size();

        cout << "Pose number: " << size << endl;
        poses_total += size;
        *(traj_size + i) = size;
        *(pose_sum + i) = poses_total;
        cout << "Summerized pose number: " << poses_total << endl;
        delete[] file;
    }    


    ///////////////////////
    /// Load loop data. ///
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
    /// Odometry optimization. ///
    cout << "Odometry Optimization... " << endl;

    int offset = 0;
    for(int i = 0; i < tr_num; i++){

        if(i > 0) offset = pose_sum[i-1];

        /// Add vertex. ///
        for(int j = 0; j < *(traj_size + i); j++){

            g2o::VertexPointXY *vertice = new g2o::VertexPointXY();
            vertice->setId(j + offset);
            vertice->setEstimate((*(trajs + i))[j]);
            optimizer.addVertex(vertice);
        }

        /// Add edges. ///
        for(int k = 0; k < *(traj_size + i) - 1; k++){

            g2o::EdgePointXY *edge = new g2o::EdgePointXY();
            edge->vertices()[0] = optimizer.vertex(k + offset);
            edge->vertices()[0] = optimizer.vertex(k + offset + 1);
            edge->setMeasurement((*(trajs + i))[k - 1] - (*(trajs + i))[k]);
            //edge->setInformation(Eigen::Matrix2d::Identity());
            optimizer.addEdge(edge);
            edge_total++;
        }      
    }

    //////////////////////////////
    /// Drifting optimization. ///
    

    /////////////////////
    /// Optimization. ///
    
    cout << poses_total << " vertices and " << edge_total << " edges" << " has been added to the graph." << endl;

    optimizer.save("loaded.g2o");
    cout << "Initialize optimization..." << endl;
    optimizer.setVerbose(true);
    optimizer.initializeOptimization();
    cout << "Running optimization..." << endl;
    optimizer.optimize(10);

    cout << "Saving optimization results..." << endl;
    optimizer.save("result.g2o");

    delete traj_size, pose_sum, lmGTID_all, lmIndex_all;

    return 0;
}
