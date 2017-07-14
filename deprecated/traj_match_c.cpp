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
#include "g2o/types/slam2d/vertex_se2.h"
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


/////////////////////////////////////////////////////////////////////////
/// The function to calculate tranform matrix between two trajactory. ///
Matrix3d calculate_transform(Vector2d v1n, Vector2d v2n, Vector2d v1o, Vector2d v2o){

    Vector2d _v1n = v1n, _v2n = v2n, _v1o = v1o, _v2o = v2o;
    double _cos, _sin, _tx, _ty;

    Matrix3d _coeff;
    Vector3d _result;
    Vector3d _solver;

    // v1 means the first cross point. v2 means the second cross point.
    // (0) means vector.x  ||  (1) means vector.y
    // n refers to new, means the new traj that need to be tranformed.
    // o refers to old, means the old traj that should be tranformed to.
    _coeff <<  
        _v1n(0) * _v1n(0) + _v1n(1) * _v1n(1),      _v1n(0),     _v1n(1),
        _v2n(0) * _v2n(0) + _v2n(1) * _v2n(1),      _v2n(0),     _v2n(1),
        _v1n(0) * _v2n(0) + _v1n(1) * _v2n(1),      _v1n(0),     _v2n(1);

    _result <<  
        _v1n(0) * _v1o(0) + _v1n(1) * _v1o(1),
        _v2n(0) * _v2o(0) + _v2n(1) * _v2o(1),
        _v1n(0) * _v2o(0) + _v1o(1) * _v2n(1);

    _solver = _coeff.colPivHouseholderQr().solve(_result);

    _cos = _solver(0);
    _tx = _solver(1);
    _ty = _solver(2);

    _sin = ( -_v1o(0) + _v1n(0) * _cos + _tx) / _v1n(1);

    cerr << "cos= " << _cos << endl;
    cerr << "tx= " << _tx << endl;
    cerr << "ty= " << _ty << endl;

    Matrix3d _tranform;
    _tranform << 
        _cos, -_sin, 	_tx,
        _sin, 	_cos, 	_ty,
        0, 	   0,	  1;

    return _tranform;
}


/////////////////////////////////////////////////////
/// This function is used to return the reference ///
/// of all combinations in a container. /////////////
vector<int> combine2(int top, int size, vector<int> memory){

    if(size <= 1) { return memory; }

    for(int i = 0; i < (size - 1); i++){

        memory.push_back( i / (size - 1) + top - size);
        memory.push_back( i % (size - 1) + 1 + top - size);
    }
    memory = combine2(top, size - 1, memory);
    return memory;
}


/////////////////////////
/// Pre_define param. ///
const int tr_num = 8;
const int loops_each_traj = 6;

int main(int argc, char **argv) {


    /////////////////////////////////////
    /// Load all trajectories data. ///
    vector<Vector2d> poses;
    int* pose_sum = new int[tr_num];
    int* traj_size = new int[tr_num];

    int count = 0;
    int diff = 0;
    for(int i = 0; i < tr_num; i++){

        stringstream intostr;
        string traj_name, traj_index;
        intostr << (i * 10 + 1);    // Need to fix for different read.
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
    vector<vector<int>> globeIndexList;   // Store loop points' corresponding global IDs in all trajectories.
    int loop_of_traj[tr_num][loops_each_traj];  	// Store loop points in each traj.
    int loopref_of_traj[tr_num][loops_each_traj]; 	// Store loop reference in each traj.

    // Temp data. //
    vector<string> temp, lmGTID, lmIndex;   // temp loop points & global ID container.
int loop_check[loops_each_traj];	// to find the reference of loop points in loopIndex.
int tr_n = 0;	// to check current reading trajactory.
int lookup = 0; 	// to check if current traj is the 1st among 10.

    while(getline(loop_file, message)){

        // Only read 1st traj in 10. //
        if( (lookup++ % 10) != 0){
            continue;
        }

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
        //   copy(loopIndex.begin(), loopIndex.end(),
        //        ostream_iterator<int>(cout, "\n"));	//for debugging.

        // Assign global ID to the corresponding loop point. //
        for(int i = 0; i < loops_each_traj; i++){

            int loop_ref = 0;
            while(loopIndex[loop_ref] != loop_check[i]) loop_ref++;
            int int_ID =  atoi(lmIndex[i].c_str());
            int globeID = (*(pose_sum + tr_n)) - (*(traj_size + tr_n)) + int_ID - 1;
            globeIndexList[loop_ref].push_back(globeID);

            vector<int>::iterator end = globeIndexList[loop_ref].end();
            cerr << "The " << i + 1 << "th loop point for lmGTID for traj " << tr_n + 1 << " is " << loopIndex[loop_ref] << endl;
            cerr << "The " << i + 1 << "th loop point for lmIndex for traj " << tr_n + 1 << " is " << *(end - 1) << endl;
        }

        for(int i = 0; i < loops_each_traj; i++){
            loop_of_traj[tr_n][i] = atoi(lmGTID[i].c_str());
            loopref_of_traj[tr_n][i] = atoi(lmIndex[i].c_str());
        }
        tr_n++;
        }
        loop_file.close();


    //////////////////////////////////////////////
    /// Move trajs to fit in the global scene. ///
    Matrix3d* trans_matrices = new Matrix3d[tr_num];
    vector<int> traj_comb;
    traj_comb = combine2(tr_num, tr_num, traj_comb);
    vector<int> moved;
    moved.push_back(0);     // Initiate the first traj.

    int match_count = 0;
    int traj_old = 0, traj_new = 0;
    for(int i = 0; i < (tr_num * (tr_num - 1) / 2); i++){

        traj_old = traj_comb[2 * i];
        traj_new = traj_comb[2 * i + 1];
        Vector2d buf_old[2];
        Vector2d buf_new[2];

        // Search two points that exists in both trajs. //
        int search = 2;
        int loop_old = 0, loop_new = 0;
        int ref_old = 0, ref_new = 0;
        for(int j = 0; j < (loops_each_traj * loops_each_traj); j++){

            loop_old = j / loops_each_traj;
            loop_new = j % loops_each_traj;

            if(loop_of_traj[traj_old][loop_old] == loop_of_traj[traj_new][loop_new] && search != 0){

                ref_old = loopref_of_traj[traj_old][loop_old] + pose_sum[traj_old] - traj_size[traj_old];
                buf_old[2 - search] = poses[ref_old];
                ref_new = loopref_of_traj[traj_new][loop_new] + pose_sum[traj_new] - traj_size[traj_new];
                buf_new[2 - search] = poses[ref_new];
                search--;
            }
        }

        // Use the matched point pair to calibrate. //
        Matrix3d transform;
        int traj_start = pose_sum[traj_new] - traj_size[traj_new];
        int traj_end = pose_sum[traj_new];

        if(search == 0) {

            cerr << "traj" << traj_old + 1 << " and traj" << traj_new + 1 << " has match pairs " << endl;
            match_count++;

            vector<int>::iterator it_old;
            vector<int>::iterator it_new;
            it_old = find(moved.begin(), moved.end(), traj_old);
            it_new = find(moved.begin(), moved.end(), traj_new);
            if(it_new != moved.end()) { 
                continue;      // If the new new traj has already been moved.
            }
            else if(it_old != moved.end()){

                // The old traj has been moved, thus can be used to calibrate.
                transform = calculate_transform(
                    buf_new[0], buf_new[1], buf_old[0], buf_old[1]);
                for(int k =  traj_start; k < traj_end; k++){
                    Vector3d temp;
                    temp << poses[k], 1;
                    temp = transform * temp;
                    poses[k](0) = temp(0);
                    poses[k](1) = temp(1);
                }
                moved.push_back(traj_new);
            }
            else{
                cout << "Unable to calibrate now..." << endl;
            }
        }

    }
    cerr << match_count << " pairs found." << endl;


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
    int* gap_ptr = pose_sum;

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

    //////////////////////////////////
    /// Cross points Optimization. ///
    // check input. //
    int loopList_size = loopIndex.size();
    for(int i = 0; i < loopList_size; i++){
        cout << "loop point" << i + 1 << " is " << loopIndex[i] << endl;
        cout << "Global IDs are ";
        vector<int>::iterator iter;
        for(iter = globeIndexList[i].begin(); iter != globeIndexList[i].end(); iter++){
            cout << *iter << " ";
        }
        cout << endl;
    }
    for(int i = 0; i < tr_num; i++){
        cout << "Loops in traj" << i + 1 << " are ";
        for(int j = 0; j < loops_each_traj; j++){
            cout << loop_of_traj[i][j] << " ";
        }
        cout << endl;
        cout << "Global id in traj" << i + 1 << " are ";
        for(int j = 0; j < loops_each_traj; j++){
            cout << loopref_of_traj[i][j] << " ";
        }
        cout << endl;
    }
    // Add edge. //
    vector<vector<int>>::iterator globe_iter = globeIndexList.begin() - 1;
    while(globe_iter++ != globeIndexList.end()){

        int list_size = (*globe_iter).size();
        vector<int> globeList_comb;
        globeList_comb = combine2(list_size, list_size, globeList_comb);
        if(globeList_comb.size() != 0){
	    cerr << globeList_comb.size() << endl;
            for(int i = 0; i < (list_size * (list_size - 1) / 2); i++){

                int j = (*globe_iter)[globeList_comb[2 * i]];
                int k = (*globe_iter)[globeList_comb[2 * i + 1]];
		
                g2o::EdgePointXY *edge = new g2o::EdgePointXY();
                edge->vertices()[0] = optimizer.vertex(j);
                edge->vertices()[1] = optimizer.vertex(k);
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

    delete traj_size, pose_sum, gap_ptr, loop_of_traj, trans_matrices;

    return 0;
}