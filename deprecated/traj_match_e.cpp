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
#include "g2o/types/slam2d/edge_se2_pointxy_offset.h"

#include "g2o/core/base_vertex.h"
#include "g2o/core/base_unary_edge.h"
//#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/core/robust_kernel_impl.h"

//#include "g2o/solvers/cholmod/linear_solver_cholmod.h"
#include "g2o/solvers/eigen/linear_solver_eigen.h"

using namespace std;
using namespace Eigen;


///////////////////////////////////////
/// Class for transform optimizing. ///
/// Vector4d = {s, theta, tx, ty} /////
class VertexSim2: public g2o::BaseVertex<4, Vector4d>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    VertexSim2() : BaseVertex<4, Vector4d> () {}
    virtual bool read ( istream& is )
    {
        return true;
    }
    virtual bool write ( ostream& os ) const
    {
        return true;
    }

    virtual void setToOriginImpl()
    {
        _estimate = Vector4d ( 1,0,0,0 );
    }

    virtual void oplusImpl ( const double* update_ )
    {
        Map<const Vector4d> update ( update_ );
        _estimate += update;
    }
};

class EdgeSim2: public g2o::BaseUnaryEdge<2, Vector2d, VertexSim2>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    EdgeSim2 ( const Vector2d origin ) : BaseUnaryEdge<2, Vector2d, VertexSim2> (), _origin ( origin ) {}
    virtual bool read ( istream& is )
    {
        return true;
    }
    virtual bool write ( ostream& os ) const
    {
        return true;
    }
    virtual void computeError()
    {
        VertexSim2* v = ( VertexSim2* ) this->vertices() [0];
        Vector4d estimate = v->estimate();
        Vector2d m = this->measurement();
        double s = estimate[0];
        double theta = estimate[1];
        double tx = estimate[2];
        double ty = estimate[3];
        double x = s* ( cos ( theta ) *m[0] - sin ( theta ) *m[1] ) + tx ;
        double y = s* ( sin ( theta ) *m[0] + cos ( theta ) *m[1] ) + ty ;
        _error = _origin - Vector2d ( x,y );
    }
private:
    Vector2d _origin;
};


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


// ------------------------------------------------- //
// ------------------------------------------------- //


/////////////////////////
/// Pre_define param. ///
const int tr_num = 80;
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
        intostr << (i + 1);    // Need to fix for different read.
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
      
	if(tr_n == tr_num) { break; }

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
                vector<int> temp_vec2;
            }
            // Assign new index if found. //
            else{
                vector<int>::iterator iter = find(
                    loopIndex.begin(), loopIndex.end(), int_ID);
                if(iter == loopIndex.end()){
                    loopIndex.push_back(int_ID);
                    vector<int> temp_vec;	// Need to initialize global index as well.
                    globeIndexList.push_back(temp_vec);
                    vector<int> temp_vec2;
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
        }

        for(int i = 0; i < loops_each_traj; i++){
            loop_of_traj[tr_n][i] = atoi(lmGTID[i].c_str());
            loopref_of_traj[tr_n][i] = atoi(lmIndex[i].c_str());
        }
        tr_n++;
        }
        loop_file.close();


// ------------------------------------------------------------------- //

    
    Matrix2d eye2d; eye2d << 1, 0, 0, 1;

    ///////////////////////////////////////////
    /// Optimize the tranform of each traj. ///

        // check input //
    for(int i = 0; i < tr_num; i++){
        cout << "Loops in traj" << i + 1 << " are ";
        for(int j = 0; j < loops_each_traj; j++){
            cout << loop_of_traj[i][j] << " ";
        }
        cout << endl;
        cout << "Local id in traj" << i + 1 << " are ";
        for(int j = 0; j < loops_each_traj; j++){
            cout << loopref_of_traj[i][j] << " ";
        }
        cout << endl;
    }
    

    for(int i = 0; i < tr_num; i++){
	
	if(i == 0) { 
	  
	  cout << "Pass this traj..." << endl;
	  continue;
	} 	// Skip the first traj.

	g2o::SparseOptimizer pre_optimizer;
	g2o::BlockSolverX::LinearSolverType *pre_linearSolver; 
	pre_linearSolver = new g2o::LinearSolverEigen<g2o::BlockSolverX::PoseMatrixType>();
	g2o::BlockSolverX* pre_solver_ptr = new g2o::BlockSolverX( pre_linearSolver );
	g2o::OptimizationAlgorithmLevenberg* pre_solver = new g2o::OptimizationAlgorithmLevenberg( pre_solver_ptr );
	pre_optimizer.setAlgorithm( pre_solver );

        cout << "Initialize pre-optimization for traj" << i + 1 << endl;
	
        int pre_edge_count = 0;
	
        // Add vertex. i.e. the transform matrix. //
        VertexSim2* v = new VertexSim2();
        v->setId(0);
        pre_optimizer.addVertex(v);

        // Add edges by loopIndex. //
        for(int j = 0; j < loops_each_traj * loopIndex.size(); j++){

            int ref_List = j / loops_each_traj;     
            // Current loop point to be compared in the all loop point list.
            int ref_This = j % loops_each_traj; 
            // Current loop point to compare in this traj.

            if(loopIndex[ref_List] == loop_of_traj[i][ref_This]){
            
                int ref_This_globe = loopref_of_traj[i][ref_This] + pose_sum[i] - traj_size[i] - 1;

		int read_check = 0;
                for(int k = 0; k < globeIndexList[ref_List].size(); k++){
		  
		    //if(read_check > 0) { break; }

                    int ref_List_globe = globeIndexList[ref_List][k];
		    
                    if(ref_List_globe >= ref_This_globe){
                        break;      // Omit those should not be loaded yet.
                    }
                    
                    EdgeSim2* transform_edge = new EdgeSim2( poses[ref_List_globe] );
		    transform_edge->setId(pre_edge_count);
                    transform_edge->setVertex(0, v);
                    transform_edge->setMeasurement( poses[ref_This_globe] );
                    //transform_edge->setInformation(eye2d);
                    transform_edge->setRobustKernel( new g2o::RobustKernelHuber() );
                    pre_optimizer.addEdge(transform_edge);
		    
                    cout << "Edge" << pre_edge_count << ": ";
		    cout << "pos1: " << poses[ref_List_globe].transpose() << " " << "with globe ID " << ref_List_globe << " ";
		    cout << "pos2: " << poses[ref_This_globe].transpose() << " " << "with globe ID " << ref_This_globe << endl;
                    pre_edge_count++;
		    read_check++;
                }
            }
        }
        cout << pre_edge_count << " edges" << " has been added to this graph." << endl;

        // Pre optimization for each new traj. //
        pre_optimizer.setVerbose(true);
        pre_optimizer.initializeOptimization();
        cout << "Running pre-optimization..." << endl;
	if(pre_optimizer.verifyInformationMatrices()){
	  pre_optimizer.optimize(30);
	}

        // Do the transformation for current traj and before. //
        cout << "Do the transformation..." << endl;

        Vector4d transform = v->estimate();
	cout << "transform quatenion is " << transform << endl << endl;
        double scale = transform[0];
        double theta = transform[1];
        double tx = transform[2];
        double ty = transform[3];

        int offset = pose_sum[i] - traj_size[i];
        for(int j = 0; j < traj_size[i]; j++){

            Vector2d temp_pose = poses[j + offset];
            double x = scale * ( cos(theta) * temp_pose[0] - sin(theta) * temp_pose[1] ) + tx;
            double y = scale * ( sin(theta) * temp_pose[0] + cos(theta) * temp_pose[1] ) + ty;
            poses[j + offset] = Vector2d(x, y);
        }

        pre_optimizer.clear();
    }
    

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
    cout << "Loop points Optimization... " << endl;
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
        // Add edge. //
    vector<vector<int>>::iterator globe_iter = globeIndexList.begin() - 1;
    while(globe_iter++ != globeIndexList.end()){

        int list_size = (*globe_iter).size();
        vector<int> globeList_comb;
        globeList_comb = combine2(list_size, list_size, globeList_comb);
        if(globeList_comb.size() != 0){
	    cerr << globeList_comb.size() << " combinations." << endl;
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

    delete traj_size, pose_sum, gap_ptr, loop_of_traj;

    return 0;
}