#include "Astar_searcher.h"

using namespace std;
using namespace Eigen;

// initialize grid map
void AstarPathFinder::initGridMap(double _resolution, Vector3d global_xyz_l, Vector3d global_xyz_u, int max_x_id, int max_y_id, int max_z_id)
{   
    // lower bound of map
    gl_xl = global_xyz_l(0);
    gl_yl = global_xyz_l(1);
    gl_zl = global_xyz_l(2);

    //upper bound of map
    gl_xu = global_xyz_u(0);
    gl_yu = global_xyz_u(1);
    gl_zu = global_xyz_u(2);
    
    //grid map coordinate lower bound
    GLX_SIZE = max_x_id;
    GLY_SIZE = max_y_id;
    GLZ_SIZE = max_z_id;
    GLYZ_SIZE  = GLY_SIZE * GLZ_SIZE;
    GLXYZ_SIZE = GLX_SIZE * GLYZ_SIZE;
    
    // 
    resolution = _resolution;
    inv_resolution = 1.0 / _resolution;    

    data = new uint8_t[GLXYZ_SIZE];
    memset(data, 0, GLXYZ_SIZE * sizeof(uint8_t));
    
    // allocate grid map
    GridNodeMap = new GridNodePtr ** [GLX_SIZE]; //pointer of pointer of pointer
    for(int i = 0; i < GLX_SIZE; i++){
        GridNodeMap[i] = new GridNodePtr * [GLY_SIZE];
        for(int j = 0; j < GLY_SIZE; j++){
            GridNodeMap[i][j] = new GridNodePtr [GLZ_SIZE];
            for( int k = 0; k < GLZ_SIZE;k++){
                Vector3i tmpIdx(i,j,k);
                Vector3d pos = gridIndex2coord(tmpIdx);
                GridNodeMap[i][j][k] = new GridNode(tmpIdx, pos);
            }
        }
    }
}

void AstarPathFinder::resetGrid(GridNodePtr ptr)
{
    ptr->id = 0;
    ptr->cameFrom = NULL;
    ptr->gScore = inf;
    ptr->fScore = inf;
}

void AstarPathFinder::resetUsedGrids()
{   
    for(int i=0; i < GLX_SIZE ; i++)
        for(int j=0; j < GLY_SIZE ; j++)
            for(int k=0; k < GLZ_SIZE ; k++)
                resetGrid(GridNodeMap[i][j][k]);
}

void AstarPathFinder::setObs(const double coord_x, const double coord_y, const double coord_z)
{
    if( coord_x < gl_xl  || coord_y < gl_yl  || coord_z <  gl_zl || 
        coord_x >= gl_xu || coord_y >= gl_yu || coord_z >= gl_zu )
        return;

    int idx_x = static_cast<int>( (coord_x - gl_xl) * inv_resolution);
    int idx_y = static_cast<int>( (coord_y - gl_yl) * inv_resolution);
    int idx_z = static_cast<int>( (coord_z - gl_zl) * inv_resolution);      

    data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] = 1;
}

vector<Vector3d> AstarPathFinder::getVisitedNodes()
{   
    vector<Vector3d> visited_nodes;
    for(int i = 0; i < GLX_SIZE; i++)
        for(int j = 0; j < GLY_SIZE; j++)
            for(int k = 0; k < GLZ_SIZE; k++){   
                if(GridNodeMap[i][j][k]->id != 0) // visualize all nodes in open and close list
                //if(GridNodeMap[i][j][k]->id == -1)  // visualize nodes in close list only
                    visited_nodes.push_back(GridNodeMap[i][j][k]->coord);
            }

    ROS_WARN("visited_nodes size : %d", visited_nodes.size());
    return visited_nodes;
}

Vector3d AstarPathFinder::gridIndex2coord(const Vector3i & index) 
{
    Vector3d pt;

    pt(0) = ((double)index(0) + 0.5) * resolution + gl_xl;
    pt(1) = ((double)index(1) + 0.5) * resolution + gl_yl;
    pt(2) = ((double)index(2) + 0.5) * resolution + gl_zl;

    return pt;
}

Vector3i AstarPathFinder::coord2gridIndex(const Vector3d & pt) 
{
    Vector3i idx;
    idx <<  min( max( int( (pt(0) - gl_xl) * inv_resolution), 0), GLX_SIZE - 1),
            min( max( int( (pt(1) - gl_yl) * inv_resolution), 0), GLY_SIZE - 1),
            min( max( int( (pt(2) - gl_zl) * inv_resolution), 0), GLZ_SIZE - 1);                  
  
    return idx;
}

Eigen::Vector3d AstarPathFinder::coordRounding(const Eigen::Vector3d & coord)
{
    return gridIndex2coord(coord2gridIndex(coord));
}

inline bool AstarPathFinder::isOccupied(const Eigen::Vector3i & index) const
{
    return isOccupied(index(0), index(1), index(2));
}

inline bool AstarPathFinder::isFree(const Eigen::Vector3i & index) const
{
    return isFree(index(0), index(1), index(2));
}

inline bool AstarPathFinder::isOccupied(const int & idx_x, const int & idx_y, const int & idx_z) const 
{
    return  (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE && idx_z >= 0 && idx_z < GLZ_SIZE && 
            (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] == 1));
}

inline bool AstarPathFinder::isFree(const int & idx_x, const int & idx_y, const int & idx_z) const 
{
    return (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE && idx_z >= 0 && idx_z < GLZ_SIZE && 
           (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] < 1));
}

inline void AstarPathFinder::AstarGetSucc(GridNodePtr currentPtr, vector<GridNodePtr> & neighborPtrSets, vector<double> & edgeCostSets)
{   
    neighborPtrSets.clear();
    edgeCostSets.clear();

    //STEP 4: finish AstarPathFinder::AstarGetSucc yourself 
    Eigen::Vector3i current_node_index = currentPtr->index;
    int current_index_x = current_node_index[0];
    int current_index_y = current_node_index[1];
    int current_index_z = current_node_index[2];

    int node_x, node_y, node_z;
    GridNodePtr next_ptr = NULL;

    // 1-step ahead expanding for x, y, z dimensions across 3 headings respectively(8 connections?)
    for(int i=-1;i <=1 ;++i)
    {
        for(int j=-1;j<=1;++j)
        {
            for(int k = -1;k<=1;++k)
            {
                if((i==0) && (j==0) && (k==0)){continue;} // stay steel

                node_x = current_index_x + i; // move in x direction by i step
                node_y = current_index_y + j; // move in y direction by j step
                node_z = current_index_z + k; // move in z direction by k step

                if ((node_x<0) || (node_y<0) || (node_z<0) || (node_x > GLX_SIZE-1) || (node_y > GLY_SIZE-1) || (node_z > GLZ_SIZE-1))
                { // move outside of the map
                  continue;
                }

                if (isOccupied(node_x, node_y, node_z)){continue;} //obstacles

                next_ptr = GridNodeMap[node_x][node_y][node_z]; // store pointer into next node pointer

                double cnm = getHeu(currentPtr, next_ptr); // edge cost between two nodes
                
                neighborPtrSets.push_back(next_ptr);
                edgeCostSets.push_back(cnm);
            }
        }
    }
}

double AstarPathFinder::getHeu(GridNodePtr node1, GridNodePtr node2)
{

    //STEP 1: finish the AstarPathFinder::getHeu , which is the heuristic function

    //Vector3d difference = (*node1.coord - *node2.coord);
    Vector3d difference = (node1->coord - node2->coord);

    // tie breaker
    //double h_multiplier = 1;
    double h_multiplier = 1 + 0.001;
    /*
    #define tie_break 1
    #if tie_break
    {
        double h_multiplier = 1 + 0.001;
    }
    else
    {
        double h_multiplier = 1;
    }
    #endif 
    */

    // euclidean_distance
    #define Euclidean 0 
    #if Euclidean
    {
        double euclidean_distance = difference.norm(); 
        return euclidean_distance * h_multiplier;
    }
    #endif

    //manhattan_distance
    #define Manhatten 0
    #if Manhatten
    {
        double manhattan_distance = difference.cwiseAbs().sum();
        return manhattan_distance * h_multiplier;
    }
    #endif
    
    //dijkstra_distance
    #define Dijkstra 0
    #if Dijkstra
    {
        double dijkstra_distance = 0.0;
        return dijkstra_distance * h_multiplier;
    }
    #endif

    // diagonal_distance
    #define Diagonal 1
    #if Diagonal
    {
        double D1 = 1;
        double D2 = sqrt(2);
        double D3 = sqrt(3);
        
        double dx = abs(node1->coord(0) - node2->coord(0));
        double dy = abs(node1->coord(1) - node2->coord(1));
        double dz = abs(node1->coord(2) - node2->coord(2));
        double dmin = min(dx, dy);
        dmin = min(dmin, dz);
        double dmax = max(dx, dy);
        dmax = max(dmax, dz);
        double dmid = dx + dy + dz - dmin - dmax;
        double diagonal_distance = (D3-D2)*dmin + (D2-D1)*dmid + D1*dmax;
        return diagonal_distance * h_multiplier;
    }
    #endif
}

void AstarPathFinder::AstarGraphSearch(Vector3d start_pt, Vector3d end_pt)
{   
    ros::Time time_1 = ros::Time::now();    

    //index of start_point and end_point
    Vector3i start_idx = coord2gridIndex(start_pt);
    Vector3i end_idx   = coord2gridIndex(end_pt);
    goalIdx = end_idx; //target point index

    //position of start_point and end_point
    start_pt = gridIndex2coord(start_idx);
    end_pt   = gridIndex2coord(end_idx);

    //Initialize the pointers of struct GridNode which represent start node and goal node
    GridNodePtr startPtr = new GridNode(start_idx, start_pt);
    GridNodePtr endPtr   = new GridNode(end_idx,   end_pt);

    //openSet is the open_list implemented through multimap in STL library
    openSet.clear();
    // currentPtr represents the node with lowest f(n) in the open_list
    GridNodePtr currentPtr  = NULL;
    GridNodePtr neighborPtr = NULL;

    //put start node in open set
    startPtr -> gScore = 0;
    startPtr -> fScore = getHeu(startPtr,endPtr);   
    //STEP 1: finish the AstarPathFinder::getHeu , which is the heuristic function
    startPtr -> id = 1; // unexpand=0， closed=-1， expanded=1
    startPtr -> coord = start_pt;  // coordinate
    startPtr->nodeMapIt = openSet.insert(make_pair(startPtr -> fScore, startPtr));  //openSet:<fscore, pointer to GridNode type>

    //STEP 2 :  some else preparatory works which should be done before while loop
    vector<GridNodePtr> neighborPtrSets;  // vector stores neighbour node pointers
    vector<double> edgeCostSets;  //cost between nodes, c(n)

    // this is the main loop
    while ( !openSet.empty() ){

        //step 3: Remove the node with lowest cost function from open set to closed set
        //find node with leaset fScore
        currentPtr = openSet.begin()->second;  //multimap<T1, T2, less> ,即按key的升序排列
        currentPtr->id = -1; // set current node to closed set
        
        // remove node with least fscore
        openSet.erase(openSet.begin());
        //Eigen::Vector3i current_idx = currentPtr->index;

        //STEP 4:  get neighbour nodes of the current node
        AstarGetSucc(currentPtr, neighborPtrSets, edgeCostSets);  
       
        //STEP 5:  For all unexpanded neigbors "m" of node "n"
        for(int i = 0; i < (int)neighborPtrSets.size(); i++)
        {
            /*
            IMPORTANT NOTE!!!
            neighborPtrSets[i]->id = -1 : expanded, equal to this node is in close set
            neighborPtrSets[i]->id = 1 : unexpanded, equal to this node is in open set        
            */
            neighborPtr = neighborPtrSets[i];
            double gn = neighborPtr->gScore + edgeCostSets[i]; //update gn 
            double fn = gn + getHeu(neighborPtr, endPtr); // update fn

            if(neighborPtr -> id == 0)
            { 
                //discover a new node, which is not in the closed set and open set
                //STEP 6:  As for a new node, do what you need do ,and then put neighbor in open set and record it
                neighborPtr -> gScore = gn;
                neighborPtr -> fScore = fn;
                neighborPtr -> cameFrom = currentPtr; // pointer to parent node 
                neighborPtr->nodeMapIt = openSet.insert(make_pair(neighborPtr -> fScore, neighborPtr));
                neighborPtr->id = 1;

                // if the neighbour node is the goal 
                if(neighborPtr->index == goalIdx)
                {
                    ros::Time time_2 = ros::Time::now();
                    terminatePtr = neighborPtr;
                    ROS_WARN("[A*]{sucess}  Time in A*  is %f ms, path cost if %f m", (time_2 - time_1).toSec() * 1000.0, neighborPtr->gScore * resolution);            
                    return;
                }
                else
                {
                    neighborPtr->id = 1;
                    continue;
                }
            }
            else if(neighborPtr->id == 1)
            { //this node is in open set and need to judge if it needs to update, the "0" should be deleted when you are coding
                //STEP 7:  As for a node in open set, update it , maintain the openset ,and then put neighbor in open set and record it      
                if (neighborPtr->gScore > gn)
                {// update node parameter in openset
                    neighborPtr->gScore = gn;
                    neighborPtr->fScore = fn;
                    neighborPtr->cameFrom = currentPtr;
                    openSet.erase(neighborPtr->nodeMapIt);
                    neighborPtr->nodeMapIt = openSet.insert(make_pair(neighborPtr->fScore,neighborPtr));
                }
            }
            else
            {//this node is in closed set
                continue;
            }
        }      
    }
    //if search fails
    ros::Time time_2 = ros::Time::now();
    if((time_2 - time_1).toSec() > 0.1)
        ROS_WARN("Time consume in Astar path finding is %f", (time_2 - time_1).toSec() );
}


vector<Vector3d> AstarPathFinder::getPath() 
{   
    vector<Vector3d> path;
    vector<GridNodePtr> gridPath;
    //STEP 8:  trace back from the curretnt nodePtr to get all nodes along the path
    GridNodePtr temp_node_ptr = terminatePtr;
    while( temp_node_ptr->cameFrom != NULL)
    {
        gridPath.push_back(temp_node_ptr);
        temp_node_ptr = temp_node_ptr->cameFrom;
    }

    for (auto ptr: gridPath)
        path.push_back(ptr->coord);
        
    reverse(path.begin(),path.end());
    return path;
}
