#ifndef _NODE_H_
#define _NODE_H_

#include <iostream>
#include <ros/ros.h>
#include <ros/console.h>
#include <Eigen/Eigen>
#include "backward.hpp"

#define inf 1>>20
struct GridNode;
typedef GridNode* GridNodePtr; // new name for GrideNode type of pointer 


//用结构体变量GridNode 表示，存储了节点的坐标、g(n)、f(n)值、父节点指针等信息。
struct GridNode
{     
    int id;        // 1--> open set, -1 --> closed set
    Eigen::Vector3d coord; // 3×1 vector of type double.
    Eigen::Vector3i dir;   // direction of expanding, 3×1 vector of type integer
    Eigen::Vector3i index;  // 3×1 vector of type integer
	
    double gScore, fScore;  // g(n), f(n)
    GridNodePtr cameFrom; // parent node address
    // declare an empty multimap --> std::multimap<std::string, std::string>mymultimap;
    // declare a multimap iterator
    std::multimap<double, GridNodePtr>::iterator nodeMapIt;  //

    GridNode(Eigen::Vector3i _index, Eigen::Vector3d _coord){  
		id = 0;
		index = _index; // index of x,y,z corrdinate in map?
		coord = _coord;
		dir   = Eigen::Vector3i::Zero(); //returns the (0, 0, 0) vector.

		gScore = inf;
		fScore = inf;
		cameFrom = NULL;  // parent node
    }

    GridNode(){}; //constructor
    ~GridNode(){}; // destructor
};


#endif
