#include "lec5_hw/visualizer.hpp"
#include "lec5_hw/trajectory.hpp"

#include <ros/ros.h>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/PoseStamped.h>

#include <cmath>
#include <iostream>
#include <vector>

struct Config
{
    std::string targetTopic;
    double clickHeight;
    std::vector<double> initialVel;
    std::vector<double> initialAcc;
    std::vector<double> terminalVel;
    std::vector<double> terminalAcc;
    double allocationSpeed;
    double allocationAcc;
    int maxPieceNum;

    Config(const ros::NodeHandle &nh_priv)
    {
        nh_priv.getParam("TargetTopic", targetTopic);
        nh_priv.getParam("ClickHeight", clickHeight);
        nh_priv.getParam("InitialVel", initialVel);
        nh_priv.getParam("InitialAcc", initialAcc);
        nh_priv.getParam("TerminalVel", terminalVel);
        nh_priv.getParam("TerminalAcc", terminalAcc);
        nh_priv.getParam("AllocationSpeed", allocationSpeed);
        nh_priv.getParam("AllocationAcc", allocationAcc);
        nh_priv.getParam("MaxPieceNum", maxPieceNum);
    }
};

double timeTrapzVel(const double dist,
                    const double vel,
                    const double acc)
{ // calculate the duration of the piece
    const double t = vel / acc;
    const double d = 0.5 * acc * t * t;

    if (dist < d + d)
    {
        return 2.0 * sqrt(dist / acc);
    }
    else
    {
        return 2.0 * t + (dist - 2.0 * d) / vel;
    }
}



Eigen::Matrix<double, 1, 12> derivative_constraints(const double T)
{ // first order derivative constraints
    Eigen::Matrix<double, 1, 12> derivative_constraints;
    derivative_constraints.setZero();
    derivative_constraints(0, 0) = 1.0;
    derivative_constraints(0, 1) = pow(T,1);
    derivative_constraints(0, 2) = pow(T,2);
    derivative_constraints(0, 3) = pow(T,3);
    derivative_constraints(0, 4) = pow(T,4);
    derivative_constraints(0, 5) = pow(T,5);
    return derivative_constraints;
}

Eigen::Matrix<double, 1, 12>  continuous_position(const double T)
{ // continuous constraints for position
    Eigen::Matrix<double, 1, 12> continuous_position;
    continuous_position.setZero();
    continuous_position(0, 0) = 1.0;
    continuous_position(0, 1) = pow(T,1);
    continuous_position(0, 2) = pow(T,2);
    continuous_position(0, 3) = pow(T,3);
    continuous_position(0, 4) = pow(T,4);
    continuous_position(0, 5) = pow(T,5);
    continuous_position(0, 6) = -1.0;
    return continuous_position;
}

Eigen::Matrix<double, 1, 12>  continuous_velocity(const double T)
{// continuous constraints for velocity
    Eigen::Matrix<double, 1, 12> continuous_velocity;
    continuous_velocity.setZero();
    continuous_velocity(0, 1) = 1.0;
    continuous_velocity(0, 2) = 2.0 * pow(T,1);
    continuous_velocity(0, 3) = 3.0 * pow(T,2);
    continuous_velocity(0, 4) = 4.0 * pow(T,3);
    continuous_velocity(0, 5) = 5.0 * pow(T,4);
    continuous_velocity(0, 7) = -1.0; // velocity continuous for next segment
    return continuous_velocity;
}

Eigen::Matrix<double, 1, 12> continuous_acceleration(const double T)
{
    Eigen::Matrix<double, 1, 12> continuous_acceleration;
    continuous_acceleration.setZero();
    continuous_acceleration(0, 2) = 2.0;
    continuous_acceleration(0, 3) = 6.0 * pow(T,1);
    continuous_acceleration(0, 4) = 12.0 * pow(T,2);
    continuous_acceleration(0, 5) = 20.0 * pow(T,3);
    continuous_acceleration(0, 8) = -2.0;
    return continuous_acceleration;
}

Eigen::Matrix<double, 1, 12> continuous_jerk(const double T)
{
    Eigen::Matrix<double, 1, 12> continuous_jerk;
    continuous_jerk.setZero();
    continuous_jerk(0, 3) = 6.0;
    continuous_jerk(0, 4) = 24.0 * pow(T,1);
    continuous_jerk(0, 5) = 60.0 * pow(T,2);
    continuous_jerk(0, 9) = -6.0;
    return continuous_jerk;
}

Eigen::Matrix<double, 1, 12> continuous_snap(const double T)
{
    Eigen::Matrix<double, 1, 12> continuous_snap;
    continuous_snap.setZero();
    
    continuous_snap(0, 4) = 24.0;
    continuous_snap(0, 5) = 120.0 * pow(T,1);
    continuous_snap(0, 10) = -24.0;
    return continuous_snap;
}

Eigen::Matrix<double, 3, 6> target_state(const double T)
{
    Eigen::Matrix<double, 3, 6> target_state;
    target_state << 1.0, T, pow(T, 2.0), pow(T, 3.0), pow(T, 4.0), pow(T, 5.0),
                   0.0, 1.0, 2.0 * T, 3.0 * pow(T, 2.0), 4.0 * pow(T, 3.0), 5.0 * pow(T, 4.0),
                   0.0, 0.0, 2.0, 6.0 * T, 12.0 * pow(T, 2.0), 20.0 * pow(T, 3.0);

    return target_state;
}




void minimumJerkTrajGen(
    // Inputs:
    const int pieceNum,
    const Eigen::Vector3d &initialPos,
    const Eigen::Vector3d &initialVel,
    const Eigen::Vector3d &initialAcc,
    const Eigen::Vector3d &terminalPos,
    const Eigen::Vector3d &terminalVel,
    const Eigen::Vector3d &terminalAcc,
    const Eigen::Matrix3Xd &intermediatePositions,
    const Eigen::VectorXd &timeAllocationVector,
    // Outputs:
    Eigen::MatrixX3d &coefficientMatrix) //size of (6*num_piece , 3) 
{
    // coefficientMatrix is a matrix with 6*piece num rows and 3 columes
    // As for a polynomial c0+c1*t+c2*t^2+c3*t^3+c4*t^4+c5*t^5,
    // each 6*3 sub-block of coefficientMatrix is
    // --              --
    // | c0_x c0_y c0_z |
    // | c1_x c1_y c1_z |
    // | c2_x c2_y c2_z |
    // | c3_x c3_y c3_z |
    // | c4_x c4_y c4_z |
    // | c5_x c5_y c5_z |
    // --              --
    // Please computed coefficientMatrix of the minimum-jerk trajectory
    // in this function

    // ------------------------ Put your solution below ------------------------

    int coeffMatRow = coefficientMatrix.rows();
    //int intermediateNum = intermediatePositions.rows();


    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> M;
    M.setZero(coeffMatRow, coeffMatRow);
    Eigen::MatrixX3d b = Eigen::MatrixXd::Zero(coeffMatRow, 3);

    //Eigen::Matrix<double, coeffMatRow, coeffMatRow> M; // Mc = b
    //M.setZero(coeffMatRow, coeffMatRow);
    //Eigen::Matrix<double, coeffMatRow, 3> b; // Mc = b 
    //b.setZero(coeffMatRow, 3); 


    M(0, 0) = 1.0;
    M(1, 1) = 1.0;
    M(2, 2) = 2.0;


    b(0,0) = initialPos(0);
    b(0,1) = initialPos(1);
    b(0,2) = initialPos(2);
    b(1,0) = initialVel(0);
    b(1,1) = initialVel(1);
    b(1,2) = initialVel(2);
    b(2,0) = initialAcc(0);
    b(2,1) = initialAcc(1);
    b(2,2) = initialAcc(2);   
    


    for (int i = 0; i < pieceNum-1; i ++)
    {
        b(3 + 6 * i, 0) = intermediatePositions(0, i);
        b(3 + 6 * i, 1) = intermediatePositions(1, i);
        b(3 + 6 * i, 2) = intermediatePositions(2, i);


        M.block(3 + 6 * i + 0, 6 * i, 1, 12) = derivative_constraints(timeAllocationVector(i));
        M.block(3 + 6 * i + 1, 6 * i, 1, 12) = continuous_position(timeAllocationVector(i));
        M.block(3 + 6 * i + 2, 6 * i, 1, 12) = continuous_velocity(timeAllocationVector(i));
        M.block(3 + 6 * i + 3, 6 * i, 1, 12) = continuous_acceleration(timeAllocationVector(i));
        M.block(3 + 6 * i + 4, 6 * i, 1, 12) = continuous_jerk(timeAllocationVector(i));
        M.block(3 + 6 * i + 5, 6 * i, 1, 12) = continuous_snap(timeAllocationVector(i));


    }

    M.block(coeffMatRow - 3, coeffMatRow - 6, 3, 6) = target_state(timeAllocationVector(pieceNum - 1));

    b(coeffMatRow-3,0) = terminalPos(0);
    b(coeffMatRow-3,1) = terminalPos(1);
    b(coeffMatRow-3,2) = terminalPos(2);
    b(coeffMatRow-2,0) = terminalVel(0);
    b(coeffMatRow-2,1) = terminalVel(1);
    b(coeffMatRow-2,2) = terminalVel(2);
    b(coeffMatRow-1,0) = terminalAcc(0);
    b(coeffMatRow-1,1) = terminalAcc(1);
    b(coeffMatRow-1,2) = terminalAcc(2);  
    

    coefficientMatrix = M.inverse() * b;


    // ------------------------ Put your solution above ------------------------
}

class ClickGen
{
private:
    Config config;

    ros::NodeHandle nh;
    ros::Subscriber targetSub;

    Visualizer visualizer;

    Eigen::Matrix3Xd positions;
    Eigen::VectorXd times;
    int positionNum;
    Trajectory<5> traj;

public:
    ClickGen(const Config &conf,
             ros::NodeHandle &nh_)
        : config(conf),
          nh(nh_),
          visualizer(nh),
          positions(3, config.maxPieceNum + 1),
          times(config.maxPieceNum),
          positionNum(0)
    {
        targetSub = nh.subscribe(config.targetTopic, 1,
                                 &ClickGen::targetCallBack, this,
                                 ros::TransportHints().tcpNoDelay());
    }

    void targetCallBack(const geometry_msgs::PoseStamped::ConstPtr &msg)
    {
        if (positionNum > config.maxPieceNum)
        {
            positionNum = 0;
            traj.clear();
        }

        positions(0, positionNum) = msg->pose.position.x;
        positions(1, positionNum) = msg->pose.position.y;
        positions(2, positionNum) = std::fabs(msg->pose.orientation.z) * config.clickHeight;

        if (positionNum > 0)
        {
            const double dist = (positions.col(positionNum) - positions.col(positionNum - 1)).norm();
            times(positionNum - 1) = timeTrapzVel(dist, config.allocationSpeed, config.allocationAcc);
        }

        ++positionNum;

        if (positionNum > 1)
        {
            const int pieceNum = positionNum - 1;
            const Eigen::Vector3d initialPos = positions.col(0);
            const Eigen::Vector3d initialVel(config.initialVel[0], config.initialVel[1], config.initialVel[2]);
            const Eigen::Vector3d initialAcc(config.initialAcc[0], config.initialAcc[1], config.initialAcc[2]);
            const Eigen::Vector3d terminalPos = positions.col(pieceNum);
            const Eigen::Vector3d terminalVel(config.terminalVel[0], config.terminalVel[1], config.terminalVel[2]);
            const Eigen::Vector3d terminalAcc(config.terminalAcc[0], config.terminalAcc[1], config.terminalAcc[2]);
            const Eigen::Matrix3Xd intermediatePositions = positions.middleCols(1, pieceNum - 1);
            const Eigen::VectorXd timeAllocationVector = times.head(pieceNum);

            Eigen::MatrixX3d coefficientMatrix = Eigen::MatrixXd::Zero(6 * pieceNum, 3);

            minimumJerkTrajGen(pieceNum,
                               initialPos, initialVel, initialAcc,
                               terminalPos, terminalVel, terminalAcc,
                               intermediatePositions,
                               timeAllocationVector,
                               coefficientMatrix);

            traj.clear();
            traj.reserve(pieceNum);
            for (int i = 0; i < pieceNum; i++)
            {
                traj.emplace_back(timeAllocationVector(i),
                                  coefficientMatrix.block<6, 3>(6 * i, 0).transpose().rowwise().reverse());
            }
        }

        visualizer.visualize(traj, positions.leftCols(positionNum));

        return;
    }
};

int main(int argc, char **argv)
{
    ros::init(argc, argv, "click_gen_node");
    ros::NodeHandle nh_;
    ClickGen clickGen(Config(ros::NodeHandle("~")), nh_);
    ros::spin();
    return 0;
}
