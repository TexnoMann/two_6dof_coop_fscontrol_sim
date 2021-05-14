#include <kdl/chain.hpp>
#include <kdl/chaindynparam.hpp>
#include <kdl/jntspaceinertiamatrix.hpp>
#include <kdl/jntarray.hpp>
#include <kdl/chainjnttojacsolver.hpp>
#include <kdl/chainjnttojacdotsolver.hpp>
#include <kdl/chainiksolverpos_lma.hpp>
#include <kdl/chainiksolvervel_pinv.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/frames.hpp>
#include <kdl_parser/kdl_parser.hpp>

#include "Eigen/Eigen"
#include "Eigen/Core"
#include "Eigen/Geometry"
#include <string>

#include "ros/ros.h"
#include "std_msgs/String.h"
#include "std_msgs/Float64.h"
#include "sensor_msgs/JointState.h"
#include "math.h"


#include <iostream>
#include <boost/array.hpp>

#include <boost/numeric/odeint.hpp>
#include <functional>

namespace pl = std::placeholders;

using namespace std;
using namespace boost::numeric::odeint;

typedef std::vector<double> state_type;
typedef std::vector<state_type> joint_state_data;

#define dt 0.01
joint_state_data jsd;

class RobotState{
  Eigen::VectorXd q;
  Eigen::VectorXd dq;
  Eigen::VectorXd ddq;
  Eigen::VectorXd pa;
  Eigen::VectorXd vw;
  Eigen::VectorXd ae;

  Eigen::MatrixXd J;
  Eigen::MatrixXd Ja;

  Eigen::MatrixXd M;
  Eigen::VectorXd C;
  Eigen::VectorXd g;

  Eigen::MatrixXd Me;
  Eigen::VectorXd Ce;
  Eigen::VectorXd ge;


  Eigen::VectorXd f;
  Eigen::VectorXd tau;
  Eigen::VectorXd ddpc;
  Eigen::VectorXd dpc;
  Eigen::VectorXd pc;

  Eigen::VectorXd SM;


  Eigen::VectorXd q_ref;


  std::shared_ptr<KDL::ChainDynParam> dyn_model;
  KDL::Vector grav;
  KDL::ChainFkSolverPos_recursive fksolverpose;
  KDL::ChainJntToJacSolver jacsolver;
  KDL::ChainJntToJacDotSolver djacsolver;

  KDL::Chain chaint;

public:
  RobotState(KDL::Chain & chain):
    grav(0, 0, -9.82), chaint(chain), fksolverpose(chain), jacsolver(chain), djacsolver(chain)
    {
// State
    q = Eigen::VectorXd::Zero(6);
    dq = Eigen::VectorXd::Zero(6);
    ddq = Eigen::VectorXd::Zero(6);
    pa = Eigen::VectorXd::Zero(6);
    vw = Eigen::VectorXd::Zero(6);
    ae = Eigen::VectorXd::Zero(6);
    q_ref = Eigen::VectorXd::Zero(6);
// Dynamic
    M = Eigen::MatrixXd::Zero(6, 6);
    C = Eigen::VectorXd::Zero(6);
    g = Eigen::VectorXd::Zero(6);
    Me = Eigen::MatrixXd::Zero(6, 6);
    Ce = Eigen::VectorXd::Zero(6);
    ge = Eigen::VectorXd::Zero(6);
// Kinematic
    J = Eigen::MatrixXd::Zero(6, 6);
    Ja = Eigen::MatrixXd::Zero(6, 6);
// Control
    f = Eigen::VectorXd::Zero(6);
    tau =  Eigen::VectorXd::Zero(6);
    ddpc = Eigen::VectorXd::Zero(6);
    dpc =  Eigen::VectorXd::Zero(6);
    pc =  Eigen::VectorXd::Zero(6);

    dyn_model = std::make_shared<KDL::ChainDynParam>(chaint, grav);


    // Control references
    q_ref<< M_PI, -M_PI/2, M_PI/2, 0, -M_PI/2, 0;
  }

  void set_q(state_type & qs){
    for(int i = 0;i<6; i++){
      q[i] = qs[i];
      dq[i] = qs[i+6];
    }
  }

  void operator()( const state_type &qs , state_type &dqsdt , const double t)
  {
    for(int i = 0;i<6; i++){
      q[i] = qs[i];
      dq[i] = qs[i+6];
      pc[i] = qs[i+12];
      dpc[i] = qs[i+18];
    }
    // Convert to kdl
    KDL::JntArray q_kdl;
    q_kdl.data = q;
    KDL::JntArray dq_kdl;
    dq_kdl.data = dq;


    KDL::JntSpaceInertiaMatrix M_kdl(6);
    KDL::JntArray C_kdl(6);
    KDL::JntArray g_kdl(6);
    KDL::Jacobian jac(6);

    dyn_model->JntToMass(dq_kdl, M_kdl);
    dyn_model->JntToCoriolis(q_kdl, dq_kdl, C_kdl);
    dyn_model->JntToGravity(q_kdl, g_kdl);
    // std::cout<<q_kdl(0)<<"\n";

    // Convert back to the eigen
    M = M_kdl.data;
    C = C_kdl.data;
    g = g_kdl.data;

    // Kinematics
    jacsolver.JntToJac(q_kdl, jac);
    J = jac.data;

    // Calc cartesian matrices
    Me = (J.transpose()).inverse()*M*(J.inverse());
    Ce = (J.transpose()).inverse()*C;
    ge = (J.transpose()).inverse()*g;
    std::cout<<ge<<"\n";
    // Calc control
    control(q_kdl, dq_kdl);

    // std::cout<<M<<"\n";
    ddq = M.inverse() * ( tau - (C + g + f));

    // Calc constrain motion
    // ddpc =

    // Comntaining in one vector
    for(int i = 0;i<6; i++){
      dqsdt[i] = dq[i];
      dqsdt[i+6] = ddq[i];
      dqsdt[i+12] = dpc[i];
      dqsdt[i+18] = ddpc[i];
    }
  }
  void control(KDL::JntArray & q_kdl, KDL::JntArray & dq_kdl){
    Eigen::VectorXd e = q_ref - q;
    Eigen::VectorXd ul = 300*e - 200*dq;
    tau = (M*ul + C + g + f);
  }
};


void write_state( const state_type &qs , const double t )
{
    jsd.push_back(qs);
    std::cout<<"Time: "<<t<<", q: "<< qs[0]<<", "<<qs[1]<<", "<<qs[2]<<", "<<qs[3]<<", "<<qs[4]<<", "<<qs[5]<<"\n";
    std::cout<<jsd.size()<<"\n";
}





int main(int argc, char **argv) {
  ros::init(argc, argv, "main_test");
  ros::NodeHandle nh;

  // ros::Publisher robot_state_publisher = n.advertise<sensor_msgs::JointState>("writer", 1000);
  // ros::Rate loop_rate((int)(1/dt));


// Get parameters for kinematic from setup.yaml file
  std::string urdf;
  std::string base_link;
  std::string tool_link;
  ros::param::get("/master/base_link_name", base_link);
  ros::param::get("/master/tool_link_name", tool_link);
  ros::param::get("ur_urdf_model", urdf);

// Ros links for visualize
  ros::Publisher pub = nh.advertise<sensor_msgs::JointState>("joint_states", 100);
  ros::Rate loop_rate(200);

// Generate kinematic model for orocos_kdl
  KDL::Tree tree;
  if (!kdl_parser::treeFromString(urdf, tree)) {
      printf("Failed to construct kdl tree\n");
      return -1;
  }
  KDL::Chain chain;
  if(!tree.getChain(base_link, tool_link, chain)){
      ROS_ERROR_STREAM("Failed to get KDL chain from tree ");
      return false;
  }
  ROS_INFO("tip_name:  %s",tool_link.c_str());
  ROS_INFO("root_name: %s",base_link.c_str());

  RobotState rs(chain);

  state_type q0 = {0,-M_PI/2,0,-M_PI/2,0,0,  0,0,0,0,0,0,  0,0,0,0,0,0,  0,0,0,0,0,0}; // initial conditions: q, dq, s, ds
  rs.set_q(q0);
  integrate(rs ,q0 , 0.0 , 5.0 , dt, write_state);


  ROS_INFO("data_size:  %d",jsd.size());
  loop_rate.sleep();

  while(ros::ok()){
    for(int i=0; i<jsd.size(); i++){
      sensor_msgs::JointState msg;
      state_type wq = jsd[i];
      msg.header.stamp = ros::Time::now();
      msg.name = {"shoulder_pan_joint", "shoulder_lift_joint", "elbow_joint", "wrist_1_joint", "wrist_2_joint", "wrist_3_joint"};
      msg.position = {wq.begin(), wq.begin() +6};

      pub.publish(msg);
      ros::spinOnce();
      loop_rate.sleep();
    }
    for(int i=0; i<1000; i++){
      sensor_msgs::JointState msg;
      state_type wq = jsd[jsd.size()-1];
      msg.header.stamp = ros::Time::now();
      msg.name = {"shoulder_pan_joint", "shoulder_lift_joint", "elbow_joint", "wrist_1_joint", "wrist_2_joint", "wrist_3_joint"};
      msg.position = {wq.begin(), wq.begin() +6};
      pub.publish(msg);
      ros::spinOnce();
      loop_rate.sleep();
    }
  }
}
