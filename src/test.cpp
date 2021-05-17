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

#include <Eigen/Eigen>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <string>

#include <math.h>


#include <iostream>
#include <boost/array.hpp>

#include <boost/numeric/odeint.hpp>
#include <functional>

#include <fstream>

#define SIGN(a) ((a >= 0) - (a <= 0))

namespace pl = std::placeholders;

using namespace std;
// using namespace boost::numeric::odeint;

typedef std::vector<double> state_type;
typedef std::vector<state_type> joint_state_data;

#define dt 0.001
joint_state_data jsd;

class RobotState{


	KDL::JntArrayVel jnt;
	Eigen::VectorXd ddq;
	Eigen::VectorXd pa;
	Eigen::VectorXd vw;
	Eigen::VectorXd ae;

	double phi, theta, psi; // Euler angles
	double dphi, dtheta, dpsi; // Euler angles
	KDL::Frame fk_H;
	KDL::Jacobian J;
	KDL::Jacobian J_dot;
	Eigen::MatrixXd Ja;
	Eigen::MatrixXd Ja_dot;
	Eigen::MatrixXd inv_Ta;
	Eigen::MatrixXd Ta_dot;

	KDL::JntSpaceInertiaMatrix M;
	KDL::JntArray C;
	KDL::JntArray g;

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
	Eigen::VectorXd p_ref;
	Eigen::VectorXd p_r;

	Eigen::MatrixXd Kpl;
	Eigen::MatrixXd Kdl;
	Eigen::MatrixXd Kpf;
	Eigen::MatrixXd Kif;
	Eigen::MatrixXd T;
	Eigen::MatrixXd Y;
	double z_I, z_I_dot;
	double Fz_ref;


	std::shared_ptr<KDL::ChainDynParam> dyn_model;
	KDL::Vector gravity;
	KDL::ChainFkSolverPos_recursive fksolverpose;
	KDL::ChainJntToJacSolver jacsolver;
	KDL::ChainJntToJacDotSolver djacsolver;

	KDL::Chain chaint;

	std::shared_ptr<std::ofstream> myfile;

public:
	RobotState(KDL::Chain & chain, std::shared_ptr<std::ofstream> file)
	: gravity(0, 0, -9.82), chaint(chain), fksolverpose(chain), jacsolver(chain), djacsolver(chain)
	{

		// State
		jnt		= KDL::JntArrayVel(6);
		ddq 	= Eigen::VectorXd::Zero(6);
		pa		= Eigen::VectorXd::Zero(6);
		vw		= Eigen::VectorXd::Zero(6);
		ae		= Eigen::VectorXd::Zero(6);
		q_ref 	= Eigen::VectorXd::Zero(6);
		p_ref 	= Eigen::VectorXd::Zero(6);
		p_r		= Eigen::VectorXd::Zero(6);

		// std::cout << std::endl;
		// std::cout << q.data.transpose() << std::endl;

		// Dynamic
		M		= KDL::JntSpaceInertiaMatrix(6);
		C		= KDL::JntArray(6);
		g		= KDL::JntArray(6);
		Me		= Eigen::MatrixXd::Zero(6, 6);
		Ce		= Eigen::VectorXd::Zero(6);
		ge		= Eigen::VectorXd::Zero(6);

		// Kinematic
		phi 	= 0;
		theta	= 0;
		psi		= 0;
		dphi 	= 0;
		dtheta	= 0;
		dpsi	= 0;
		fk_H	= KDL::Frame::Identity();
		J		= KDL::Jacobian(6);
		J_dot	= KDL::Jacobian(6);
		Ja		= Eigen::MatrixXd::Zero(6, 6);
		Ja_dot	= Eigen::MatrixXd::Zero(6, 6);
		inv_Ta	= Eigen::MatrixXd::Identity(6, 6);
		Ta_dot	= Eigen::MatrixXd::Zero(6, 6);

		// Control
		f		= Eigen::VectorXd::Zero(6);
		tau		= Eigen::VectorXd::Zero(6);
		ddpc	= Eigen::VectorXd::Zero(6);
		dpc		= Eigen::VectorXd::Zero(6);
		pc		= Eigen::VectorXd::Zero(6);
		Kpl		= Eigen::MatrixXd::Identity(6, 6);
		Kdl		= Eigen::MatrixXd::Identity(6, 6);
		Kpf		= Eigen::MatrixXd::Identity(6, 6);
		Kif		= Eigen::MatrixXd::Identity(6, 6);
		T		= Eigen::MatrixXd::Identity(6, 6);
		Y		= Eigen::MatrixXd::Identity(6, 6);
		z_I		= 0;
		z_I_dot = 0;

		Kpl(0, 0) = 100; Kpl(1, 1) = 100; Kpl(2, 2) = 100;
		Kpl(3, 3) = 100; Kpl(4, 4) = 100; Kpl(5, 5) = 100;

		Kdl(0, 0) = 50; Kdl(1, 1) = 50; Kdl(2, 2) = 50;
		Kdl(3, 3) = 50; Kdl(4, 4) = 50; Kdl(5, 5) = 50;

		Kpf(0, 0) = 0.1; Kpf(1, 1) = 0.1; Kpf(2, 2) = 0.1;
		Kpf(3, 3) = 0.1; Kpf(4, 4) = 0.1; Kpf(5, 5) = 0.1;

		Kif(0, 0) = 0.2; Kif(1, 1) = 0.2; Kif(2, 2) = 0.2;
		Kif(3, 3) = 0.2; Kif(4, 4) = 0.2; Kif(5, 5) = 0.2;

		T(0, 0) = 1; T(1, 1) = 1; T(2, 2) = 0; // Disable z
		T(3, 3) = 1; T(4, 4) = 1; T(5, 5) = 1;

		Y(0, 0) = 0; Y(1, 1) = 0; Y(2, 2) = 1; // Enable z
		Y(3, 3) = 0; Y(4, 4) = 0; Y(5, 5) = 0;

		dyn_model = std::make_shared<KDL::ChainDynParam>(chaint, gravity);
		myfile = file;

		// Control references
		q_ref << M_PI, -M_PI/2, M_PI/2, 0, -M_PI/2, 0;
		p_ref << 0.372, 0.134, 0.1, M_PI/2, 2.0944, M_PI;
		Fz_ref = 20;
		std::cout << p_ref.transpose() << std::endl;
		std::cout << "---" << std::endl;
		std::cout << T << std::endl;
		std::cout << "---" << std::endl;
		std::cout << Y << std::endl;
		std::cout << "---" << std::endl;
		*myfile << "t,ex,ey,ez,ephi,etheta,epsi,fz,fz_dot,z" << std::endl;
	}

	~RobotState()
	{}

	void set_q(state_type & qs){
		for(int i = 0; i < 6; i++){
			jnt.q(i) = qs[i];
			jnt.qdot(i) = qs[i+6];
		}
	}

	// Euler ZYZ calculus
	// std::cout << phi << ", " << theta << ", " << psi << std::endl;
	// double thetat = atan2(sqrt(1 - fk_H.M(2, 2)*fk_H.M(2, 2)), fk_H.M(2, 2));
	// double phit = abs(thetat) > 1e-5 ? atan2(fk_H.M(2, 1), -fk_H.M(2, 0)) : -100;
	// double psit = abs(thetat) > 1e-5 ? atan2(fk_H.M(1, 2), fk_H.M(0, 2)) : -100;
	// std::cout << phit << ", " << thetat << ", " << psit << std::endl;

	void operator()( const state_type &qs , state_type &dqsdt , const double t)
	{
		for(int i = 0; i < 6; i++) {
			jnt.q(i) = qs[i];
			jnt.qdot(i) = qs[i+6];
			// pc[i] = qs[i+12];
			// dpc[i] = qs[i+18];
		}
		z_I = qs[24];

		// Calc Dynamic and Kinematics matrices using KDL
		dyn_model->JntToMass(jnt.q, M);
		dyn_model->JntToCoriolis(jnt.q, jnt.qdot, C);
		dyn_model->JntToGravity(jnt.q, g);
		jacsolver.JntToJac(jnt.q, J);
		djacsolver.JntToJacDot(jnt, J_dot);
		fksolverpose.JntToCart(jnt.q, fk_H);
		fk_H.M.GetEulerZYZ(psi, theta, phi);

		calcAnalitycalJacobian();

		pc[0] = fk_H.p.data[0];
		pc[1] = fk_H.p.data[1];
		pc[2] = fk_H.p.data[2];
		pc[3] = phi;
		pc[4] = theta;
		pc[5] = psi;

		calcConstraints();

		// Calc control
		control(jnt.q, jnt.qdot, t);

		ddq = M.data.inverse() * (tau - C.data - g.data - J.data.transpose()*f);

		// Comntaining in one vector
		for(int i = 0;i<6; i++){
			dqsdt[i] = jnt.qdot(i);
			dqsdt[i+6] = ddq[i];
			dqsdt[i+12] = dpc[i];
			dqsdt[i+18] = ddpc[i];
		}
		dqsdt[24] = z_I_dot;
	}

	void calculateCertasianDynMatrices()
	{
		// Calc cartesian matrices
		Me = (J.data.transpose()).inverse()*M.data*(J.data.inverse());
		Ce = (J.data.transpose()).inverse()*C.data;
		ge = (J.data.transpose()).inverse()*g.data;
		// std::cout<<ge<<"\n";
	}

	// Non iverted
	// Ta(3, 3) = cos(psi)*sin(theta);
	// Ta(4, 3) = sin(psi)*sin(theta);
	// Ta(5, 3) =          cos(theta);

	// Ta(3, 4) = -sin(psi);
	// Ta(4, 4) =  cos(psi);
	// Ta(5, 4) =         0;

	// Ta(3, 5) = 0;
	// Ta(4, 5) = 0;
	// Ta(5, 5) = 1;

	// Ja = Ta.inverse()*J.data;

	void calcAnalitycalJacobian()
	{
		// Inverted Jacobian
		double inv_det = 1/sin(theta);
		inv_Ta(3, 3) = inv_det*cos(psi);
		inv_Ta(4, 3) = -sin(psi);
		inv_Ta(5, 3) = -inv_det*cos(psi)*cos(theta);

		inv_Ta(3, 4) = inv_det*sin(psi);
		inv_Ta(4, 4) = cos(psi);
		inv_Ta(5, 4) = -inv_det*sin(psi)*cos(theta);

		inv_Ta(3, 5) = 0;
		inv_Ta(4, 5) = 0;
		inv_Ta(5, 5) = 1;

		Ja = inv_Ta*J.data;

		// Jacobian derivative
		Eigen::VectorXd p_dot(6);
		p_dot	= Ja*jnt.qdot.data;
		dphi 	= p_dot[3];
		dtheta 	= p_dot[4];
		dpsi 	= p_dot[5];

		Ta_dot(3, 3) = cos(psi)*cos(theta)*dtheta - sin(psi)*sin(theta)*dpsi;
		Ta_dot(4, 3) = cos(theta)*sin(psi)*dtheta + cos(psi)*sin(theta)*dpsi;
		Ta_dot(5, 3) = -sin(theta)*dtheta;

		Ta_dot(3, 4) = -cos(psi)*dpsi;
		Ta_dot(4, 4) = -sin(psi)*dpsi;
		Ta_dot(5, 4) = 0;

		Ta_dot(3, 5) = 0;
		Ta_dot(4, 5) = 0;
		Ta_dot(5, 5) = 0;

		Ja_dot = inv_Ta*(J_dot.data - Ta_dot*Ja);
	}

	void calcConstraints()
	{
		f[2] = (0 - pc[2]) > 0 ? 1e+4*(0 - pc[2]) : 0;
	}

	void control(KDL::JntArray & q_kdl, KDL::JntArray & dq_kdl, const double t)
	{

		// TODO: move it to the heaer
		Eigen::VectorXd e 		= Eigen::VectorXd::Zero(6);
		Eigen::VectorXd dF 		= Eigen::VectorXd::Zero(6);
		Eigen::VectorXd dF_I 	= Eigen::VectorXd::Zero(6);
		Eigen::VectorXd ul 		= Eigen::VectorXd::Zero(6);
		Eigen::VectorXd uf 		= Eigen::VectorXd::Zero(6);

		// Motion and Force control
		z_I_dot = -Kif(2, 2)*(Fz_ref - f[2]);
		std::cout << z_I_dot << std::endl;
		dF[2] = -(Fz_ref - f[2]);
		// dF[2] = -(Fz_ref);
		dF_I[2] = z_I;

		p_r = p_ref + Kpl.inverse()*(Kpf*dF + dF_I);
		e = p_r - pc;
		e[3] = abs(e[3]) > 1.2*M_PI ? e[3] - SIGN(e[3])*2*M_PI : e[3];
		e[4] = abs(e[4]) > 1.2*M_PI ? e[4] - SIGN(e[4])*2*M_PI : e[4];
		e[5] = abs(e[5]) > 1.2*M_PI ? e[5] - SIGN(e[5])*2*M_PI : e[5];
		ul = Kpl*e - Kdl*Ja*jnt.qdot.data;
		uf = f;

		tau = M.data*Ja.inverse()*(ul - Ja_dot*jnt.qdot.data) + J.data.transpose()*Y*uf + C.data + g.data;

		// tau = M.data*f + C.data + g.data;

		*myfile << t << "," << e[0] << "," << e[1] << "," << e[2] << "," << e[3] << "," << e[4] << "," << e[5] << "," << f[2] << "," << z_I << "," << pc[2] << std::endl;

		// std::cout << z_I << std::endl;
		std::cout << pc.transpose() << std::endl;
		std::cout << f[2] << std::endl;
		// std::cout << jnt.q.data.transpose() << std::endl;
		// std::cout << e.transpose() << std::endl;
		// std::cout << J.data.determinant() << std::endl;
		// std::cout << "---" << std::endl;
	}

};


void write_state(const state_type &qs, const double t)
{
	jsd.push_back(qs);
	std::cout << "Time: " << t << ", q: " << qs[0] << ", " << qs[1] << ", " << qs[2] << ", " <<qs[3] << ", " << qs[4] << ", " << qs[5] << "\n";
	std::cout << jsd.size() << "\n";
}

int main(int argc, char **argv) {

	// Get parameters for kinematic from setup.yaml file
	std::string urdf_file = "../ur5e/ur5e.urdf";
	std::string base_link = "base_link";
	std::string tool_link = "tool0";
	KDL::Vector gravity(0, 0, -0.82);

	// Generate kinematic model for orocos_kdl
	KDL::Tree tree;
	if (!kdl_parser::treeFromFile(urdf_file, tree)) {
		std::cout << "Failed to construct kdl tree\n" << std::endl;
		return -1;
	}
	KDL::Chain chain;
	if(!tree.getChain(base_link, tool_link, chain)){
		std::cout << "Failed to get KDL chain from tree\n" << std::endl;
		return 0;
	}
	std::cout << std::endl;

	double	q1 = 0.0,
			q2 = -1.12611,
			q3 =  1.94035,
			q4 = -1.33784,
			q5 = -M_PI/2,
			q6 = 0.0;
	double df0 = 0; // delta force integral
	// 0, -M_PI/2, M_PI/2, -M_PI/2, 0, 0
	state_type q0 = {q1,q2,q3,q5,q5,q6,  0,0,0,0,0,0,  0,0,0,0,0,0,  0,0,0,0,0,0, df0}; // initial conditions: q, dq, s, ds
	// std::ofstream file("data.csv");
	std::shared_ptr<std::ofstream> file = std::make_shared<std::ofstream>("data.csv");

	RobotState robot(chain, file);
	robot.set_q(q0);

	// boost::numeric::odeint::integrate(robot, q0, 0.0 , 20.0, dt);
	boost::numeric::odeint::integrate_const(boost::numeric::odeint::runge_kutta4<state_type>(), robot, q0, 0.0, 20.0, dt);

	file->close();
	std::cout << "OK" << std::endl;
}
