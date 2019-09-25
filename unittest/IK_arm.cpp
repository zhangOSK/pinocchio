
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;

#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include <urdf_parser/urdf_parser.h>

int main(int /* argc */, char ** /* argv */)
{
  pinocchio::Model model;
  //const std::string filename = "/opt/openrobots/share/talos_data/urdf/talos_arm.urdf";
  const std::string filename = "/home/ang/Downloads/talos_data/urdf/talos_arm.urdf";
  pinocchio::urdf::buildModel(filename, model);
  pinocchio::Data data(model);
  const int    JOINT_ID = 7;
  Eigen::Vector3d xdes;        xdes << 0.048, 0.3, -0.3; 

  //Eigen::VectorXd q     = pinocchio::neutral(model);
  Eigen::VectorXd q(7);
  q << 0.25847, 0.173046, -0.0002, -0.525366, 0, 0, 0.1; 
  const double eps      = 1e-3;
  const int IT_MAX      = 1000;
  const double DT       = 5e-3;

  pinocchio::Data::Matrix6x J(6,model.nv); J.setZero();
  unsigned int svdOptions = Eigen::ComputeThinU | Eigen::ComputeThinV;
  Eigen::BDCSVD<pinocchio::Data::Matrix3x> svd(3,model.nv,svdOptions);
  Eigen::Vector3d err;

  /*
  //--------------------------get init lw------------------------------
  pinocchio::forwardKinematics(model,data,q);
  const Eigen::Vector3d & x_base   = data.oMi[7].translation();
  const Eigen::Matrix3d & R_base   = data.oMi[7].rotation();
  Eigen::Vector3d base111;
  base111 = R_base.transpose()*x_base;
  std::cout<<"\n hand7 init= \n "<< base111(0) << ", "<< base111(1) << ", "<< base111(2) << ", "<< std::endl;
  */
  // Init lw Pos: 0.0459326, 0.332457, -0.335712

  for (int i=0;;i++)
  {
    pinocchio::forwardKinematics(model,data,q);
    const Eigen::Vector3d & x   = data.oMi[JOINT_ID].translation();
    const Eigen::Matrix3d & R   = data.oMi[JOINT_ID].rotation();
    err = R.transpose()*(x)-xdes; //! xdes is given in shoulder frame
    
    if(err.norm() < eps)
    {
      std::cout << "Convergence achieved!" << std::endl;
      break;
    }
    if (i >= IT_MAX)
    {
      std::cout << "\nWarning: the iterative algorithm has not reached convergence to the desired precision" << std::endl;
      break;
    }
    pinocchio::jointJacobian(model,data,q,JOINT_ID,J);
    //------------update gain-----------------------
    double norm =  err.norm();
    double coeff_a = 0.1;
    double coeff_b = 0.1;
    double coeff_c = 125e3;
    double res = coeff_a * exp( -coeff_b*norm ) + coeff_c;
    //-----------------------------------------------
    const Eigen::VectorXd v     = - svd.compute(J.topRows<3>()).solve(err*res); //adaptive gain
    q = pinocchio::integrate(model,q,v*DT);
    if(!(i%10)) std::cout << "error = " << err.transpose()  
                          << " norm L2:" << err.transpose()*err << std::endl;
  }

  //std::cout << "\nresult: " << q.transpose() << std::endl;
  //std::cout << "\nfinal error: " << err.transpose() 
    //        << " norm L2:" << err.transpose()*err << std::endl;

  const Eigen::Vector3d & x_arm0   = data.oMi[7].translation();
  const Eigen::Matrix3d & R_arm0   = data.oMi[7].rotation();
  Eigen::Vector3d base = R_arm0.transpose()*x_arm0;
  std::cout <<"Init lw = 0.0459326, 0.332457, -0.335712"<<std::endl;
  std::cout <<"xdes = "<<xdes(0)<<", "<<xdes(1)<<", "<<xdes(2)<<std::endl;
  std::cout<<"\n hand final = "<< base(0) << ", "<<base(1)<<", "<<base(2) <<std::endl;
  
  return 0;
  
}
