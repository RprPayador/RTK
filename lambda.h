#ifndef _LAMBDA_H_
#define _LAMBDA_H_

#include "RTK_Structs.h"
#include <Eigen/Dense>

/* LAMBDA 整数模糊度固定与基线恢复接口 */
bool FixAmbiguity(RAWDAT* Raw, PPRESULT* Base, PPRESULT* Rov,
                  const Eigen::Vector3d& b_hat,
                  const Eigen::VectorXd& a_hat,
                  const Eigen::Matrix3d& Qbb,
                  const Eigen::MatrixXd& Qaa,
                  const Eigen::MatrixXd& Qba,
                  double ratioThres);

#endif
