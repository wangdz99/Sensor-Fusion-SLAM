// Author:   Tong Qin               qintonguav@gmail.com
// 	         Shaozu Cao 		    saozu.cao@connect.ust.hk

//
// TODO: implement analytic Jacobians for LOAM residuals in this file
// 

#include <eigen3/Eigen/Dense>

//
// TODO: Sophus is ready to use if you have a good undestanding of Lie algebra.
// 
#include <sophus/so3.hpp>

#include <ceres/ceres.h>
#include <ceres/rotation.h>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl_conversions/pcl_conversions.h>


struct LidarEdgeFactor
{
	LidarEdgeFactor(Eigen::Vector3d curr_point_, Eigen::Vector3d last_point_a_,
					Eigen::Vector3d last_point_b_, double s_)
		: curr_point(curr_point_), last_point_a(last_point_a_), last_point_b(last_point_b_), s(s_) {}

	template <typename T>
	bool operator()(const T *q, const T *t, T *residual) const
	{

		Eigen::Matrix<T, 3, 1> cp{T(curr_point.x()), T(curr_point.y()), T(curr_point.z())};
		Eigen::Matrix<T, 3, 1> lpa{T(last_point_a.x()), T(last_point_a.y()), T(last_point_a.z())};
		Eigen::Matrix<T, 3, 1> lpb{T(last_point_b.x()), T(last_point_b.y()), T(last_point_b.z())};

		//Eigen::Quaternion<T> q_last_curr{q[3], T(s) * q[0], T(s) * q[1], T(s) * q[2]};
		Eigen::Quaternion<T> q_last_curr{q[3], q[0], q[1], q[2]};
		Eigen::Quaternion<T> q_identity{T(1), T(0), T(0), T(0)};
		q_last_curr = q_identity.slerp(T(s), q_last_curr);
		Eigen::Matrix<T, 3, 1> t_last_curr{T(s) * t[0], T(s) * t[1], T(s) * t[2]};

		Eigen::Matrix<T, 3, 1> lp;
		lp = q_last_curr * cp + t_last_curr;

		Eigen::Matrix<T, 3, 1> nu = (lp - lpa).cross(lp - lpb);
		Eigen::Matrix<T, 3, 1> de = lpa - lpb;

		residual[0] = nu.x() / de.norm();
		residual[1] = nu.y() / de.norm();
		residual[2] = nu.z() / de.norm();

		return true;
	}

	static ceres::CostFunction *Create(const Eigen::Vector3d curr_point_, const Eigen::Vector3d last_point_a_,
									   const Eigen::Vector3d last_point_b_, const double s_)
	{
		return (new ceres::AutoDiffCostFunction<
				LidarEdgeFactor, 3, 4, 3>(
			new LidarEdgeFactor(curr_point_, last_point_a_, last_point_b_, s_)));
	}

	Eigen::Vector3d curr_point, last_point_a, last_point_b;
	double s;
};

// Eigen::Matrix3d skew(Eigen::Vector3d& mat_in){
// 	Eigen::Matrix3d skew_mat;
// 	skew_mat.setZero();
// 	skew_mat(0,1) = -mat_in(2);
// 	skew_mat(0,2) =  mat_in(1);
// 	skew_mat(1,2) = -mat_in(0);
// 	skew_mat(1,0) =  mat_in(2);
// 	skew_mat(2,0) = -mat_in(1);
// 	skew_mat(2,1) =  mat_in(0);
// 	return skew_mat;
// }


class LidarEdgeCostFn : public ceres::SizedCostFunction<3, 4, 3> {
	public:
		LidarEdgeCostFn(Eigen::Vector3d curr_point_, Eigen::Vector3d last_point_a_,
						Eigen::Vector3d last_point_b_, double s_)
			: curr_point(curr_point_), last_point_a(last_point_a_), last_point_b(last_point_b_), s(s_){
				std::cout << s << std::endl;
			}
		virtual ~LidarEdgeCostFn() {}
		// friend Eigen::Matrix3d skew(Eigen::Vector3d& mat_in);
		virtual bool Evaluate(double const* const* parameters,
							double* residuals,
							double** jacobians) const {
		
		// std::cout << parameters[] <<std::endl;
		Eigen::Map<const Eigen::Quaternion<double>> q_last_curr(parameters[0]);
		Eigen::Map<const Eigen::Vector3d> t_last_curr(parameters[0] + 4);
		Eigen::Vector3d lp;
		lp = q_last_curr * curr_point + t_last_curr;

		Eigen::Vector3d nu = (lp - last_point_a).cross(lp - last_point_b);
		Eigen::Vector3d de = last_point_a - last_point_b;

		residuals[0] = nu.x() / de.norm();
		residuals[0] = nu.y() / de.norm();
		residuals[0] = nu.z() / de.norm();

		if (jacobians && jacobians[0]){
			Eigen::Matrix3d skew_lp; // = skew(lp);
			skew_lp.setZero();
			skew_lp(0,1) = -lp(2);
			skew_lp(0,2) =  lp(1);
			skew_lp(1,2) = -lp(0);
			skew_lp(1,0) =  lp(2);
			skew_lp(2,0) = -lp(1);
			skew_lp(2,1) =  lp(0); 
			Eigen::Matrix<double, 3, 6> dp_by_so3;
			dp_by_so3.block<3, 3>(0, 0) = -skew_lp;
			(dp_by_so3.block<3, 3>(0, 3)).setIdentity();
			Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>> J_se3(jacobians[0]);
			J_se3.setZero();
			Eigen::Vector3d re = last_point_b -last_point_a;
			Eigen::Matrix3d skew_re; // = skew(re);
			skew_re.setZero();
			skew_re(0,1) = -re(2);
			skew_re(0,2) =  re(1);
			skew_re(1,2) = -re(0);
			skew_re(1,0) =  re(2);
			skew_re(2,0) = -re(1);
			skew_re(2,1) =  re(0); 
			J_se3.block<3, 6>(0, 0) = skew_re * dp_by_so3 /de.norm();
		}
		
		return true;
		}

		private:
		Eigen::Vector3d curr_point, last_point_a, last_point_b;
		double s;
 };

struct LidarPlaneFactor
{
	LidarPlaneFactor(Eigen::Vector3d curr_point_, Eigen::Vector3d last_point_j_,
					 Eigen::Vector3d last_point_l_, Eigen::Vector3d last_point_m_, double s_)
		: curr_point(curr_point_), last_point_j(last_point_j_), last_point_l(last_point_l_),
		  last_point_m(last_point_m_), s(s_)
	{
		ljm_norm = (last_point_j - last_point_l).cross(last_point_j - last_point_m);
		ljm_norm.normalize();
	}

	template <typename T>
	bool operator()(const T *q, const T *t, T *residual) const
	{

		Eigen::Matrix<T, 3, 1> cp{T(curr_point.x()), T(curr_point.y()), T(curr_point.z())};
		Eigen::Matrix<T, 3, 1> lpj{T(last_point_j.x()), T(last_point_j.y()), T(last_point_j.z())};
		//Eigen::Matrix<T, 3, 1> lpl{T(last_point_l.x()), T(last_point_l.y()), T(last_point_l.z())};
		//Eigen::Matrix<T, 3, 1> lpm{T(last_point_m.x()), T(last_point_m.y()), T(last_point_m.z())};
		Eigen::Matrix<T, 3, 1> ljm{T(ljm_norm.x()), T(ljm_norm.y()), T(ljm_norm.z())};

		//Eigen::Quaternion<T> q_last_curr{q[3], T(s) * q[0], T(s) * q[1], T(s) * q[2]};
		Eigen::Quaternion<T> q_last_curr{q[3], q[0], q[1], q[2]};
		Eigen::Quaternion<T> q_identity{T(1), T(0), T(0), T(0)};
		q_last_curr = q_identity.slerp(T(s), q_last_curr);
		Eigen::Matrix<T, 3, 1> t_last_curr{T(s) * t[0], T(s) * t[1], T(s) * t[2]};

		Eigen::Matrix<T, 3, 1> lp;
		lp = q_last_curr * cp + t_last_curr;

		residual[0] = (lp - lpj).dot(ljm);

		return true;
	}

	static ceres::CostFunction *Create(const Eigen::Vector3d curr_point_, const Eigen::Vector3d last_point_j_,
									   const Eigen::Vector3d last_point_l_, const Eigen::Vector3d last_point_m_,
									   const double s_)
	{
		return (new ceres::AutoDiffCostFunction<
				LidarPlaneFactor, 1, 4, 3>(
			new LidarPlaneFactor(curr_point_, last_point_j_, last_point_l_, last_point_m_, s_)));
	}

	Eigen::Vector3d curr_point, last_point_j, last_point_l, last_point_m;
	Eigen::Vector3d ljm_norm;
	double s;
};

struct LidarPlaneNormFactor
{

	LidarPlaneNormFactor(Eigen::Vector3d curr_point_, Eigen::Vector3d plane_unit_norm_,
						 double negative_OA_dot_norm_) : curr_point(curr_point_), plane_unit_norm(plane_unit_norm_),
														 negative_OA_dot_norm(negative_OA_dot_norm_) {}

	template <typename T>
	bool operator()(const T *q, const T *t, T *residual) const
	{
		Eigen::Quaternion<T> q_w_curr{q[3], q[0], q[1], q[2]};
		Eigen::Matrix<T, 3, 1> t_w_curr{t[0], t[1], t[2]};
		Eigen::Matrix<T, 3, 1> cp{T(curr_point.x()), T(curr_point.y()), T(curr_point.z())};
		Eigen::Matrix<T, 3, 1> point_w;
		point_w = q_w_curr * cp + t_w_curr;

		Eigen::Matrix<T, 3, 1> norm(T(plane_unit_norm.x()), T(plane_unit_norm.y()), T(plane_unit_norm.z()));
		residual[0] = norm.dot(point_w) + T(negative_OA_dot_norm);
		return true;
	}

	static ceres::CostFunction *Create(const Eigen::Vector3d curr_point_, const Eigen::Vector3d plane_unit_norm_,
									   const double negative_OA_dot_norm_)
	{
		return (new ceres::AutoDiffCostFunction<
				LidarPlaneNormFactor, 1, 4, 3>(
			new LidarPlaneNormFactor(curr_point_, plane_unit_norm_, negative_OA_dot_norm_)));
	}

	Eigen::Vector3d curr_point;
	Eigen::Vector3d plane_unit_norm;
	double negative_OA_dot_norm;
};


class LidarPlaneNormCostFn : public ceres::SizedCostFunction<1, 4, 3> {
	public:

	LidarPlaneNormCostFn(Eigen::Vector3d curr_point_, Eigen::Vector3d plane_unit_norm_, double negative_OA_dot_norm_) 
			: curr_point(curr_point_), plane_unit_norm(plane_unit_norm_), negative_OA_dot_norm(negative_OA_dot_norm_) {}
	virtual ~LidarPlaneNormCostFn() {}
	virtual bool Evaluate(double const* const* parameters,
						double* residuals,
						double** jacobians) const {
					
		Eigen::Map<const Eigen::Quaternion<double>> q_w_curr(parameters[0]);
		Eigen::Map<const Eigen::Vector3d> t_w_curr(parameters[0] + 4);
		Eigen::Vector3d point_w = q_w_curr * curr_point + t_w_curr;

		residuals[0] = plane_unit_norm.dot(point_w) + negative_OA_dot_norm;

		if (jacobians && jacobians[0]){
			Eigen::Matrix3d skew_point_w;   //= skew(point_w);
			skew_point_w.setZero();
			skew_point_w(0,1) = -point_w(2);
			skew_point_w(0,2) =  point_w(1);
			skew_point_w(1,2) = -point_w(0);
			skew_point_w(1,0) =  point_w(2);
			skew_point_w(2,0) = -point_w(1);
			skew_point_w(2,1) =  point_w(0); 
			Eigen::Matrix<double, 3, 6> dp_by_so3;
			dp_by_so3.block<3, 3>(0, 0) = -skew_point_w;
			(dp_by_so3.block<3, 3>(0, 3)).setIdentity();
			Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor>> J_se3(
				jacobians[0]);
			J_se3.setZero();
			J_se3.block<1, 6>(0, 0) = plane_unit_norm.transpose() * dp_by_so3;
		}

		return true;
	
	}

	private:
	Eigen::Vector3d curr_point;
	Eigen::Vector3d plane_unit_norm;
	double negative_OA_dot_norm;
	

};

struct LidarDistanceFactor
{

	LidarDistanceFactor(Eigen::Vector3d curr_point_, Eigen::Vector3d closed_point_) 
						: curr_point(curr_point_), closed_point(closed_point_){}

	template <typename T>
	bool operator()(const T *q, const T *t, T *residual) const
	{
		Eigen::Quaternion<T> q_w_curr{q[3], q[0], q[1], q[2]};
		Eigen::Matrix<T, 3, 1> t_w_curr{t[0], t[1], t[2]};
		Eigen::Matrix<T, 3, 1> cp{T(curr_point.x()), T(curr_point.y()), T(curr_point.z())};
		Eigen::Matrix<T, 3, 1> point_w;
		point_w = q_w_curr * cp + t_w_curr;


		residual[0] = point_w.x() - T(closed_point.x());
		residual[1] = point_w.y() - T(closed_point.y());
		residual[2] = point_w.z() - T(closed_point.z());
		return true;
	}

	static ceres::CostFunction *Create(const Eigen::Vector3d curr_point_, const Eigen::Vector3d closed_point_)
	{
		return (new ceres::AutoDiffCostFunction<
				LidarDistanceFactor, 3, 4, 3>(
			new LidarDistanceFactor(curr_point_, closed_point_)));
	}

	Eigen::Vector3d curr_point;
	Eigen::Vector3d closed_point;
};