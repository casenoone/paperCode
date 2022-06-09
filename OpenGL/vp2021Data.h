#ifndef FOAMVORTEXDATA_H
#define FOAMVORTEXDATA_H

#include <Eigen/Dense>
#include "../titmouse2d/src/Lagrangian/ParticleSystemData2.h"

#include "../titmouse2d/src/Geometry/RegularPolygon.h"

class Vp2021Data : public ParticleSystemData2 {
public:
	Vp2021Data() = default;

	ArrayD& gamma();

public:

	RegularPolygonPtr panelSet;

	//2021 paper 边界处理
	//消去切向分量
	Eigen::VectorXd slip_strength;

	//消去切向分量
	Eigen::MatrixXd slip_matrix;

	//保存发射粒子的位置
	Array<Vector2D> newParticles;

	//tracer粒子的速度和位置
	Array<Vector2D> tracePosition;
	Array<Vector2D> traceVelocity;

	//存放发射出去的涡粒子的信息
	Array<Vector2D> vortexPosition;
	Array<Vector2D> vortexVelocity;


private:

	//二维情况下，涡量是标量
	ArrayD _gamma;


};

using Vp2021DataPtr = std::shared_ptr<Vp2021Data>;

inline ArrayD& Vp2021Data::gamma() {
	return _gamma;
}

#endif
