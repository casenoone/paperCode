#pragma once

#include <Eigen/Dense>
#include "../titmouse2d/src/Lagrangian/ParticleSystemData2.h"
#include "../titmouse2d/src/Geometry/RegularPolygon.h"


class VortexParticleData : public ParticleSystemData2 {
public:
	VortexParticleData() = default;

	ArrayD& vorticities();

public:
	RegularPolygonPtr panelSet;


	//Array<Panel> panelSet;

	//消去法向分量
	Eigen::VectorXd strength;

	//消去切向分量
	Eigen::VectorXd slipStrength;

	//消去法向分量
	Eigen::MatrixXd A;

	//消去切向分量
	Eigen::MatrixXd B;

	Array<Vector2D> newParticles;

	//tracer粒子的速度和位置
	Array<Vector2D> tracePosition;
	Array<Vector2D> traceVelocity;

private:

	//二维情况下，涡量是标量
	ArrayD _vorticities;


};

using VortexParticleDataPtr = std::shared_ptr<VortexParticleData>;

inline ArrayD& VortexParticleData::vorticities() {
	return _vorticities;
}
