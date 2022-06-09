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

	//��ȥ�������
	Eigen::VectorXd strength;

	//��ȥ�������
	Eigen::VectorXd slipStrength;

	//��ȥ�������
	Eigen::MatrixXd A;

	//��ȥ�������
	Eigen::MatrixXd B;

	Array<Vector2D> newParticles;

	//tracer���ӵ��ٶȺ�λ��
	Array<Vector2D> tracePosition;
	Array<Vector2D> traceVelocity;

private:

	//��ά����£������Ǳ���
	ArrayD _vorticities;


};

using VortexParticleDataPtr = std::shared_ptr<VortexParticleData>;

inline ArrayD& VortexParticleData::vorticities() {
	return _vorticities;
}
