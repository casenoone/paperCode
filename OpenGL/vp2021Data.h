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

	//2021 paper �߽紦��
	//��ȥ�������
	Eigen::VectorXd slip_strength;

	//��ȥ�������
	Eigen::MatrixXd slip_matrix;

	//���淢�����ӵ�λ��
	Array<Vector2D> newParticles;

	//tracer���ӵ��ٶȺ�λ��
	Array<Vector2D> tracePosition;
	Array<Vector2D> traceVelocity;

	//��ŷ����ȥ�������ӵ���Ϣ
	Array<Vector2D> vortexPosition;
	Array<Vector2D> vortexVelocity;


private:

	//��ά����£������Ǳ���
	ArrayD _gamma;


};

using Vp2021DataPtr = std::shared_ptr<Vp2021Data>;

inline ArrayD& Vp2021Data::gamma() {
	return _gamma;
}

#endif
