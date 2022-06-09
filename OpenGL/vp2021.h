#pragma once


#include "vp2021Data.h"
#include "../titmouse2d/src/Lagrangian/ParticleSystemSolver2.h"
#include "../titmouse2d/src/OtherMethod/SWE/ShallowWaveSolver2.h"
#include "../titmouse2d/src/boundingbox2.h"

#include <Eigen/Dense>


//����涨��
//x : pos.x
//y : height
//z : pos.y


class Vp2021Solver : public ParticleSystemSolver2 {

public:
	Vp2021Solver();

	//ʱ�����
	virtual void timeIntegration(double timeIntervalInSeconds)override;

	//ʱ������
	virtual void onAdvanceTimeStep(double timeIntervalInSeconds)override;

	Vp2021DataPtr& vp2021Data();

	//���������Ӷ��������ӵ��յ��ٶ�
	Vector2D computeUSingle(const Vector2D& pos, int i)const;

	//���������Ӷ��������ӵ��յ��ٶ�
	Vector2D computeUSingle1(const Vector2D& pos, int i)const;

	//�����ƶ��߽磨ʹ��siggraph2020�ı߽���ⷽ����
	void setMovingBoudnary(RegularPolygonPtr surfaces);

	//����׷������
	void emitTracerParticles();

	//���ƶ��߽��Ϸ���������������no-slip�߽�����
	void emitParticlesFromPanels(double timeIntervalInSeconds);

	//����vortex ring
	void emitVortexRing();

private:
	//���tracer����
	void tracerParticlesSolve();

	//ֻҪ�߽���״���䣬�߽����Ͳ����
	//�ƶ��߽����߽���ַ���ʱʹ�õľ���SIGGRAPH2020��
	void constructMovingBoundaryMatrix();

	//��װ�ƶ��߽��Ӧ�ľ�������
	Vector2D computeUnitVelocityFromPanels(const Vector2D& pos, int index);

	//���׷�����ӵı߽紩͸����
	void tarcerCollisionSolve(Vector2D& pos);

private:
	Vp2021DataPtr _vp2021Data;
};




