#pragma once


#include "vp2021Data.h"
#include "../titmouse2d/src/Lagrangian/ParticleSystemSolver2.h"
#include "../titmouse2d/src/OtherMethod/SWE/ShallowWaveSolver2.h"
#include "../titmouse2d/src/boundingbox2.h"

#include <Eigen/Dense>


//坐标规定：
//x : pos.x
//y : height
//z : pos.y


class Vp2021Solver : public ParticleSystemSolver2 {

public:
	Vp2021Solver();

	//时间积分
	virtual void timeIntegration(double timeIntervalInSeconds)override;

	//时步更新
	virtual void onAdvanceTimeStep(double timeIntervalInSeconds)override;

	Vp2021DataPtr& vp2021Data();

	//计算涡粒子对其他粒子的诱导速度
	Vector2D computeUSingle(const Vector2D& pos, int i)const;

	//计算涡粒子对其他粒子的诱导速度
	Vector2D computeUSingle1(const Vector2D& pos, int i)const;

	//设置移动边界（使用siggraph2020的边界求解方法）
	void setMovingBoudnary(RegularPolygonPtr surfaces);

	//发射追踪粒子
	void emitTracerParticles();

	//从移动边界上发射涡粒子以满足no-slip边界条件
	void emitParticlesFromPanels(double timeIntervalInSeconds);

	//发射vortex ring
	void emitVortexRing();

private:
	//求解tracer粒子
	void tracerParticlesSolve();

	//只要边界形状不变，边界矩阵就不会变
	//移动边界求解边界积分方程时使用的矩阵（SIGGRAPH2020）
	void constructMovingBoundaryMatrix();

	//组装移动边界对应的矩阵所用
	Vector2D computeUnitVelocityFromPanels(const Vector2D& pos, int index);

	//求解追踪粒子的边界穿透问题
	void tarcerCollisionSolve(Vector2D& pos);

private:
	Vp2021DataPtr _vp2021Data;
};




