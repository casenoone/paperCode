#pragma once


#include "vp2005Data.h"
#include "../titmouse2d/src/Lagrangian/ParticleSystemSolver2.h"

#include "../titmouse2d/src/boundingbox2.h"


#include <Eigen/Dense>


class VortexParticleSolver : public ParticleSystemSolver2 {

public:
	VortexParticleSolver();

	virtual ~VortexParticleSolver();

	virtual void timeIntegration(double timeIntervalInSeconds)override;


	virtual void onAdvanceTimeStep(double timeIntervalInSeconds)override;

	VortexParticleDataPtr& vortexParticleData();

	Vector2D computeUSingle(const Vector2D& pos, int i)const;

	void setData(int numberOfParticles,
		Array<Vector2D>& pos,
		int resolutionX,
		int resolutionY);

	void setPanels(RegularPolygonPtr surfaces);

	void emitParticles();

	void emitParticlesFromPanels(double timeIntervalInSeconds);



private:

	//求解tracer粒子
	void tracerParticlesSolve();

	Vector2D computeUnitVelocityFromPanels(int index, const Vector2D& midPoint);

	//只要边界形状不变，边界矩阵就不会变
	//这个函数只调用一次
	void computeBoundaryMatrix();

	//我也不想起这么长的名字
	Vector2D computeUnitVelocityFromPanels(const Vector2D& pos, int index);

	Vector2D computeSingleVelocityFromPanels(int index);

	void vortexSheetSolve(double timeIntervalInSeconds);

	void slipVortexSheetSolve(double timeIntervalInSeconds);

	void onBeginAdvanceTimeStep();

	void onEndAdvanceTimeStep();


private:
	VortexParticleDataPtr _vortexParticleData;

};


