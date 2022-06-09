#include "vp2021.h"

#include "../titmouse2d/src/ConstVar.h"
#include "../titmouse2d/src/random.h"
#include "../titmouse2d/src/boundingbox2.h"
#include <iostream>
#include <cmath>

//每次每个边界控制点上发射的粒子数
static const int numOfStep = 5;

static const double vorticity_eps = 0.08;
static double fv_eps = 0.001;

static double temp_r = 0.25;

Vp2021Solver::Vp2021Solver() {
	_particleSystemData = std::make_shared<Vp2021Data>();
	_vp2021Data = std::make_shared<Vp2021Data>();
	_vp2021Data = std::dynamic_pointer_cast<Vp2021Data>(_particleSystemData);
}

//这里额外处理
void Vp2021Solver::timeIntegration(double timeIntervalInSeconds) {

	auto& vortex_pos = _vp2021Data->vortexPosition;
	auto& vortex_vel = _vp2021Data->vortexVelocity;
	auto n = vortex_pos.dataSize();
	Array<Vector2D> tempP(n);

	//求解涡粒子的更新
	for (int i = 0; i < n; ++i) {
		Vector2D vortexVel;
		for (int j = 0; j < n; ++j) {
			if (i != j) {
				//vortexVel += computeUSingle(vortex_pos[i], j);
				vortexVel += computeUSingle1(vortex_pos[i], j);
			}
		}
		Vector2D votexSheetVel;
		tempP[i] = vortex_pos[i] + timeIntervalInSeconds * (votexSheetVel + vortexVel);
	}
	vortex_pos = tempP;


	//求解tracer粒子的速度场
	tracerParticlesSolve();


	//在这里更新tracer粒子
	auto tracerPos = _vp2021Data->tracePosition;
	auto traceVel = _vp2021Data->traceVelocity;
	auto tracer_n = tracerPos.dataSize();
	for (int i = 0; i < tracer_n; ++i) {
		tracerPos[i] += (traceVel[i]) * timeIntervalInSeconds;
		tarcerCollisionSolve(tracerPos[i]);
	}
}


void Vp2021Solver::onAdvanceTimeStep(double timeIntervalInSeconds) {
	beginAdvanceTimeStep();
	timeIntegration(timeIntervalInSeconds);
	emitParticlesFromPanels(timeIntervalInSeconds);
	endAdvanceTimeStep();
}


Vp2021DataPtr& Vp2021Solver::vp2021Data() {
	return _vp2021Data;
}

//2005
Vector2D Vp2021Solver::computeUSingle(const Vector2D& pos, int i)const {

	auto position = _vp2021Data->vortexPosition;
	auto gamma = _vp2021Data->gamma();
	auto r2 = (pos - position[i]).getLengthSquared();
	auto uv = Vector2D(position[i].y - pos.y, pos.x - position[i].x);
	return gamma[i] * uv / (kPiD * r2) * 0.5 * (1.0 - pow(kE, -r2 / (vorticity_eps * vorticity_eps)));
}

//2021
Vector2D Vp2021Solver::computeUSingle1(const Vector2D& pos, int i)const {

	auto position = _vp2021Data->vortexPosition;
	auto gamma = _vp2021Data->gamma();
	auto r2 = (pos - position[i]).getLengthSquared();
	auto uv = Vector2D(position[i].y - pos.y, pos.x - position[i].x);
	return gamma[i] * uv / (kPiD * r2 + temp_r * temp_r);
}


void Vp2021Solver::setMovingBoudnary(RegularPolygonPtr surfaces) {
	_vp2021Data->panelSet = surfaces;
	constructMovingBoundaryMatrix();
}


void Vp2021Solver::emitVortexRing() {
	auto& pos = _vp2021Data->vortexPosition;
	auto& vel = _vp2021Data->vortexVelocity;
	auto& gamma = _vp2021Data->gamma();
	int n = 4;
	pos.reSize(n);
	vel.reSize(n);
	gamma.reSize(n);
	Vector2D A(0.2, 1.2);
	Vector2D B(0.2, 1.1);
	Vector2D C(0.2, 1.0);
	Vector2D D(0.2, 0.9);

	pos[0] = A;
	pos[1] = B;
	pos[2] = C;
	pos[3] = D;

	gamma[0] = 0.6;
	gamma[1] = 0.6;
	gamma[2] = -0.6;
	gamma[3] = -0.6;
}


void Vp2021Solver::emitParticlesFromPanels(double timeIntervalInSeconds) {

	static int step = 0;
	auto data = _vp2021Data;
	auto& pos = _vp2021Data->vortexPosition;
	auto panelVel = data->panelSet->velocity;
	auto panel = data->panelSet;
	auto& emitParticle = data->newParticles;

	auto emitNum = emitParticle.dataSize();

	//更改粒子发射位置使得与panel同步
	for (int i = 0; i < emitNum; ++i) {
		emitParticle[i] += panelVel * timeIntervalInSeconds;
	}

	if (step % 4 == 0 && pos.dataSize() < 100000 && (panel->center().x - panel->r()) < 2.0) {
		auto& pos = data->vortexPosition;
		auto& vel = data->vortexVelocity;
		auto panels = data->panelSet;
		auto gamma = data->slip_strength;

		Vector2D tempPos;

		auto startNum = pos.dataSize();

		for (int i = 0; i < emitNum; ++i) {
			pos.push(emitParticle(i));
			vel.push(Vector2D(0.0, 0.0));
			double temp_gamma = 0;
			data->gamma().push(temp_gamma);
		}


		Eigen::MatrixXd& B = _vp2021Data->slip_matrix;
		Eigen::VectorXd& x = _vp2021Data->slip_strength;

		auto panelSize = data->panelSet->size();
		auto panelVelocity = data->panelSet->velocity;
		//组装b
		Eigen::VectorXd b(panelSize * 2);

		for (int i = 0; i < panelSize; ++i) {
			auto temp_pos = panels->midPoint(i);
			Vector2D tempVel;
			for (int j = 0; j < pos.dataSize(); ++j) {
				//tempVel += computeUSingle(temp_pos, j);
				tempVel += computeUSingle1(temp_pos, j);
			}
			auto v1 = (panelVelocity - tempVel);
			b[i] = v1.x;
			b[i + panelSize] = v1.y;
		}

		x = B * b;
		auto& vor = data->gamma();
		int j = 0;
		for (int i = startNum; i < pos.dataSize(); ++i) {
			vor[i] = x[j++];
		}
	}

	step++;
}




//计算在pos处引发的诱导速度
Vector2D Vp2021Solver::computeUnitVelocityFromPanels(const Vector2D& pos, int index) {

	auto position = _vp2021Data->newParticles;
	auto r2 = (pos - position[index]).getLengthSquared();
	auto uv = Vector2D(position[index].y - pos.y, pos.x - position[index].x);
	return uv / (kPiD * r2) * 0.5 * (1.0 - pow(kE, -r2 / (vorticity_eps * vorticity_eps)));
}


//移动边界所使用的tracer粒子
void Vp2021Solver::emitTracerParticles() {
	static int step = 0;
	auto data = _vp2021Data;
	auto& tracerPos = _vp2021Data->tracePosition;
	auto& tracerVel = _vp2021Data->traceVelocity;
	auto n = tracerPos.dataSize();
	auto panels = data->panelSet;

	int emitNum = 1000;
	//int emitNum = 1;
	tracerPos.reSize(emitNum);
	tracerVel.reSize(emitNum);
	Vector2D tempPos;
	for (int i = 0; i < emitNum; ++i) {
		tempPos.x = random_double(0.2, 1.0);
		tempPos.y = random_double(0.5, 1.5);
		tracerPos[i] = tempPos;
	}
}


//求解tracer粒子的运动
void Vp2021Solver::tracerParticlesSolve() {
	auto& tracerPos = _vp2021Data->tracePosition;
	auto& tracerVel = _vp2021Data->traceVelocity;
	auto n = tracerPos.dataSize();

	auto vor_n = _vp2021Data->vortexPosition.dataSize();
	for (int i = 0; i < n; ++i) {
		Vector2D tempVel;
		Vector2D pos = tracerPos(i);
		for (int j = 0; j < vor_n; ++j) {
			//tempVel += computeUSingle(pos, j);
			tempVel += computeUSingle1(pos, j);
		}
		tracerVel(i) = tempVel;
	}
}





//只执行一次
void Vp2021Solver::constructMovingBoundaryMatrix() {
	auto panels = _vp2021Data->panelSet;
	auto& emittedParticles = _vp2021Data->newParticles;
	auto panelSize = panels->size();
	auto emitNum = panelSize * numOfStep;
	emittedParticles.reSize(emitNum);

	Vector2D tempPos;
	for (int i = 0; i < panelSize; ++i) {
		for (int j = 0; j < numOfStep; ++j) {
			auto line = panels->lookAt(i).end - panels->lookAt(i).start;
			auto lambda = random_double(0, 1);
			tempPos = panels->lookAt(i).start + lambda * line;
			tempPos += panels->lookAt(i).normal * random_double(0.01, 0.02);
			emittedParticles(i * numOfStep + j) = tempPos;
		}
	}

	Eigen::MatrixXd A;

	A.resize(panelSize * 2, emitNum);
	for (int i = 0; i < panelSize; ++i) {
		auto pos = panels->midPoint(i);
		for (int j = 0; j < emitNum; ++j) {
			auto vel = computeUnitVelocityFromPanels(pos, j);
			A(i, j) = vel.x;
			A(i + panelSize, j) = vel.y;
		}
	}



	Eigen::MatrixXd& B = _vp2021Data->slip_matrix;
	Eigen::MatrixXd A_trans = A.transpose();

	auto I = Eigen::MatrixXd::Identity(emitNum, emitNum);
	B = (A_trans * A + 3 * I).inverse() * A_trans;
}

void Vp2021Solver::tarcerCollisionSolve(Vector2D& pos) {
	if (pos.x > 2)
		pos.x = 2;
	if (pos.x < 0)
		pos.x = 0;
	if (pos.y > 2)
		pos.y = 2;
	if (pos.y < 0)
		pos.y = 0;
}



