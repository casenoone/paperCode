#include "vp2005.h"
#include "../titmouse2d/src/ConstVar.h"
#include "../titmouse2d/src/random.h"
#include "../titmouse2d/src/boundingbox2.h"
#include "../titmouse2d/src/LinearSystem/JacobiSolver.hpp"
#include "../titmouse2d/src/LinearSystem/ConjugateGradientSolver.hpp"

#include <iostream>
#include <cmath>
#include <windows.h>
#include <Eigen/SVD>
#include <Eigen/LU>



Vector2D vel_to_world(const Vector2D vel_local, const Vector2D n_i, const Vector2D normal_i) {
	double cos_theta = normal_i.y;
	double sin_theta = normal_i.x;
	Vector2D vel = Vector2D(cos_theta * vel_local.x + sin_theta * vel_local.y,
		-sin_theta * vel_local.x + cos_theta * vel_local.y);
	return vel;
}


const double vorticity_eps = 0.08;

const Vector2D vs_vec = Vector2D(5, 0.0);

double vs_tau = 40.0;

double fv_eps = 0.00001;

//newPosition赋值给positions可以改用移动构造函数，提高效率

VortexParticleSolver::VortexParticleSolver() {
	_particleSystemData = std::make_shared<VortexParticleData>();
	_vortexParticleData = std::make_shared<VortexParticleData>();
	_vortexParticleData = std::dynamic_pointer_cast<VortexParticleData>(_particleSystemData);
}

VortexParticleSolver::~VortexParticleSolver() {

}

void VortexParticleSolver::timeIntegration(double timeIntervalInSeconds) {

	auto n = _vortexParticleData->numberOfParticles();
	auto& positions = _vortexParticleData->positions();
	auto& vel = _vortexParticleData->velocities();

	for (int i = 0; i < n; ++i) {
		Vector2D vortexVel;
		for (int j = 0; j < n; ++j) {
			if (i != j) {
				vortexVel += computeUSingle(positions[i], j);
			}
		}
		Vector2D votexSheetVel;
		votexSheetVel = computeSingleVelocityFromPanels(i);

		_newVelocities[i] = vs_vec + votexSheetVel;
		_newPositions[i] = positions[i] + timeIntervalInSeconds * (_newVelocities[i] + vortexVel);
	}

	tracerParticlesSolve();

	//在这里更新tracer粒子
	auto tracerPos = _vortexParticleData->tracePosition;
	auto traceVel = _vortexParticleData->traceVelocity;
	auto tracer_n = tracerPos.dataSize();

	for (int i = 0; i < tracer_n; ++i) {
		tracerPos[i] += traceVel[i] * timeIntervalInSeconds;
	}

}


//求解tracer粒子的运动
void VortexParticleSolver::tracerParticlesSolve() {
	auto& tracerPos = _vortexParticleData->tracePosition;
	auto& tracerVel = _vortexParticleData->traceVelocity;
	auto n = tracerPos.dataSize();

	auto vor_n = _vortexParticleData->numberOfParticles();
	for (int i = 0; i < n; ++i) {
		Vector2D tempVel;
		Vector2D pos = tracerPos(i);
		for (int j = 0; j < vor_n; ++j) {
			tempVel += computeUSingle(pos, j);
		}
		tracerVel(i) = tempVel;
	}
}

void VortexParticleSolver::onAdvanceTimeStep(double timeIntervalInSeconds) {

	//消去法向分量
	vortexSheetSolve(timeIntervalInSeconds);
	timeIntegration(timeIntervalInSeconds);
	onEndAdvanceTimeStep();
	emitParticlesFromPanels(timeIntervalInSeconds);
	//消去切向分量
	slipVortexSheetSolve(timeIntervalInSeconds);
}


VortexParticleDataPtr& VortexParticleSolver::vortexParticleData() {
	return _vortexParticleData;
}


void VortexParticleSolver::onBeginAdvanceTimeStep() {
	ParticleSystemSolver2::beginAdvanceTimeStep();
}

void VortexParticleSolver::onEndAdvanceTimeStep() {
	ParticleSystemSolver2::endAdvanceTimeStep();
}



Vector2D VortexParticleSolver::computeUSingle(const Vector2D& pos, int i)const {

	auto position = _vortexParticleData->positions();
	auto gamma = _vortexParticleData->vorticities();
	auto r2 = (pos - position[i]).getLengthSquared();
	auto uv = Vector2D(position[i].y - pos.y, pos.x - position[i].x);
	return gamma[i] * uv / (kPiD * r2) * 0.5 * (1.0 - pow(kE, -r2 / (vorticity_eps * vorticity_eps)));
}

void VortexParticleSolver::setData(int numberOfParticles,
	Array<Vector2D>& pos,
	int resolutionX,
	int resolutionY) {

	ParticleSystemSolver2::setData(numberOfParticles, pos, resolutionX, resolutionY);
	_vortexParticleData->vorticities().reSize(numberOfParticles);
	auto vorticity = _vortexParticleData->vorticities();
}

void VortexParticleSolver::setPanels(RegularPolygonPtr surfaces) {
	_vortexParticleData->panelSet = surfaces;
	computeBoundaryMatrix();
}


void VortexParticleSolver::emitParticles() {
	auto data = _vortexParticleData;
	auto& vel = data->velocities();
	auto& n = data->numberOfParticles();
	auto& pos = data->positions();
	auto& gamma = data->vorticities();
	data->numberOfParticles() = 4;

	_newVelocities.reSize(n);
	_newPositions.reSize(n);
	vel.reSize(n);


	Vector2D A(0.2, 1.2);
	Vector2D B(0.2, 1.1);
	Vector2D C(0.2, 1.0);
	Vector2D D(0.2, 0.9);

	pos.push(A);
	pos.push(B);
	pos.push(C);
	pos.push(D);

	gamma.push(1.0);
	gamma.push(1.0);
	gamma.push(-1.0);
	gamma.push(-1.0);

	int tracerNum = 1000;
	data->tracePosition.reSize(tracerNum);
	data->traceVelocity.reSize(tracerNum);

	for (int i = 0; i < tracerNum; ++i) {
		Vector2D tempPos;
		tempPos.x = random_double(0, 0.4);
		tempPos.y = random_double(0.7, 1.4);
		data->tracePosition[i] = tempPos;
	}
}


void VortexParticleSolver::emitParticlesFromPanels(double timeIntervalInSeconds) {
	auto data = _vortexParticleData;
	auto n = data->numberOfParticles();
	auto& pos = data->positions();
	auto& vel = data->velocities();
	auto panels = data->panelSet;

	auto gamma = data->slipStrength;
	int emitNum = panels->size();
	Vector2D tempPos;
	static int step = 0;

	if (step % 3 == 0) {
		for (int i = 0; i < emitNum; ++i) {

			auto line = panels->lookAt(i).end - panels->lookAt(i).start;
			double panelLength = line.getLength();
			auto lambda = random_double(0, 1);
			//tempPos = panels->lookAt(i).start + lambda * line;
			tempPos = panels->midPoint(i);
			tempPos += panels->lookAt(i).normal * random_double(0.01, 0.02);

			pos.push(tempPos);
			vel.push(Vector2D(0.0, 0.0));
			double temp_vorticity = 0.00;

			data->vorticities().push(temp_vorticity);
			_newVelocities.push(Vector2D(0.0, 0.0));
			_newPositions.push(Vector2D());

			data->numberOfParticles()++;
		}

	}
	step++;
}




//计算index点处的速度值
//这里的index是流体粒子的index
Vector2D VortexParticleSolver::computeSingleVelocityFromPanels(int index) {
	auto& pos = _vortexParticleData->positions();
	auto panel = _vortexParticleData->panelSet;
	auto panelSize = panel->size();
	auto& gama1 = _vortexParticleData->strength;

	Vector2D result;

	for (int i = 0; i < panelSize; ++i) {

		auto start = panel->lookAt(i).start;
		auto end = panel->lookAt(i).end;
		auto normal = panel->lookAt(i).normal;

		Vector2D temp2;

		double beta = 0.0;
		auto vec_r1 = start - pos[index];
		auto vec_r2 = end - pos[index];
		auto r1 = vec_r1.getLength();
		auto r2 = vec_r2.getLength();
		auto temp1 = vec_r1.dot(vec_r2) / (r1 * r2);
		beta = acos(temp1);
		if (isnan(beta)) {
			beta = kPiD;
		}

		temp2.x = beta / (2.0 * kPiD);
		temp2.y = log((r2 + fv_eps) / (r1 + fv_eps)) / (2.0 * kPiD);

		double tempgamma = 0.0;
		if (gama1.size() > 0)
			tempgamma = gama1[i];
		temp2 = temp2 * tempgamma;
		temp2 = vel_to_world(temp2, start, normal);
		result += temp2;
	}


	return result;
}



//这里的index是panel的index
Vector2D VortexParticleSolver::computeUnitVelocityFromPanels(int index, const Vector2D& midPoint) {
	auto panel = _vortexParticleData->panelSet;
	auto start = panel->lookAt(index).start;
	auto end = panel->lookAt(index).end;
	auto normal = panel->lookAt(index).normal;

	Vector2D result;

	//计算beta值.
	double beta = 0.0;
	auto vec_r1 = start - midPoint;
	auto vec_r2 = end - midPoint;
	auto r1 = vec_r1.getLength();
	auto r2 = vec_r2.getLength();
	auto temp1 = vec_r1.dot(vec_r2) / (r1 * r2);
	beta = acos(temp1);
	if (isnan(beta)) {
		beta = kPiD;
	}
	//eq(16) eq(17)
	result.x = beta / (2.0 * kPiD);
	result.y = log((r2 + fv_eps) / (r1 + fv_eps)) / (2.0 * kPiD);
	result = vel_to_world(result, start, normal);
	return result;
}






//只执行一次
void VortexParticleSolver::computeBoundaryMatrix() {
	auto panels = _vortexParticleData->panelSet;
	auto panelSize = panels->size();
	Eigen::MatrixXd& A = _vortexParticleData->A;
	Eigen::MatrixXd& B = _vortexParticleData->B;

	auto sizex = panelSize + 1;
	auto sizey = panelSize;
	A.resize(sizex, sizey);
	B.resize(sizex, sizey);

	for (int j = 0; j < panelSize; ++j) {
		auto normal = panels->lookAt(j).normal;
		for (int i = 0; i < panelSize; ++i) {
			auto mid_j = panels->midPoint(j);
			auto u_ji = computeUnitVelocityFromPanels(i, mid_j).dot(normal);
			A(j, i) = -u_ji;
		}
	}

	for (int i = 0; i < panelSize; ++i) {
		A(panelSize, i) = 1.0;
	}

	for (int j = 0; j < panelSize; ++j) {
		auto normal = panels->lookAt(j).normal;
		Vector2D tengent(normal.y, -normal.x);
		for (int i = 0; i < panelSize; ++i) {
			auto mid_j = panels->midPoint(j);
			auto u_ji = computeUnitVelocityFromPanels(i, mid_j).dot(tengent);
			B(j, i) = -u_ji;
		}
	}

	for (int i = 0; i < panelSize; ++i) {
		B(panelSize, i) = 1.0;
	}
}





//在这里求解vortex sheet strength
void VortexParticleSolver::vortexSheetSolve(double timeIntervalInSeconds) {
	auto data = _vortexParticleData;
	auto n = data->numberOfParticles();
	auto panels = data->panelSet;
	auto panelSize = panels->size();
	auto& pos = data->positions();
	Eigen::MatrixXd& A = _vortexParticleData->A;
	Eigen::VectorXd& x = _vortexParticleData->strength;

	//组装b
	Eigen::VectorXd b(panelSize + 1);
	for (int i = 0; i < panelSize; ++i) {
		auto normal = panels->lookAt(i).normal;
		auto pos = panels->midPoint(i);
		auto vec = vs_vec;

		Vector2D temp;
		for (int j = 0; j < n; ++j) {
			temp += computeUSingle(pos, j);
		}
		b[i] = (vec + temp).dot(normal);
	}

	b[panelSize] = 0;

	x = A.colPivHouseholderQr().solve(b);

}


//在这里求解vortex sheet strength
void VortexParticleSolver::slipVortexSheetSolve(double timeIntervalInSeconds) {
	auto data = _vortexParticleData;
	auto n = data->numberOfParticles();
	auto panels = data->panelSet;
	auto panelSize = panels->size();
	auto& pos = data->positions();
	auto& vorticity = data->vorticities();

	Eigen::MatrixXd& B = _vortexParticleData->B;
	Eigen::VectorXd& x1 = _vortexParticleData->slipStrength;



	Eigen::VectorXd b1(panelSize + 1);
	for (int i = 0; i < panelSize; ++i) {
		auto normal = panels->lookAt(i).normal;
		Vector2D tengent(normal.y, -normal.x);

		auto pos = panels->midPoint(i);
		auto vec = vs_vec;

		Vector2D temp;
		for (int j = 0; j < n; ++j) {
			temp += computeUSingle(pos, j);
		}
		b1[i] = (vec + temp).dot(tengent);
	}

	b1[panelSize] = 0;

	x1 = B.colPivHouseholderQr().solve(b1);

	auto line = panels->lookAt(0).end - panels->lookAt(0).start;
	double panelLength = line.getLength();

	for (int i = 0; i < panelSize; ++i) {
		vorticity[n - i - 1] = (vs_tau * x1[panelSize - i - 1] * timeIntervalInSeconds * panelLength);
	}
}