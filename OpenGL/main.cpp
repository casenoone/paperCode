#include <iostream>

#include <fstream>
#include <cmath>
#include <string>
#include <vector>

#include "../titmouse2d/src/Geometry/Sphere2.h"
#include "vp2021.h"
#include "../titmouse2d/src/Geometry/RegularPolygon.h"
#include <GL/glut.h>

int sim_step = 0;

const float SCREEN_SIZE = 600;
const float DRAW_SIZE = SCREEN_SIZE / 200 * 10;


static void key(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 27: //按ESC
	case 'q': //按q键程序都会退出
		exit(0);
		break;
	}

	glutPostRedisplay();  //产生重绘回调
}

void drawPoint(double x, double y)
{
	//在后缓存绘制图形，就一个点
	glPointSize(3.05f);//缺省是1
	glBegin(GL_POINTS);
	glColor3f(1, 128.0 / 255, 51.0 / 255);
	glVertex3f((x - 1) * DRAW_SIZE, (y - 1) * DRAW_SIZE, 0);
	glEnd();
}


void drawPoint(double x, double y, double color_x, double color_y, double color_z)
{
	//在后缓存绘制图形，就一个点
	glPointSize(4.05f);//缺省是1
	glBegin(GL_POINTS);
	glColor3f(color_x / 255, color_y / 255, color_z / 255);
	glVertex3f((x - 1) * DRAW_SIZE, (y - 1) * DRAW_SIZE, 0);
	glEnd();
}


void drawLine(double x1, double y1, double x2, double y2) {

	glLineWidth(1);//设置线段宽度
	glBegin(GL_LINES);
	glColor3f(1.0, 0.0, 0.0);
	glVertex2f((x1 - 1) * DRAW_SIZE, (y1 - 1) * DRAW_SIZE); //定点坐标范围
	glVertex2f((x2 - 1) * DRAW_SIZE, (y2 - 1) * DRAW_SIZE);
	glEnd();
	glFlush();
}



auto vpSolver = std::make_shared<Vp2021Solver>();
RegularPolygonPtr obj1 = std::make_shared<RegularPolygon>(21, Vector2D(0.1, 1), 0.06);

double dt = 0.01;
int bubble_num = 0;

static void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glLoadIdentity();
	gluLookAt(0, 0, 100, 0, 0, 0, 0, 1, 0);
	vpSolver->onAdvanceTimeStep(dt);

	obj1->updatePosition(dt);

	//可视化涡粒子
	auto pos = vpSolver->vp2021Data()->vortexPosition;
	auto pos_n = pos.dataSize();

	/*for (int i = 0; i < pos_n; ++i) {
		drawPoint(pos[i].x, pos[i].y);
	}*/

	//可视化tracer粒子
	auto tracer_pos = vpSolver->vp2021Data()->tracePosition;
	int tracer_n = tracer_pos.dataSize();

	for (int i = 0; i < tracer_n; ++i) {
		if ((tracer_pos[i] - obj1->center()).getLength() > obj1->r())
			drawPoint(tracer_pos[i].x, tracer_pos[i].y);
	}

	//可视化移动边界
	for (auto i = obj1->_data.begin(); i != obj1->_data.end(); ++i) {
		auto start = i->start;
		auto end = i->end;
		drawLine(start.x, start.y, end.x, end.y);
	}



	glutSwapBuffers();
}

static void idle(void) {
	glutPostRedisplay();
}

static void resize(int width, int height)
{
	const float ar = (float)width / (float)height;
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	double ratio = 0.1;
	glOrtho(-width * 0.5 * ratio, width * 0.5 * ratio, -height * 0.5 * ratio, height * 0.5 * ratio, 2.0, 100.0); //这里使用正视投影
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}


int main(int argc, char** argv)
{

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(SCREEN_SIZE, SCREEN_SIZE);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("titmouse2d");

	glClearColor(6 / 255.0, 133 / 255.0, 135 / 255.0, 1);
	glShadeModel(GL_FLAT);

	vpSolver->setMovingBoudnary(obj1);
	vpSolver->emitTracerParticles();
	obj1->velocity = Vector2D(3, 0.0);





	glutKeyboardFunc(key);       //键盘按下去时
	glutIdleFunc(idle);          //空闲时
	glutReshapeFunc(resize);     //改变窗口大小时
	glutDisplayFunc(display);    //绘制窗口显示时

	glutMainLoop();





	return 0;
}





