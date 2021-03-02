#include "stdafx.h"

// standard
#include <assert.h>
#include <math.h>

// glut
#include <GL/glut.h>

//================================
// global variables
//================================
// screen size
int g_screenWidth  = 0;
int g_screenHeight = 0;

static int pointindex = 0;

static GLfloat t = 0; //This variable is time variable
static GLfloat M[16] = { 0 }; //A matrix that used to store data

void Normalize(GLfloat NM[7]) { //Normalization is used to offset the error that may be accumelate.
	GLfloat vectorlen = sqrt(NM[0] * NM[0] + NM[1] * NM[1] + NM[2] * NM[2] + NM[3] * NM[3]);
	if (vectorlen != 0) {//to prevent dividing by 0
		NM[0] /= vectorlen;
		NM[1] /= vectorlen;
		NM[2] /= vectorlen;
		NM[3] /= vectorlen;
	}
}
static GLfloat quaternionpoint[7][7] = { { 1, 0, 0, 0, -7, 7, -5 },   //The quaternion and their world cordinate.
										  { 0, 1, 0, 0, -5, 5, -10 },  
										  { 0, 0, 1, 0, -3, 3, -15 },  
										  { 0, 0, 0, 1, -1, 1, -20 }, 
										  { 0, 0, 1, 0, 3, -3, -20 },  
										  { 0, 1, 0, 0, 5, -5, -15 },
										  {1, 0, 0, 0, 7, -7, -10} };

static GLfloat eulerpoint[7][6] = { { 90, 0, 45, -7, 7, -5 },		//The euler point angle and their world cordinate.
									 { 70, 20, 65, -5, 5, -10 },	
									 { 50, 40, 85, -3, 3, -15 },	
									 { 30, 60, 105, -1, 1, -20 },	
									 { 50, 40, 85, 3, -3, -20 },		
									 { 70, 20, 65, 5, -5, -15 },		
									 { 90, 0, 45, 7, -7, -10 } };	
static int a = 0.5f;
static GLfloat CRspline[16] = { -a,2 - a, a - 2, a,
							2 * a,a - 3, 3 - 2 * a, -a,
							-a, 0.0f, a, 0.0f,
							0.0f, 1.0f, 0.0f, 0.0f };
static int v = 6.0f;
static GLfloat Bspline[16] = { -1.0f / v, 3.0f / v, -3.0f / v, 1.0f / v,
							3.0f / v, -6.0f / v, 3.0f / v, 0.0f / v,
						   -3.0f / v, 0.0f / v, 3.0f / v, 0.0f / v,
							1.0f / v, 4.0f / v ,1.0f / v, 0.0f / v };//CRSpline and BSpline that provided on ppt
GLfloat blendfunction(GLfloat T[4], GLfloat M[16], GLfloat G[4]) {
	GLfloat B[4] = { 0 };
	//Firstly we multiply T and spline matrix. This is dot product of matrix, so row times column
	B[0] = T[0] * M[0] + T[1] * M[4] + T[2] * M[8] + T[3] * M[12];
	B[1] = T[0] * M[1] + T[1] * M[5] + T[2] * M[9] + T[3] * M[13];
	B[2] = T[0] * M[2] + T[1] * M[6] + T[2] * M[10] + T[3] * M[14];
	B[3] = T[0] * M[3] + T[1] * M[7] + T[2] * M[11] + T[3] * M[15];
	//Then multiply with point
	GLfloat Q = B[0] * G[0] + B[1] * G[1] + B[2] * G[2] + B[3] * G[3];
	//Q = T*M*G
	return Q;
}
void QuaternionR(GLfloat QM[7], GLfloat R[16]) {
	GLfloat w = QM[0];
	GLfloat x = QM[1];
	GLfloat y = QM[2];
	GLfloat z = QM[3];
	R[0] = 1.0f - 2.0f * y * y - 2.0f * z * z; //The quaternion rotation matrix, I write in column order.
	R[1] = 2.0f * x * y + 2.0f * w * z;
	R[2] = 2.0f * x * z - 2.0f * w * y;
	R[3] = 0.0f;
	R[4] = 2.0f * x * y - 2.0f * w * z;
	R[5] = 1.0f - 2.0f * x * x - 2.0f * z * z;
	R[6] = 2.0f * y * z + 2.0f * w * x;
	R[7] = 0.0f;
	R[8] = 2.0f * x * z + 2.0f * w * y;
	R[9] = 2.0f * y * z - 2.0f * w * x;
	R[10] = 1.0f - 2.0f * x * x - 2.0f * y * y;
	R[11] = 0.0f;
	R[12] = QM[4];
	R[13] = QM[5];
	R[14] = QM[6];
	R[15] = 1.0f;
}



void qinter(GLfloat pquaternion[7][7], GLfloat Splinematrix[16]) {//quaternion interpolation, each loop read 4 points and are caculated by splinematrix and t variable.
	GLfloat TM[4] = {t*t*t, t*t, t, 1};
	GLfloat Mtemp[7];
	int i = 0;
	while (i<7){
		GLfloat qmatrix[4] = { pquaternion[pointindex][i],pquaternion[pointindex + 1][i],pquaternion[pointindex + 2][i],pquaternion[pointindex + 3][i]};
		Mtemp[i] = blendfunction(TM, Splinematrix, qmatrix);
		i++;
	}

	Normalize(Mtemp);
	QuaternionR(Mtemp, M);
}
void euinter(GLfloat peuler[7][6], GLfloat Splinematrix[16]) { //similar with quaternion interpolation but add a euler to quaternion step.
	GLfloat TM[4] = { t * t * t, t * t, t, 1 };
	GLfloat Mtemp[7];
	int i = 0;
	while (i<6){
		GLfloat eumatrix[4] = { peuler[pointindex][i],peuler[pointindex + 1][i],peuler[pointindex + 2][i],peuler[pointindex + 3][i]};
		Mtemp[i] = blendfunction(TM, Splinematrix, eumatrix);
		i++;
	}
	//The following is translate Euler to quaternion
	GLfloat a = Mtemp[0] / 2.0f;
	GLfloat b = Mtemp[1] / 2.0f;
	GLfloat c = Mtemp[2] / 2.0f;
	Mtemp[6] = Mtemp[5];
	Mtemp[5] = Mtemp[4];
	Mtemp[4] = Mtemp[3];
	Mtemp[0] = cos(c) * cos(b) * cos(c) + sin(c) * sin(b) * sin(a);
	Mtemp[1] = sin(c) * cos(b) * cos(c) - cos(c) * sin(b) * sin(a);
	Mtemp[2] = cos(c) * sin(b) * cos(c) + sin(c) * cos(b) * sin(a);
	Mtemp[3] = cos(c) * cos(b) * sin(c) - sin(c) * sin(b) * cos(a);
	Normalize(Mtemp);
	QuaternionR(Mtemp, M);
}

//================================
// init
//================================
void init( void ) {
	// init something before main loop...
}
void animation() {//Here is the 4 mode of animation
	//qinter(quaternionpoint, CRspline);
	//qinter(quaternionpoint,Bspline);
	//euinter(eulerpoint,CRspline);
	euinter(eulerpoint,Bspline);

	glLoadMatrixf(M);
	glutSolidTeapot(1.0);
}
//================================
// render
//================================
void render( void ) {
	// clear buffer
	glClearColor (0.0, 0.0, 0.0, 0.0);
	glClearDepth (1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	
	
	// render state
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);

	// enable lighting
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	// light source attributes
	GLfloat LightAmbient[]	= { 0.4f, 0.4f, 0.4f, 1.0f };
	GLfloat LightDiffuse[]	= { 0.3f, 0.3f, 0.3f, 1.0f };
	GLfloat LightSpecular[]	= { 0.4f, 0.4f, 0.4f, 1.0f };
	GLfloat LightPosition[] = { 5.0f, 5.0f, 5.0f, 1.0f };

	glLightfv(GL_LIGHT0, GL_AMBIENT , LightAmbient );
	glLightfv(GL_LIGHT0, GL_DIFFUSE , LightDiffuse );
	glLightfv(GL_LIGHT0, GL_SPECULAR, LightSpecular);
	glLightfv(GL_LIGHT0, GL_POSITION, LightPosition);

	// surface material attributes
	GLfloat material_Ka[]	= { 0.11f, 0.06f, 0.11f, 1.0f };
	GLfloat material_Kd[]	= { 0.43f, 0.47f, 0.54f, 1.0f };
	GLfloat material_Ks[]	= { 0.33f, 0.33f, 0.52f, 1.0f };
	GLfloat material_Ke[]	= { 0.1f , 0.0f , 0.1f , 1.0f };
	GLfloat material_Se		= 10;

	glMaterialfv(GL_FRONT, GL_AMBIENT	, material_Ka);
	glMaterialfv(GL_FRONT, GL_DIFFUSE	, material_Kd);
	glMaterialfv(GL_FRONT, GL_SPECULAR	, material_Ks);
	glMaterialfv(GL_FRONT, GL_EMISSION	, material_Ke);
	glMaterialf (GL_FRONT, GL_SHININESS	, material_Se);

	// modelview matrix
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
	

	// render objects
	animation();

	// disable lighting
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHTING);

	// swap back and front buffers
	glutSwapBuffers();
}

//================================
// keyboard input
//================================
void keyboard( unsigned char key, int x, int y ) {
}

//================================
// reshape : update viewport and projection matrix when the window is resized
//================================
void reshape( int w, int h ) {
	// screen size
	g_screenWidth  = w;
	g_screenHeight = h;	
	
	// viewport
	glViewport( 0, 0, (GLsizei)w, (GLsizei)h );

	// projection matrix
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	gluPerspective(45.0, (GLfloat)w/(GLfloat)h, 1.0, 2000.0);
}


//================================
// timer : triggered every 16ms ( about 60 frames per second )
//================================
void timer( int value ) {	
	// render
	glutPostRedisplay();
	int pointnumber = 7;
	t = t + 0.02;
	if (t >= 1) { //Since t from 0 to 1, use if-loop to reset time.
		t = 0;
		if (pointindex < pointnumber - 4) {
			pointindex++;
		}
		else {
			pointindex = 0;
		}
	}
	
	glutTimerFunc( 16, timer, 0 );
}

//================================
// main
//================================
int main( int argc, char** argv ) {
	// create opengL window
	glutInit( &argc, argv );
	glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGB |GLUT_DEPTH );
	glutInitWindowSize( 600, 600 ); 
	glutInitWindowPosition( 100, 100 );
	glutCreateWindow( "lab1" );

	// init
	init();
	
	// set callback functions
	glutDisplayFunc( render );
	glutReshapeFunc( reshape );
	glutKeyboardFunc( keyboard );
	glutTimerFunc( 16, timer, 0 );
	
	// main loop
	glutMainLoop();

	return 0;
}