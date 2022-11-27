/*==================================================================================
* COSC363 Computer Graphics (2022)
* A basic ray tracer with exntended features
*	Anti-Aliasing
*	Multiple Lights
*	Transparenct Objects
*	Fog
*===================================================================================
*/
#include <iostream>
#include <cmath>
#include <vector>
#include <glm/glm.hpp>
#include "TextureBMP.h"
#include "Sphere.h"
#include "Plane.h"
#include "SceneObject.h"
#include "Ray.h"
#include <GL/freeglut.h>
#include <algorithm> 
#include <math.h>
using namespace std;

#define ambient 0.2f;

const bool antiAliasing = true;
const int ALIASING_DEPTH = 1;

const float floorHeight = -15.0f;		//yp of the plane
const float floorWidth = 800;			//width of the plan along x axis)

const bool multipleLights = true;
const bool lighting = true;
const glm::vec3 lightPos = glm::vec3(0, floorHeight + 29, -55.0);		// just below ceiling
const glm::vec3 lightPos2 = glm::vec3(-10.0, floorHeight + 6, -30.0);	// back wall

const bool fog = false;
const glm::vec3 fogColor = glm::vec3(0.8);

TextureBMP texture;
const int wsize = 1260;		// window = size x size pixels
const float EDIST = 32;		// the distance of the image plane from the camera/origin
const int NUMDIV = wsize*2;	// the number of cells (subdivisions of the image plane) along x and y directions.
const int MAX_STEPS = 6;	// the number of levels (depth) of recursion

const float XMIN = -20.0;	// the boundary values of the image plane defined such that the view axis passes through its centre.
const float XMAX = 20.0;	//
const float YMIN = -10.0;	//
const float YMAX = 10.0;	//

vector<SceneObject*> sceneObjects;


//---The most important function in a ray tracer! ---------------------------------- 
//   Computes the colour value obtained by tracing a ray and finding its 
//     closest point of intersection with objects in the scene.
//----------------------------------------------------------------------------------
glm::vec3 trace(Ray ray, int step)
{
	glm::vec3 backgroundCol(0.0);						
		
	glm::vec3 color(0);
	SceneObject* obj;

    ray.closestPt(sceneObjects);	// Compare the ray with all objects in the scene to find the closest object
	
	// no intersection
	if (fog && ray.index == -1) {
		return backgroundCol + fogColor;		
	} else if (ray.index == -1)
		return backgroundCol;

	obj = sceneObjects[ray.index];					//object on which the closest point of intersection is found
	
	// ray hit the floor plane (index 0 in scene objects) 
	if (ray.index == 0)
	{	// chequered pattern
		glm::vec3 chequerColourPrimary = glm::vec3(1);
		glm::vec3 chequerColourSecondary = glm::vec3(0);
		bool checkeredW = false;
		bool checkeredH = false;
		int stripeWidth = 8;
		int stripeHeight = 4;
		
		int iz = (ray.hit.z) / stripeWidth;
		int k = iz % 2; //2 colors
		checkeredW = k == 0;
		
		int ix = (ray.hit.x-(1000*stripeHeight)) / stripeHeight;
		int j = ix % 2; //2 colors	
		checkeredH = j == 0;

		color = chequerColourPrimary;
		if (checkeredW && !checkeredH || checkeredH && !checkeredW) color = chequerColourSecondary;
		obj->setColor(color);

		// texture mapping 
		if (false) {		
			float x1 = -15, x2 = 5, z1 = -60, z2 = -90;// x1/z1 is bottom left of made up plane, x2/z2 top right
			float texcoords = (ray.hit.x - x1) / (x2 - x1);
			float texcoordt = (ray.hit.z - z1) / (z2 - z1);
			if (texcoords > 0 && texcoords < 1 &&
				texcoordt > 0 && texcoordt < 1)
			{
				color = texture.getColorAt(texcoords, texcoordt);
				obj->setColor(color);
			}
		}
	}

	glm::vec3 ambientColor = obj->getColor()* ambient;// ambient colour
	
	if (lighting) {
		
		bool shadow1 = false;
		bool shadow2 = false;
		SceneObject* obj1;
		SceneObject* obj2;
		
		// calculations for light 1
		glm::vec3 color1 = obj->lighting(lightPos, -ray.dir, ray.hit, multipleLights);// Object's colour after applying effects of lighting
		glm::vec3 lightVec = lightPos - ray.hit;
		Ray shadowRay(ray.hit, lightVec);
		shadowRay.closestPt(sceneObjects);
		if (shadowRay.index > -1) { // vector from the obj toward the light source intersects with another object (could potentially be between the obj and the light)		
			if (shadowRay.dist < glm::length(lightVec)) { // if the blocking obj is closer than the light it is blocking the light and therefore a shadow is cast
				obj1 = sceneObjects[shadowRay.index];			
				if (obj1->isTransparent()) {
					color1 *= 0.75*(obj1->getTransparencyCoeff());// trace(Ray(shadowRay.hit, shadowRay.hit*2.0f), step + 1);
				}
				else {
					shadow1 = true;
				}		
			}
		}
		
		// calculations for light 2
		glm::vec3 color2 = glm::vec3(0);
		if (multipleLights) {
			color2 = obj->lighting(lightPos2, -ray.dir, ray.hit, multipleLights);//obj->getColor();						//Object's colour after applying effects of ligth
			glm::vec3 lightVec2 = lightPos2 - ray.hit;
			Ray shadowRay2(ray.hit, lightVec2);
			shadowRay2.closestPt(sceneObjects);
			if (shadowRay2.index > -1) {//vector from the obj toward the light source intersects with another object (could potentially be between the obj and the light)
				if (shadowRay2.dist < glm::length(lightVec2)) { //if the blocking obj is closer than the light it is blocking the light and therefore a shadow is cast
					obj2 = sceneObjects[shadowRay2.index];
					if (obj2->isTransparent()) {
						color2 *= 0.75 *(obj2->getTransparencyCoeff());//trace(Ray(shadowRay.hit, shadowRay.hit*2.0f), step + 1);
					}
					else {
						shadow2 = true;
					}
				}
			}	
			const float x = 1.0f;
			if (shadow1 && shadow2) {
				color = ambientColor;
			}
			else if (shadow1) {
				color = ambientColor + color2 * x;
			}
			else if (shadow2) {
				color = ambientColor + color1 * x;
			}
			else {
				color = ambientColor + color1 + color2;
			}
		}
		else {
			if (shadow1) {
				color = ambientColor;
			}
			else {
				color = color1;
			}
		}		
	}
	else {
		color = ambientColor;
	}
	
	

	if (obj->isReflective() && step < MAX_STEPS)
	{
		float rho = obj->getReflectionCoeff();
		glm::vec3 normalVec = obj->normal(ray.hit);
		glm::vec3 reflectedDir = glm::reflect(ray.dir, normalVec);
		Ray reflectedRay(ray.hit, reflectedDir);
		glm::vec3 reflectedColor = trace(reflectedRay, step + 1);
		color = color + (rho * reflectedColor);
	}

	if (obj->isTransparent() && step < MAX_STEPS) {
		float rho = obj->getTransparencyCoeff();
		Ray transparentRay(ray.hit, ray.hit * 2.0f);//ray on surface of transparent object, continuing forward (potentially through itself)
		glm::vec3 transparentColor = trace(transparentRay, step + 1);
		//color = ambientColor + (rho * transparentColor);
		if (step != 0) {
			//color -= ambientColor;
		}
		color *= 1 - rho;
		color += (rho * transparentColor);
	}

	

	if (fog) {
		float z1 = -60, z2 = -400, t=0, x=180;
		if (true) { // exponential fog
			const float exponent = min(ray.hit.z - z1, 0.0f);
			t = 1 - exp(exponent / x);
		}
		else {		// linear fog
			t = max(0.0f, min((ray.hit.z - z1) / (z2 - z1), 1.0f));//linear fog 
		}
		return ((1.0f - t) * color + (t * fogColor));
	}
	return color;
}

//---------------------------------------------------------------------------------------
//given a cell (bottom left point (xp,yp), width (cellX), and height(cellY)) returns the calculated colour
glm::vec3 cellColour(float xp, float yp, float cellX, float cellY, glm::vec3 eye) {
	glm::vec3 dir(xp + cellX, yp + cellY, -EDIST);
	Ray ray = Ray(eye, dir);
	glm::vec3 col = trace(ray, 1);
	return col;
}

//---------------------------------------------------------------------------------------
//	given a cell (bottom left point (xp,yp), width (cellX), and height(cellY)) 
//	if depth == 1:
//		calculates average col of cell from 4 sub-cells
//	if depth == 2:
//		calculates average col of cell from 16 sub-cells (takes average of the 4 sub-cells from depth 1)
//---------------------------------------------------------------------------------------
glm::vec3 antiAliasingColour(float xp, float yp, float cellX, float cellY, glm::vec3 eye, float depth=ALIASING_DEPTH) {
	if (depth == 0) {//return colour of the cell
		return cellColour(xp, yp, cellX * 0.5, cellY * 0.5, eye);
	}
	glm::vec3 LL = antiAliasingColour(xp, yp, cellX * 0.5, cellY * 0.5, eye, depth-1);
	glm::vec3 LR = antiAliasingColour(xp + cellX * 0.5, yp, cellX * 0.5, cellY * 0.5, eye, depth - 1);
	glm::vec3 TL = antiAliasingColour(xp, yp + cellY * 0.5, cellX * 0.5, cellY * 0.5, eye, depth - 1);
	glm::vec3 TR = antiAliasingColour(xp + cellX * 0.5, yp + cellY * 0.5, cellX * 0.5, cellY * 0.5, eye, depth - 1);
	glm::vec3 averageColour = glm::vec3(LL + LR + TL + TR) / 4.0f;
	return averageColour;
}

//---The main display module -----------------------------------------------------------
// In a ray tracing application, it just displays the ray traced image by drawing
// each cell as a quad.
//---------------------------------------------------------------------------------------
void display()
{
	float xp, yp;  //grid point
	float cellX = (XMAX - XMIN) / NUMDIV;  //cell width
	float cellY = (YMAX - YMIN) / NUMDIV;  //cell height
	glm::vec3 eye(0., 0., 0.);

	glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

	glBegin(GL_QUADS);  //Each cell is a tiny quad.

	for(int i = 0; i < NUMDIV; i++)	//Scan every cell of the image plane
	{
		xp = XMIN + i * cellX;
		for(int j = 0; j < NUMDIV; j++)
		{
			yp = YMIN + j * cellY;
		
			glm::vec3 col;
			if (antiAliasing) {
				col = antiAliasingColour(xp, yp, cellX, cellY, eye);
			}
			else {
				glm::vec3 dir(xp + 0.5 * cellX, yp + 0.5 * cellY, -EDIST);		//direction of the primary ray
				Ray ray = Ray(eye, dir);
				col = trace(ray, 1); //Trace the primary ray and get the colour value
			}

			//draws the square pixel
			glColor3f(col.r, col.g, col.b);
			glVertex2f(xp, yp);				//Draw each cell with its color value
			glVertex2f(xp + cellX, yp);
			glVertex2f(xp + cellX, yp + cellY);
			glVertex2f(xp, yp + cellY);
        }
    }

    glEnd();
    glFlush();
}

//---
//constructs cube from planes and adds them to sceneObjects
void cube(glm::vec3 pos, bool transparency=false, bool reflect=false, float baseRadius = 10.0f, float height = 10.0f, glm::vec3 color = glm::vec3(0.2), bool specular = false) {
	const float reflectionCoeff = 0.8f;
	const float shininess = 2048;
	const float transparencyCoeff = 0.8f;

	//cube has 8 vertices/points (4 on floor plane, 4 on roof plane)
	const glm::vec3 BFL = glm::vec3(pos.x - baseRadius * 0.5, pos.y, pos.z + baseRadius * 0.5);//base front left
	const glm::vec3 BFR = glm::vec3(pos.x + baseRadius * 0.5, pos.y, pos.z + baseRadius * 0.5);//front right
	const glm::vec3 BBR = glm::vec3(pos.x + baseRadius * 0.5, pos.y, pos.z - baseRadius * 0.5);//back right
	const glm::vec3 BBL = glm::vec3(pos.x - baseRadius * 0.5, pos.y, pos.z - baseRadius * 0.5);//back left etc..

	const glm::vec3 TFL = glm::vec3(pos.x - baseRadius * 0.5, pos.y + height, pos.z + baseRadius * 0.5);//top front left
	const glm::vec3 TFR = glm::vec3(pos.x + baseRadius * 0.5, pos.y + height, pos.z + baseRadius * 0.5);//
	const glm::vec3 TBR = glm::vec3(pos.x + baseRadius * 0.5, pos.y + height, pos.z - baseRadius * 0.5);//
	const glm::vec3 TBL = glm::vec3(pos.x - baseRadius * 0.5, pos.y + height, pos.z - baseRadius * 0.5);//

	Plane* PlaneBase = new Plane(BFL, BBL, BBR, BFR);
	PlaneBase->setColor(color);
	if (reflect) PlaneBase->setReflectivity(reflect, reflectionCoeff);
	PlaneBase->setShininess(shininess);
	PlaneBase->setSpecularity(specular);
	if(transparency) PlaneBase->setTransparency(transparency, transparencyCoeff);
	sceneObjects.push_back(PlaneBase);

	Plane* PlaneFront = new Plane(BFL, BFR, TFR, TFL);
	PlaneFront->setColor(color);
	if (reflect) PlaneFront->setReflectivity(reflect, reflectionCoeff);
	PlaneFront->setShininess(shininess);
	PlaneFront->setSpecularity(specular);
	if (transparency) PlaneFront->setTransparency(transparency, transparencyCoeff);
	sceneObjects.push_back(PlaneFront);

	Plane* PlaneBack = new Plane(BBR, BBL, TBL, TBR);
	PlaneBack->setColor(color);
	if (reflect) PlaneBack->setReflectivity(reflect, reflectionCoeff);
	PlaneBack->setShininess(shininess);
	PlaneBack->setSpecularity(specular);
	if (transparency) PlaneBack->setTransparency(transparency, transparencyCoeff);
	sceneObjects.push_back(PlaneBack);

	Plane* PlaneLeft = new Plane(BBL, BFL, TFL, TBL);
	PlaneLeft->setColor(color);
	if (reflect) PlaneLeft->setReflectivity(reflect, reflectionCoeff);
	PlaneLeft->setShininess(shininess);
	PlaneLeft->setSpecularity(specular);
	if (transparency) PlaneLeft->setTransparency(transparency, transparencyCoeff);
	sceneObjects.push_back(PlaneLeft);

	Plane* PlaneRight = new Plane(BFR, BBR, TBR, TFR);
	PlaneRight->setColor(color);
	if (reflect) PlaneRight->setReflectivity(reflect, reflectionCoeff);
	PlaneRight->setShininess(shininess);
	PlaneRight->setSpecularity(specular);
	if (transparency) PlaneRight->setTransparency(transparency, transparencyCoeff);
	sceneObjects.push_back(PlaneRight);

	Plane* PlaneTop = new Plane(TBR, TBL, TFL, TFR);
	PlaneTop->setColor(color);
	if (reflect) PlaneTop->setReflectivity(reflect, reflectionCoeff);
	PlaneTop->setShininess(shininess);
	PlaneTop->setSpecularity(specular);
	if (transparency) PlaneTop->setTransparency(transparency, transparencyCoeff);
	sceneObjects.push_back(PlaneTop);
}

//---This function initializes the scene ------------------------------------------- 
//   Specifically, it creates scene objects (spheres, planes, cones, cylinders etc)
//     and add them to the list of scene objects.
//   It also initializes the OpenGL 2D orthographc projection matrix for drawing the
//     the ray traced image.
//----------------------------------------------------------------------------------
void initialize()
{
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(XMIN, XMAX, YMIN, YMAX);
	glClearColor(0, 0, 0, 1);

	texture = TextureBMP("Butterfly.bmp");
	
	// floor plane
	if (true) {
		Plane* plane = new Plane(glm::vec3(-floorWidth * 0.5, floorHeight, 1400), //Point A
		glm::vec3(floorWidth * 0.5, floorHeight, 1400), //Point B
		glm::vec3(floorWidth * 0.5, floorHeight, -1400), //Point C
		glm::vec3(-floorWidth * 0.5, floorHeight, -1400)); //Point D
		plane->setColor(glm::vec3(0.8, 0.8, 0));
		plane->setSpecularity(false);
		sceneObjects.push_back(plane);
	}
	
	// stack of cubes
	if (true) {
		const float radius = 7.0f;
		glm::vec3 stackPos = glm::vec3(9.5, floorHeight, -70);
		glm::vec3 c1 = glm::vec3(0, 0, -0);				// cube1
		cube(stackPos + c1, false, false, radius, radius, glm::vec3(1, 0, 0),true);
		glm::vec3 c2 = glm::vec3(0.25, radius, -0.25);	// cube2
		cube(stackPos + c2, false, false, radius, radius, glm::vec3(0, 1, 0));
		glm::vec3 c3 = glm::vec3(0.5, radius*2, -0.5);	// cube3
		cube(stackPos + c3, false, false, radius, radius, glm::vec3(0, 0, 1));
	}

	// walls
	if (true) {
		// left 
		glm::vec3 pos = glm::vec3(95, floorHeight, -30);
		cube(pos, false, false, 150, 150, glm::vec3(1, 0, 0), false);

		// right
		glm::vec3 pos2 = glm::vec3(-95, floorHeight, -30);
		cube(pos2, false, false, 150, 150, glm::vec3(0, 1, 0), false);

		// far
		glm::vec3 pos3 = glm::vec3(0, floorHeight, -150);
		cube(pos3, false, false, 100, 100, glm::vec3(1, 1, 1), false);

		// near / behind camera
		glm::vec3 pos5 = glm::vec3(0, floorHeight, 50);
		cube(pos5, false, false, 100, 100, glm::vec3(0, 0, 0), false);

		// roof
		glm::vec3 pos4 = glm::vec3(0, floorHeight+30, -30);
		cube(pos4, false, false, 150, 150, glm::vec3(0, 0, 1), false);
	}

	// spheres
	if (true) {
		const float sphereRadius = 7.0f;

		// reflective sphere	
		Sphere* reflectiveSphere = new Sphere(glm::vec3(-11, floorHeight + sphereRadius, -70.0), sphereRadius);
		reflectiveSphere->setColor(glm::vec3(0.0));
		reflectiveSphere->setShininess(1024);
		reflectiveSphere->setSpecularity(true);
		reflectiveSphere->setReflectivity(true, 0.7f);
		sceneObjects.push_back(reflectiveSphere);

		// sphere under light
		Sphere* sphere = new Sphere(glm::vec3(lightPos.x, floorHeight + sphereRadius / 1.5, lightPos.z), sphereRadius / 1.5);
		sphere->setColor(glm::vec3(1,1,1));
		sphere->setShininess(128);
		sphere->setSpecularity(true);
		sceneObjects.push_back(sphere);
	}

	// mirror
	if (false) {
		const float mirrorSize = 12.0f;
		const float mirrorZ = -99.0f;
		const float mirrorY = 10.0f;
		Plane* mirror = new Plane(glm::vec3(-mirrorSize * 0.5, floorHeight + mirrorY, mirrorZ),		//Point A
			glm::vec3(mirrorSize * 0.5, floorHeight + mirrorY, mirrorZ),							//Point B
			glm::vec3(mirrorSize * 0.5, floorHeight + mirrorSize + mirrorY, mirrorZ),				//Point C
			glm::vec3(-mirrorSize * 0.5, floorHeight + mirrorSize + mirrorY, mirrorZ));				//Point D
		mirror->setColor(glm::vec3(0.05, 0.05, 0.05));
		mirror->setShininess(256);
		mirror->setSpecularity(false);
		mirror->setReflectivity(true, 0.8);
		sceneObjects.push_back(mirror);
	}
}

int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB );
    glutInitWindowSize(wsize*2, wsize);
    glutInitWindowPosition(20, 20);
	glutCreateWindow("C++ OpenGL Raytracing");

    glutDisplayFunc(display);
    initialize();

    glutMainLoop();
    return 0;
}
