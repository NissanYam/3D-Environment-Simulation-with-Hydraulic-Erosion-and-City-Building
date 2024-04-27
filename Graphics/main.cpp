
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <time.h>
#include <math.h>
#include "glut.h"
#include <cmath>
using namespace std;

const double PI = 3.14;
const int GSZ = 150;
const int H = 600;
const int W = 600;
double angle = 0;

typedef struct {
	double x, y, z;
} POINT3;
typedef struct {
	int x, z;
}POINT3_INT;
typedef struct {
	double R, G, B;
}RGB;
POINT3 eye = {2,24,20 };
double sight_angle = PI;
POINT3 sight_dir = {sin(sight_angle),0,cos(sight_angle)}; // in plane X-Z
double speed = 0;
double angular_speed = 0;
double ground[GSZ][GSZ] = { 0 };
double startingGround[GSZ][GSZ] = { 0 };
double tmp[GSZ][GSZ];
double delta = 0.0004;
bool start_rain = false;
bool build_city = false;
double waterDepthFactor = 0.3;
bool HydraulicErosionAct = false;
int maxDistFronCenter = 10;
int minimumRiverDist = 1;
int maximumRiverDist = 5;
int bottomR = 1;
int max_height_building = 5;
int min_height_building = 2;
int manximum_number_of_rain_points = 1000;
vector<vector<bool>> visited(GSZ, vector<bool>(GSZ, false)); // Matrix to mark visited points
vector<vector<int>> buildingHeights(GSZ, vector<int>(GSZ, 0)); // All heights initialized to 0
vector<vector<bool>> build(GSZ, vector<bool>(GSZ, false)); // Matrix to mark points to build
vector<vector<RGB>> rgbVector(GSZ, vector<RGB>(GSZ)); //A vector holding the color of the roof at point i,j

// Texture definitions

const int TW = 512; // must be a power of 2
const int TH = 512;
unsigned char tx0[TH][TW][3]; // RGB

void UpdateGround();
void Smooth();
void createWorld();
void copyStartingWorld();
void menu(int choice);

// fill up the texture matrix
void SetupTexture(int tnum)
{
	int i, j;
	int rnd;

	if (tnum == 0) // bricks
	{
		for (i = 0; i < TH; i++)
			for (j = 0; j < TW; j++)
			{
				rnd = rand() % 25;
				if (i % (TH / 2) < 20 ||
					(i < TH / 2 && j % (TW / 2) < 20) ||
					(i > TH / 2 && j % (TW / 4) < 20) && (j < TW / 3 || j>2 * TW / 3) && j > 20) //  gray
				{
					tx0[i][j][0] = 180 + rnd;
					tx0[i][j][1] = 180 + rnd;
					tx0[i][j][2] = 180 + rnd;

				}
				else if (i < TH / 2 && j < TW / 2) // brick
				{
					tx0[i][j][0] = 150 + rnd;
					tx0[i][j][1] = 60 + rnd;
					tx0[i][j][2] = 0 + rnd;
				}
				else if (i < TH / 2 && j > TW / 2)
				{
					tx0[i][j][0] = 160 + rnd;
					tx0[i][j][1] = 70 + rnd;
					tx0[i][j][2] = 0 + rnd;
				}
				else if (j > TW / 4 && j < 3 * TW / 4)
				{
					tx0[i][j][0] = 130 + rnd;
					tx0[i][j][1] = 40 + rnd;
					tx0[i][j][2] = 0 + rnd;
				}
				else
				{
					tx0[i][j][0] = 160 + rnd;
					tx0[i][j][1] = 70 + rnd;
					tx0[i][j][2] = 0 + rnd;
				}
			}
	}
	else if (tnum == 1) // Road
	{
		for (i = 0; i < TH; i++)
			for (j = 0; j < TW; j++)
			{
				rnd = rand() % 25;
				tx0[i][j][0] = 110 + rnd;
				tx0[i][j][1] = 110 + rnd;
				tx0[i][j][2] = 110 + rnd;
			}
	}
}
void init()
{
	//                 R     G    B
	glClearColor(0.5,0.7,1,0);// color of window background

	glEnable(GL_DEPTH_TEST);

	int i, j;

	srand(time(0));
	createWorld();
	copyStartingWorld();

	// define textures
	SetupTexture(1); // road
	glBindTexture(GL_TEXTURE_2D, 1);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, TW, TH, 0, GL_RGB, GL_UNSIGNED_BYTE, tx0);
}
/*Reset building construction variables*/
void resetConstruction()
{
	for (int i = 0; i < GSZ; i++)
		for (int j = 0; j < GSZ; j++)
		{
			build[i][j] = false;
			buildingHeights[i][j] = 0;
		}
}
void createWorld()
{
	delta = 0.04;///// Set the initial delta value for terrain generation
	int i,j;
	// Initialize the ground heights to zero
	for (i = 0; i < GSZ; i++)
	{
		for (j = 0; j < GSZ; j++)
		{
			ground[i][j] = 0;
		}
	}
	for (i = 0; i < 3000; i++)
		UpdateGround();
	Smooth();
	for (i = 0; i < 1000; i++)
		UpdateGround();
	delta = 0.0004;//// Reset delta value for subsequent updates
}
void UpdateGround()
{
	if (rand() % 2 == 0)
		delta = -delta;
	int x1, y1, x2, y2;
	x1 = rand() % GSZ;
	y1 = rand() % GSZ;
	x2 = rand() % GSZ;
	y2 = rand() % GSZ;
	double a, b;
	if (x1 != x2)
	{
		a = (y2 - y1) / ((double)x2 - x1);
		b = y1 - a * x1;
		for(int i=0;i<GSZ;i++)
			for (int j = 0;j < GSZ;j++)
			{
				if (i < a * j + b) ground[i][j] += delta;
				else ground[i][j] -= delta;
			}
	}
}
void Smooth()
{

	for(int i=1;i<GSZ-1;i++)
		for (int j = 1;j < GSZ - 1;j++)
		{
			tmp[i][j] = (ground[i-1][j-1]+ ground[i-1 ][j]+ ground[i - 1][j + 1]+
				ground[i][j - 1] + ground[i ][j] + ground[i ][j + 1]+
				ground[i + 1][j - 1] + ground[i + 1][j] + ground[i + 1][j + 1]) / 9.0;
		}

	for (int i = 1;i < GSZ - 1;i++)
		for (int j = 1;j < GSZ - 1;j++)
			ground[i][j] = tmp[i][j];

}
void SetColor(double h)
{
	h = fabs(h) / 6;
	// sand
	if (h < 0.02)
		glColor3d(0.8, 0.7, 0.5);
	else	if(h<0.3)// grass
		glColor3d(0.4+0.8*h,0.6-0.6*h,0.2+ 0.2 * h);
	else if(h<0.9) // stones
		glColor3d(0.4 + 0.1 * h,  0.4+0.1*h, 0.2 + 0.1 * h);
	else // snow
		glColor3d(  h,  h, 1.1 * h);
}
/*
Makes a layer of water under the entire surface
*/
void DrawWaterLayer() {
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	int i, j;
	float baseBlue = 0.7;
	float baseAlpha = 0.7; // Constant base for transparency
	for (i = 1; i < GSZ; i++)
		for (j = 1; j < GSZ; j++)
		{
			if (startingGround[i][j] > 0)
			{
				glBegin(GL_POLYGON);

				glColor4f(0, 0.4, baseBlue + startingGround[i][j]/10, baseAlpha + startingGround[i][j] * 0.3);
				glVertex3d(j - GSZ / 2, startingGround[i][j] - waterDepthFactor * startingGround[i][j], i - GSZ / 2);

				glColor4f(0, 0.4, baseBlue + startingGround[i - 1][j]/10, baseAlpha + startingGround[i - 1][j] * 0.3);
				glVertex3d(j - GSZ / 2, startingGround[i - 1][j] - waterDepthFactor * startingGround[i - 1][j], i - 1 - GSZ / 2);

				glColor4f(0, 0.4, baseBlue + startingGround[i - 1][j - 1]/10, baseAlpha + startingGround[i - 1][j - 1] * 0.3);
				glVertex3d(j - 1 - GSZ / 2, startingGround[i - 1][j - 1] - waterDepthFactor * startingGround[i - 1][j - 1], i - 1 - GSZ / 2);

				glColor4f(0, 0.4, baseBlue + startingGround[i][j - 1]/10, baseAlpha + startingGround[i][j - 1] * 0.3);
				glVertex3d(j - 1 - GSZ / 2, startingGround[i][j - 1] - waterDepthFactor * startingGround[i][j - 1], i - GSZ / 2);

				glEnd();
			}
		}

	// water + transparency
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glColor4d(0, 0.4, baseBlue, baseAlpha);
	glBegin(GL_POLYGON);
	glVertex3d(-GSZ / 2, 0, -GSZ / 2);
	glVertex3d(-GSZ / 2, 0, GSZ / 2);
	glVertex3d(GSZ / 2, 0, GSZ / 2);
	glVertex3d(GSZ / 2, 0, -GSZ / 2);
	glEnd();

	glDisable(GL_BLEND);
}
void DrawFloor()
{
	int i,j;

	glColor3d(0, 0, 0.3);

	for(i=1;i<GSZ;i++)
		for (j = 1;j < GSZ;j++)
		{
			glBegin(GL_POLYGON);
			SetColor(ground[i][j]);
			glVertex3d(j-GSZ/2, ground[i][j], i-GSZ/2);
			SetColor(ground[i-1][j]);
			glVertex3d(j - GSZ / 2, ground[i - 1][j], i - 1 - GSZ / 2);
			SetColor(ground[i-1][j-1]);
			glVertex3d(j - 1 - GSZ / 2, ground[i - 1][j - 1], i - 1 - GSZ / 2);
			SetColor(ground[i][j-1]);
			glVertex3d(j - 1 - GSZ / 2, ground[i ][j - 1], i - GSZ / 2);
			glEnd();
		}
	DrawWaterLayer();
}
void copyStartingWorld()
{
	for (int i = 0; i < GSZ; i++)
		for (int j = 0; j < GSZ; j++)
			startingGround[i][j] = ground[i][j];
}
/*
Finds the lowest point around x,z
Returns true if we found it and false if not
*/
bool findLowNearPoint(int* x, int* z) {
	int x_best = *x, z_best = *z;
	double lowestHeight = ground[*x][*z];
	bool foundLower = false;

	for (int i = -1; i <= 1; ++i) {
		for (int j = -1; j <= 1; ++j) {
			int new_x = *x + i, new_z = *z + j;
			if (new_x >= 0 && new_x < GSZ && new_z >= 0 && new_z < GSZ) {
				double newHeight = ground[new_x][new_z];
				// Update to find strictly lower height
				if (newHeight < lowestHeight) {
					lowestHeight = newHeight;
					x_best = new_x;
					z_best = new_z;
					foundLower = true;
				}
			}
		}
	}

	if (foundLower) {
		*x = x_best;
		*z = z_best;
	}
	return foundLower;
}
/*
Finding the route through which the water drops from the rain pass
*/
POINT3_INT* getRainRoute(int x, int z, int* size) {
	*size = 1;
	POINT3_INT* route = (POINT3_INT*)malloc(sizeof(POINT3_INT) * (*size));
	if (!route) {
		fprintf(stderr, "Memory allocation failed\n");
		return NULL;
	}
	route[0].x = x;
	route[0].z = z;

	while (true) {
		if (!findLowNearPoint(&x, &z)) {
			break;  // No lower point found, end of route
		}

		// Check for cycles
		bool isCycle = false;
		for (int i = 0; i < *size; i++) {
			if (route[i].x == x && route[i].z == z) {
				isCycle = true;
				break;
			}
		}
		if (isCycle) break;  // End the route if a cycle is detected

		// Add new point to the route
		(*size)++;
		POINT3_INT* new_route = (POINT3_INT*)realloc(route, sizeof(POINT3_INT) * (*size));
		if (!new_route) {
			fprintf(stderr, "Memory reallocation failed\n");
			free(route);
			return NULL;
		}
		route = new_route;
		route[(*size) - 1].x = x;
		route[(*size) - 1].z = z;

		// Stop if the point is below the water line
		if (ground[route[(*size) - 1].x][route[(*size) - 1].z] < 0 ||
			ground[route[(*size) - 1].x][route[(*size) - 1].z] <= startingGround[x][z] - waterDepthFactor *startingGround[x][z]) {
			break;
		}
	}

	return route;
}
/*
Hydraulic fracturing process
You get random places where the rain process starts
You get a path where the raindrops pass and each point on the route lowers its height
*/
void HydraulicErosion() {
	int x, z, size, i;
	do {
		x = rand() % GSZ;
		z = rand() % GSZ;
	} while (ground[x][z] <= 0 || 
		ground[x][z] < startingGround[x][z] - waterDepthFactor * startingGround[x][z]);
	int number_of_rain_points = rand() % manximum_number_of_rain_points;
	for (int l = 0; l < number_of_rain_points; l++)
	{
		size = 1;
		POINT3_INT* route = getRainRoute(x, z, &size);
		if (route == NULL) {
			fprintf(stderr, "Failed to generate rain route.\n");
			return;
		}


		for (i = 0; i < size - 1; i++) {
			x = route[i].x;
			z = route[i].z;
			ground[x][z] -= delta;  // Apply erosion
		}
		// Deposit sediment at the last point in the route
		int last_point_x = route[size - 1].x;
		int last_point_z = route[size - 1].z;
		if (ground[last_point_x][last_point_z] > 0 
			||
			ground[last_point_x][last_point_z] > startingGround[last_point_x][last_point_z] - waterDepthFactor * startingGround[last_point_x][last_point_z])
		{
			ground[last_point_x][last_point_z] += delta;
		}
		else
		{
			ground[last_point_x][last_point_z] -= delta;
		}

		free(route);  // Free the dynamically allocated memory
	}
}
/*
Draws a wall of a house
*/
void DrawWall(int x_start, int z_start, int x_end, int z_end, double average_height)
{
	glBegin(GL_POLYGON);
	glVertex3d(z_start - GSZ / 2, ground[x_start][z_start], x_start - GSZ / 2 );
	glVertex3d(z_end - GSZ / 2, ground[x_end][z_end], x_end - GSZ / 2);
	glVertex3d(z_end - GSZ / 2, average_height, x_end - GSZ / 2);
	glVertex3d(z_start - GSZ / 2, average_height, x_start - GSZ / 2);
	glEnd();
}
/*
Draws windows on a wall of a house
average_height Maximum height points of the building
*/
void DrawWindows(int x_start, int z_start, int x_end, int z_end,double average_height,int num_of_floors, double offset)
{
	int num_of_windows = 3;
	double width_wall = bottomR * 2 + 1;
	double width_window = width_wall / (num_of_windows * 2 + 1);
	
	double y_bottom = (ground[x_start][z_start] + ground[x_end][z_end]) / 2;
	double y_top = average_height;
	double height_floor = (average_height - y_bottom) / (double)(num_of_floors * 2 + 1);
	for (int f = 0; f < num_of_floors * 2 + 1; f++)
	{
		if (f % 2 != 0)// draw floor
		{
			double y_start_floor = y_bottom + f * height_floor;
			double y_end_floor = y_start_floor + height_floor;

			for (int w = 0; w < num_of_windows *2 +1 -2 ; w++)
			{
				if (w % 2 != 0 && y_end_floor < average_height)
				{
					double x_start_window, x_end_window, z_start_window, z_end_window;
					if (x_start == x_end)
					{
						x_start_window = x_start + offset;
						x_end_window = x_start + offset;

						z_start_window = z_start + w * width_window;
						z_end_window = z_start_window + width_window;
					}
					else if (z_start == z_end)
					{
						z_start_window = z_start + offset;
						z_end_window = z_start + offset;

						x_start_window = x_start + w * width_window;
						x_end_window = x_start_window + width_window;
					}
					glBegin(GL_POLYGON);
					glVertex3d(z_start_window - GSZ / 2, y_start_floor, x_start_window - GSZ / 2);
					glVertex3d(z_end_window - GSZ / 2, y_start_floor, x_end_window - GSZ / 2);
					glVertex3d(z_end_window - GSZ / 2, y_end_floor, x_end_window - GSZ / 2);
					glVertex3d(z_start_window - GSZ / 2, y_end_floor, x_start_window - GSZ / 2);
					glEnd();
				}
			}
		}
	}
}
/*
Draws a roof of a house
average_height Maximum height points of the building
*/
void DrawRoof(int x, int z, double average_height, double roofHeight)
{
	glBegin(GL_POLYGON);
	glVertex3d(z - bottomR - GSZ / 2, average_height, x - bottomR - GSZ / 2);
	glVertex3d(z + bottomR - GSZ / 2, average_height, x - bottomR - GSZ / 2);
	glVertex3d(z - GSZ / 2, average_height + roofHeight, x - GSZ / 2);
	glEnd();

	glBegin(GL_POLYGON);
	glVertex3d(z - bottomR - GSZ / 2, average_height, x + bottomR - GSZ / 2 );
	glVertex3d(z + bottomR - GSZ / 2, average_height, x + bottomR - GSZ / 2);
	glVertex3d(z - GSZ / 2, average_height + roofHeight, x - GSZ / 2 );
	glEnd();

	glBegin(GL_POLYGON);
	glVertex3d(z - bottomR - GSZ / 2, average_height, x - bottomR - GSZ / 2 );
	glVertex3d(z - bottomR - GSZ / 2, average_height, x + bottomR - GSZ / 2 );
	glVertex3d(z - GSZ / 2, average_height + roofHeight, x - GSZ / 2 );
	glEnd();

	glBegin(GL_POLYGON);
	glVertex3d(z + bottomR - GSZ / 2, average_height, x - bottomR - GSZ / 2 );
	glVertex3d(z + bottomR - GSZ / 2, average_height, x + bottomR - GSZ / 2 );
	glVertex3d(z - GSZ / 2, average_height + roofHeight, x - GSZ / 2 );
	glEnd();

}
/*
Draws a house at point x,z
average_height Maximum height points of the building
*/
void DrawBuilding(int x, int z)
{
	int floors = buildingHeights[x][z];
	double average_height = ((ground[x - bottomR][z - bottomR]
		+ ground[x - bottomR][z + bottomR]
		+ ground[x + bottomR][z + bottomR]
		+ ground[x + bottomR][z - bottomR]) / 4) + floors;
	RGB rgb = rgbVector[x][z];
	glColor3d((double)255 / 255, (double)191 / 255, (double)115 / 255);
	DrawWall(x - bottomR, z - bottomR, x - bottomR, z + bottomR, average_height);
	glColor3d((double)0 / 255, (double)0 / 255, (double)0 / 255);
	DrawWindows(x - bottomR, z - bottomR, x - bottomR, z + bottomR, average_height, floors,-0.01);


	glColor3d((double)255 / 255, (double)191 / 255, (double)115 / 255);
	DrawWall(x - bottomR, z - bottomR, x + bottomR, z - bottomR, average_height);
	glColor3d((double)0 / 255, (double)0 / 255, (double)0 / 255);
	DrawWindows(x - bottomR, z - bottomR, x + bottomR, z - bottomR, average_height, floors,-0.01);


	glColor3d((double)255 / 255, (double)191 / 255, (double)115 / 255);
	DrawWall(x + bottomR, z - bottomR, x + bottomR, z + bottomR, average_height);
	glColor3d((double)0 / 255, (double)0 / 255, (double)0 / 255);
	DrawWindows(x + bottomR, z - bottomR, x + bottomR, z + bottomR, average_height, floors,0.01);


	glColor3d((double)255 / 255, (double)191 / 255, (double)115 / 255);
	DrawWall(x - bottomR, z + bottomR, x + bottomR, z + bottomR, average_height);
	glColor3d((double)0 / 255, (double)0 / 255, (double)0 / 255);
	DrawWindows(x - bottomR, z + bottomR, x + bottomR, z + bottomR, average_height, floors,0.01);

	glColor3d(rgb.R , rgb.G, rgb.B);
	DrawRoof(x,z, average_height,1.5);
}
void DrawRoad(int x, int z)
{
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE); // or GL_MODULATE
	glBegin(GL_POLYGON);
	glTexCoord2d(0, 0);	glVertex3d(z - GSZ / 2, ground[x][z] + 0.001, x - GSZ / 2);
	glTexCoord2d(0, 2);	glVertex3d(z - GSZ / 2, ground[x - 1][z] + 0.001, x - 1 - GSZ / 2);
	glTexCoord2d(4, 2);	glVertex3d(z - 1 - GSZ / 2, ground[x - 1][z - 1] + 0.001, x - 1 - GSZ / 2);
	glTexCoord2d(4, 0);	glVertex3d(z - 1 - GSZ / 2, ground[x][z - 1] + 0.001, x - GSZ / 2);
	glEnd();

	glDisable(GL_TEXTURE_2D);
}
/*
This method checks whether the point is near a sea / lake / river
*/
bool isValidPoint(int centerX, int centerZ) {
	for (int dx = -maximumRiverDist; dx <= maximumRiverDist; dx++) {
		for (int dz = -maximumRiverDist; dz <= maximumRiverDist; dz++) {
			int distSquared = dx * dx + dz * dz;
			if (distSquared >= minimumRiverDist * minimumRiverDist && distSquared <= maximumRiverDist * maximumRiverDist) {
				int x = centerX + dx;
				int z = centerZ + dz;

				if (x < 0 || z < 0 || x >= GSZ || z >= GSZ)
					continue;

				// Checking if the point is valid for building
				if (ground[x][z] < (startingGround[x][z] - waterDepthFactor * startingGround[x][z]) || ground[x][z] < 0) {
					return true;
				}
			}
		}
	}
	return false;
}
/*
This function checks whether in the vicinity of x,z at a distance of bottomR
has a height lower than 0.6 (beginning of the sea)
or
Did start in this environment hydraulic demolition
or
Have you already built a building in this environment?*/
bool validPlaceToBuild(int x, int z)
{
	for (int i = -bottomR; i <= bottomR; i++)
		for (int j = -bottomR; j <= bottomR; j++)
			if (x + i >= 0 && x + i < GSZ && z + j >= 0 && z + j < GSZ)
				if (ground[x + i][z + j] <= 0.2
					|| ground[x + i][z + j] <= startingGround[x][z] - waterDepthFactor * startingGround[x][z]
					|| build[x + i][z + j] == true)
					return false;
	return true;
}
bool FloodFillIterative(int centerX, int centerZ)
{
	bool returnVal = false;
	vector <POINT3_INT> myStack;
	POINT3_INT current = { centerX,centerZ };
	myStack.push_back(current);
	RGB rgb = { ((double)(rand() % 256) / 256.0) ,((double)(rand() % 256) / 256.0) ,((double)(rand() % 256) / 256.0) };
	//int counter = 0;
	while (!myStack.empty())
	{
		//counter++;
		int radius = bottomR + 1;
		// 1. extract last element from stack
		current = myStack.back();
		myStack.pop_back();
		// 2. Checking conditions
		if (isValidPoint(current.x, current.z) && validPlaceToBuild(current.x,current.z))
		{
			for (int i = -radius; i <= radius; i++)
				for (int j = -radius; j <= radius; j++)
					if (current.x - i >= 0 && current.x - i < GSZ && current.z - j >= 0 && current.z - j < GSZ)
					{
						visited[current.x - i][current.z - j] = true;
						build[current.x - i][current.z - j] = true;
					}
			buildingHeights[current.x][current.z] = (rand() % max_height_building) + min_height_building;
			rgbVector[current.x][current.z] = rgb;
			returnVal = true;
			printf("x = %d ,z = %d buildingHeights = %d groud[x][z] = %.4f\n", current.x, current.z,  buildingHeights[current.x][current.z], ground[current.x][current.z]);
		}
		else {/*If it is a valid point then I will search in points with a distance from the building
If it is a point that is not vaild, I will search in neighboring points*/
			radius = 1;
		}
		for (int dx = -radius; dx <= radius; dx++)
		{
			for (int dz = -radius; dz <= radius; dz++)
			{
				int distSquared = dx * dx + dz * dz;
				if (distSquared >= radius * radius)
				{
					int x = current.x + dx;
					int z = current.z + dz;
					if (x >= 0 && x < GSZ && z >= 0 && z < GSZ)
					{
						if (!visited[x][z])
						{
							int distFromCenter = sqrt((x - centerX) * (x - centerX) + (z - centerZ) * (z - centerZ));
							if (startingGround[x][z] > 0)
								if (ground[x][z] > startingGround[x][z] - waterDepthFactor * startingGround[x][z])
									if (distFromCenter < maxDistFronCenter)
									{
										POINT3_INT point = { x,z };
										myStack.push_back(point);
										visited[x][z] = true; // Mark as visited
									}
						}
					}
				}
			}
		}
		
	}
	//printf("counter of stack in FF = %d \n", counter);
	return returnVal;
}


void DrawWheel(int n)
{
	double alpha, teta = 2 * PI / n;
	double x, y;

	glColor3d(0, 0, 0);

	glBegin(GL_LINES);
	for (alpha = 0; alpha <= 2 * PI; alpha += teta)
	{
		x = cos(alpha);
		y = sin(alpha);

		glVertex2d(x, y);
		glVertex2d(0, 0);

	}
	glEnd();

	glLineWidth(2);
	teta = 2 * PI / (n * 10);
	glBegin(GL_LINE_STRIP);
	for (alpha = 0; alpha <= 2 * PI; alpha += teta)
	{
		x = cos(alpha);
		y = sin(alpha);

		glVertex2d(x, y);

	}
	glEnd();
	glLineWidth(1);
}

void DrawBicycle()
{
	// front wheel
	glPushMatrix(); // the transformations that are inside PUSH/POP will not be applied oustside  PUSH/POP
	glTranslated(-0.12, 0, 0);
	glScaled(0.06, 0.06, 1);
	glRotated(angle, 0, 0, 1);
	DrawWheel(12);
	glPopMatrix();

	// rear wheel
	glPushMatrix(); // the transformations that are inside PUSH/POP will not be applied oustside  PUSH/POP
	glTranslated(0.12, 0, 0);
	glScaled(0.06, 0.06, 1);
	glRotated(angle, 0, 0, 1);
	DrawWheel(12);
	glPopMatrix();

	glPushMatrix(); // the transformations that are inside PUSH/POP will not be applied oustside  PUSH/POP

	glScaled(0.02, 0.02, 1);
	glRotated(angle, 0, 0, 1);
	DrawWheel(12);
	glPopMatrix();

	// frame
	glColor3d(1, 0.4, 0);
	glBegin(GL_LINE_LOOP);
	glVertex2d(0, 0);
	glVertex2d(0.12, 0);
	glVertex2d(0.04, 0.09);
	glEnd();

	glBegin(GL_LINE_STRIP);
	glVertex2d(0, 0);
	glVertex2d(-0.10, 0.08);
	glVertex2d(-0.12, 0);
	glVertex2d(-0.09, 0.12);
	glVertex2d(0.04, 0.09);
	glEnd();
	// handle
	glBegin(GL_LINE_STRIP);
	glVertex2d(-0.09, 0.12);
	glVertex2d(-0.08, 0.14);
	glVertex2d(-0.12, 0.14);
	glEnd();
	// saddle
	glBegin(GL_LINE_STRIP);
	glVertex2d(0.04, 0.09);
	glVertex2d(0.05, 0.11);
	glVertex2d(0.02, 0.11);
	glVertex2d(0.07, 0.11);
	glEnd();

}
void dispalyWorld()
{
	DrawFloor();
	if (start_rain == true)
	{
		resetConstruction();
		HydraulicErosionAct = true;
		for (int i = 0; i < 50; i++)
		{
			HydraulicErosion();
		}
	}
	if (build_city == true && HydraulicErosionAct == true)
	{
		for (int i = 0; i < 10; i++)
		{
			int counter = 0;
			int x, z;
			do {
				if (counter % 80 == 0)
				{
					HydraulicErosion();
				}
				else if (counter == 1000)
				{
					break;
				}
				do {
					x = rand() % GSZ;
					z = rand() % GSZ;
				} while (ground[x][z] < 0.2
					|| ground[x][z] < startingGround[x][z] - waterDepthFactor * startingGround[x][z]);
				counter++;
			} while (!FloodFillIterative(x, z));
		}
		build_city = false;
	}
	for (int i = 0; i < GSZ; i++)
	{
		for (int j = 0; j < GSZ; j++)
		{
			if (build[i][j] == true && buildingHeights[i][j] > 0) {
				DrawBuilding(i, j);
				glPushMatrix();
				if (i + 2 < GSZ && j+1 < GSZ && j-1 > 0)
				{
					if (ground[i + 2][j] > startingGround[i + 2][j] - waterDepthFactor * startingGround[i + 2][j])
					{
						glTranslated((double)(j - GSZ / 2) + 0.01, ground[i + 1][j] + 0.1, (double)(i - GSZ / 2) + 1.5);
						glScaled(3, 5, 3);
						DrawBicycle();
						glPopMatrix();
						DrawRoad(i + 2, j);
						DrawRoad(i + 2, j + 1);
						DrawRoad(i + 2, j - 1);
					}
				}
			}
			if (build[i][j] == true && buildingHeights[i][j] == 0) {
				//DrawRoad(i, j);
			}
		}
	}
}
void display()
{
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT); // clean frame buffer and Z-buffer
	glViewport(0, 0, W, H);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); // unity matrix of projection

	glFrustum(-1, 1, -1, 1, 0.75, 300);
	gluLookAt(eye.x, eye.y, eye.z,  // eye position
		eye.x+ sight_dir.x, eye.y-0.3, eye.z+sight_dir.z,  // sight dir
		0, 1, 0);


	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity(); // unity matrix of model

	dispalyWorld();

	glutSwapBuffers(); // show all
}
void displayTop()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // clean frame buffer and Z-buffer
	glViewport(0, 0, W, H);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); // unity matrix of projection

	glFrustum(-1, 1, -1, 1, 0.75, 300);
	gluLookAt(eye.x, 70, eye.z,  // eye position
		eye.x + sight_dir.x, 60, eye.z + sight_dir.z,  // sight dir
		0, 1, 0);


	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity(); // unity matrix of model
	dispalyWorld();

	glutSwapBuffers(); // show all
}
void idle() 
{
	int i, j;
	double dist;
	angle += 0.1;

	// ego-motion  or locomotion
	sight_angle += angular_speed;
	// the direction of our sight (forward)
	sight_dir.x = sin(sight_angle);
	sight_dir.z = cos(sight_angle);
	// motion
	eye.x += speed * sight_dir.x;
	eye.y += speed * sight_dir.y;
	eye.z += speed * sight_dir.z;


	glutPostRedisplay();
}


void SpecialKeys(int key, int x, int y)
{
	switch (key)
	{
	case GLUT_KEY_LEFT: // turns the user to the left
		angular_speed += 0.0004;
		break;
	case GLUT_KEY_RIGHT:
		angular_speed -= 0.0004;
		break;
	case GLUT_KEY_UP: // increases the speed
		speed+= 0.005;
		break;
	case GLUT_KEY_DOWN:
		speed -= 0.005;
		break;
	case GLUT_KEY_PAGE_UP:
		eye.y += 0.1;
		break;
	case GLUT_KEY_PAGE_DOWN:
		eye.y -= 0.1;
		break;

	}
}
void menu(int choice)
{
	switch (choice)
	{
	case 1: // regular view
		glutDisplayFunc(display);
		break;
	case 2: // top view
		glutDisplayFunc(displayTop);
		break;
	case 3:
		if (HydraulicErosionAct == true)
		{
			start_rain = false;
			build_city = true;
		}
		break;
	case 4:
		start_rain = false;
		HydraulicErosionAct = false;
		createWorld();
		copyStartingWorld();
		resetConstruction();
		break;
	case 5:
		start_rain = false;
		HydraulicErosionAct = false;
		for (int i = 0; i < GSZ; i++)
			for (int j = 0; j < GSZ; j++)
				ground[i][j] = startingGround[i][j];
		resetConstruction();
		break;
	}
}
void Mouse(int button, int state, int x, int y) {
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
		start_rain = !start_rain;
	}
}
void main(int argc, char* argv[]) 
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE|GLUT_DEPTH);
	glutInitWindowSize(W, H);
	glutInitWindowPosition(400, 100);
	glutCreateWindow("Final Project");

	glutDisplayFunc(display);
	glutIdleFunc(idle);
	glutSpecialFunc(SpecialKeys);
	glutMouseFunc(Mouse);
	glutCreateMenu(menu);
	glutAddMenuEntry("Regular view", 1);
	glutAddMenuEntry("Top view", 2);
	glutAddMenuEntry("Build city", 3);
	glutAddMenuEntry("Reset world", 4);
	glutAddMenuEntry("Reset to starting world", 5);
	glutAttachMenu(GLUT_RIGHT_BUTTON);
	init();
	glutMainLoop();
}