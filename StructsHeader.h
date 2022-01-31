#ifndef STRUCTSHEADER_H_
#define STRUCTSHEADER_H_
#pragma once
/* --------------------------- Structs ---------------------------*/

typedef struct Point
{
	double x;
	double y;
	double z;
	double Vx;
	double Vy;
	double Vz;
	int idCluster;
} Point;

typedef struct Cluster
{
	double x;
	double y;
	double z;
	double diameter;
	int id;
	int numPointsInCluster;
} Cluster;

#endif