#ifndef FUNCTIONPROTOTYPESHEADER_H_
#define FUNCTIONPROTOTYPESHEADER_H_
#pragma once

#include "GlobalProjectHeader.h"
#include "StructsHeader.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <math.h>



/* --------------------------- File I/O Functions ----------------------------------*/
void checkFile(const FILE* f);
void RandomInitializeFileData(const char* fileName, int numofPoints, int numberofCoords);
void writeToFile(int numofClusters, Cluster *clusterArray, double elapsedTime, double quality);
void clusterFindError();
Point* readPointsDataSetFromFile(const char* fileName, double* t, double* dt, int* limit, double* quality_measure, int* numofPoints, int* num_clusters);


/* --------------------------- Printing Functions ----------------------------------*/
void printArray(Point* allPoints, int size, int processnumber);
void printGoodClusters(double quality_found, double time, Cluster* clustersArray, int num_clusters, int id);
void printClusters(Cluster* clustersArray, int num_clusters, int id);



/* ---------------------------- (MPI,OMP) K-Means Algorithm Functions ----------------------------------*/
MPI_Datatype createMPIDataTypePoint();
MPI_Datatype createMPIDataTypeCluster();
Point* dividePointsDataSet(int sizeOfAllPoints, Point *arrPoints, int* numofPointsPerRank, Point *PointsPerProcess, int myRank, int numOfProcs, MPI_Datatype MPI_POINT_TYPE);
Cluster* kmeansFindFinalClusters(double* qualityfound, double quality_measure, int* numofPoints, int* num_clusters, double t, double dt, int limit, Point* allpoints, int* numofPointsPerProcess, Point*  pointsPerProcess, int* myid, int* numprocs, MPI_Datatype  MPI_POINT_TYPE, MPI_Datatype MPI_CLUSTER_TYPE, Cluster* allClusters);
Cluster* initializeClustersDatafields(int numOfClusters, Point *pointsArray, Cluster* allClusters, MPI_Datatype MPI_CLUSTER_TYPE, int myRank);

void initialInit(Cluster* arr, int numClusters);
void parallel_kmeans(Cluster* clustersArr, int numClusters, Point* pointsPerProcess, int num_of_points, int limit, int* myid);
void broadcastData(int* numofpoints, int* numclusters, double* t, double *dt, int *limit, double* quality);
void collectPoints(int myRank, Point *allPoints, int numOfAllPoints, int numPointsPerProcess, Point *PointsForProcess, MPI_Datatype MPI_POINT_TYPE, int numOfProcs);
void collectPointsCenters(Point* pointsPerProcess, Cluster* clustersArr, int num_of_points, int numClusters);
void calculateQuality(int numberOfPoints, int numberOfClusters, Point* pointsArray, Cluster* clusterArray, double *quality);
void calculateAverage(int numClusters, Cluster* arrClusters);
void toContinueKMeans(int* continuekflag);
void giveClustersDiameter(Point* pointsArray, int numberOfPoints, Cluster* clusterArray, int numOfClusters);
void pointsLocationRemainderOMP(Point* pointsArray, double theTime,int remainderOffset, int totalPointsSize);
int assign_point_membership(Cluster* clusters, int numClusters, int *myid, Point point);
double clusterDiameterFinder(Point* pointsArray, int numberOfPoints, int clusterId, int numOfClusters);
double distanceClusters(Cluster* c1, Cluster* c2);
double distancePoints(Point* point1, Point* point2);
double euclid_dist_2(Point p, Cluster c);


/*---------------------------- Miscellaneous Functions ----------------------------*/
void freeAll(Point* allpoints, Cluster* clustersArray, Point* pointsBufferForEachProcess);
void checkDynamicAllocation(const void* ptr);
void InitInputFileHelper(FILE* f, int numberOfValues, int numberofCoordinates);
double doubleRandomGenerator(double min, double max);


/*-----------------------------Cuda Functions -------------------------------------*/
cudaError_t pointsLocationWholeCuda(int allPointsSize, double theTime, Point* pointsArray, unsigned int threadsperblock, int* remainderOffset);

#endif

