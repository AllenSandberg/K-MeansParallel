#include "FunctionPrototypesHeader.h"
#include "GlobalProjectHeader.h"
#include "StructsHeader.h"



void checkFile(const FILE* f)
{
	if (f == NULL)
	{
		puts("Could not open the file , please check file path");
		MPI_Finalize();
		exit(1);
	}
}

void RandomInitializeFileData(const char* fileName, int numofPoints, int numofCoordinates)
{
	FILE* f = fopen(FILE_INPUT, "w");
	checkFile(f);
	fprintf(f, "%d %d %lf %lf %d %lf\n", NUMPOINTS, NUMCLUSTERS, T, DT, LIMIT, QUALITY_MEASURE);
	InitInputFileHelper(f, numofPoints, NUMOFCOORDINATES);
}


void InitInputFileHelper(FILE* f, int numberOfPoints, int coords)
{
	for (int i = 0; i < numberOfPoints; i++)
	{
		double x = doubleRandomGenerator(((-1) * (numberOfPoints)), (numberOfPoints));
		double y = doubleRandomGenerator(((-1) * (numberOfPoints)), (numberOfPoints));
		double z = doubleRandomGenerator(((-1)* (numberOfPoints)), (numberOfPoints));
		double Vx = doubleRandomGenerator(((-1)*(numberOfPoints / 100)), (numberOfPoints / 10));
		double Vy = doubleRandomGenerator(((-1)*(numberOfPoints / 10)), (numberOfPoints / 100));
		double Vz = doubleRandomGenerator(((-1)*(numberOfPoints / 50)), (numberOfPoints / 10));

		fprintf(f, "%lf %lf %lf \n ", x, y, z);
		fprintf(f, "%lf %lf %lf \n", Vx, Vy, Vz);
	}
	fclose(f);
}

double doubleRandomGenerator(double min, double max)
{
	double range = (max - min);
	double div = RAND_MAX / range;
	return min + (rand() / div);
}




void writeToFile(int numofClusters, Cluster *clusterArray, double elapsedTime, double quality)
{
	int i;

	FILE *fp = fopen(FILE_OUT, "w");

	if (fp == NULL) {
		printf("the file is not open\n");
		exit(1);
	}

	fprintf(fp, "First occurrence at t = %lf with q = %lf\nCenters of the clusters:\n", elapsedTime, quality);
	for (i = 0; i < numofClusters; i++)
	{
		fprintf(fp, " Cluster Centre ID : [ %d ]  ( %lf , %lf,  %lf )\n", clusterArray[i].id, clusterArray[i].x, clusterArray[i].y, clusterArray[i].z);
	}

	fclose(fp);
}
//(const char* fileName, double* t, double* dt, int* limit, double* quality_measure, int* numofPoints, int* num_clusters)
Point* readPointsDataSetFromFile(const char* fileName, double* t, double* dt, int* limit, double* quality_measure, int* numofPoints, int* num_clusters)
{
	FILE* f = fopen(fileName, "r");
	checkFile(f);

	fscanf(f, "%d %d %lf %lf %d %lf ", numofPoints, num_clusters, t, dt, limit, quality_measure);
	//printf("\n After first line has been read \n NumofPoints  %d  NumOfClusters %d  T %lf  dt %lf limit %d  quality_measure %lf \n", *numofPoints, *num_clusters, *t, *dt, *limit, *quality_measure);

	Point* allpoints = (Point*)calloc(*numofPoints, sizeof(Point));
	checkDynamicAllocation(allpoints);

	for (int i = 0; i < *numofPoints; i++)
	{
		/* Read from file Xi , Yi , Zi to allpoints[i]*/
		fscanf(f, "%lf %lf %lf ", &allpoints[i].x, &allpoints[i].y, &allpoints[i].z);
		/* Read from file Vxi , Vyi , Vzi to allpoints[i]*/
		fscanf(f, "%lf %lf %lf ", &allpoints[i].Vx, &allpoints[i].Vy, &allpoints[i].Vz);
	}
	fclose(f);
	return allpoints;
}

void clusterFindError()
{
	FILE* f = fopen(FILE_OUT, "w");
	fprintf(f, "\n\n\n Error in finding Good Clusters for input data ! .  \n\n\n");
}

