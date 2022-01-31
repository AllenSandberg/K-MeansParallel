#include "FunctionPrototypesHeader.h"
#include "GlobalProjectHeader.h"
#include "StructsHeader.h"


void printArray(Point* allPoints, int size, int processnumber)
{
	printf("\n\n\n");
	fflush(stdout);
	for (int i = 0; i < size; i++)
	{
		printf("processnumber :  %d  allPoints[%d].x = %lf    allPoints[%d].y = %lf     allPoints[%d].z = %lf \n\n", processnumber, i, allPoints[i].x, i, allPoints[i].y, i, allPoints[i].z);
		fflush(stdout);
		printf("processnumber :  %d  allPoints[%d].Vx = %lf     allPoints[%d].Vy =  %lf      allPoints[%d].Vz =  %lf     ----->  This Point's membership to Cluster ID [%d]  \n\n\n\n", processnumber, i, allPoints[i].Vx, i, allPoints[i].Vy, i, allPoints[i].Vz, allPoints[i].idCluster);
		fflush(stdout);
	}
	printf("\n\n\n");
	fflush(stdout);
}


void printGoodClusters(double quality_found, double time, Cluster* clustersArray, int num_clusters, int id)
{
	printf("\n First Occurence t = %lf with q = %lf \n\n ", time, quality_found);
	fflush(stdout);
	printf("Centers of the clusters: \n \n");
	fflush(stdout);

	printClusters(clustersArray, num_clusters, id);
}


void printClusters(Cluster* clustersArray, int num_clusters, int myid)
{
	printf("\n\n\n\n");
	fflush(stdout);

	for (int i = 0; i < num_clusters; i++)
	{
		printf(" \n  Process [ %d ]  Cluster [ %d ] %lf  %lf  %lf \n ", myid, i + 1, clustersArray[i].x, clustersArray[i].y, clustersArray[i].z);
		fflush(stdout);
	}

	printf("\n\n\n\n");
	fflush(stdout);

}