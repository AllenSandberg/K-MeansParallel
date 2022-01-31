#include "GlobalProjectHeader.h"
#include "StructsHeader.h"
#include "FunctionPrototypesHeader.h"
#define _CRT_SECURE_NO_WARNINGS

//VALID

int main(int argc, char *argv[])
{
	int  numberofprocesses, processid;
	int limit; // maximum number of iterations for K-Means Algorithm.
	int input_numofPoints;
	int input_numofClusters;
	double t = 0; // (T in project guidelines) . 
	double dt = 0;  // dt - delta-interval block .
	double inputQM = 0; // Input quality measure to stop (from given Input Data File).
	double quality_found = 0;


	Point* allpoints = NULL; // whole dataset of input points.

	//dataset chunk - subgroup of points for each process.
	// size =  (input_numofPoints / numprocesses) +  (input_numofPoints % numprocesses) for process MASTER.
	// size =  (input_numofPoints / numprocesses) for each Slave process [1..n-1] .
	Point* subgroup_ofpointsForEachProcess = NULL;
	Cluster* clustersArray = NULL; // Array of Clusters.
	int numofPointsPerRank;


	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &processid);
	MPI_Comm_size(MPI_COMM_WORLD, &numberofprocesses);

	MPI_Datatype MPI_POINT_TYPE = createMPIDataTypePoint();
	MPI_Datatype MPI_CLUSTER_TYPE = createMPIDataTypeCluster();

	/* Master Job is to read File Data - [N    K    T   dT   LIMIT   QM] & [Points' Data Values]. */
	if (processid == MASTER) {
		//RandomInitializeFileData(FILE_INPUT, NUMPOINTS, NUMOFCOORDINATES);
		allpoints = readPointsDataSetFromFile(FILE_INPUT, &t, &dt, &limit, &inputQM, &input_numofPoints, &input_numofClusters);
		printf("\n input_numofPoints = %d , input_numofClusters = %d , t = %lf , dt = %lf , limit = %d , inputQM = %lf \n ", input_numofPoints, input_numofClusters, t, dt, limit, inputQM);
		fflush(stdout);
		//printArray(allpoints, input_numofPoints, MASTER);
		printf("\n ----- Main calling readPointsFromFile function -------  quality_measure %lf   dt = %lf  \n ", quality_found, dt);
		fflush(stdout);
	}

	// Master Contains input_numofPoints,input_numofClusters,t,dt,limit & quality_measure. It broadcasts this first line of values to all other processes [1..numberofprocesses-1].  (Broadcast to Slaves 1st line of file)
	broadcastData(&input_numofPoints, &input_numofClusters, &t, &dt, &limit, &inputQM);

	// Scatter allpoints to all Processes - Master & Slaves. Now each process contains its own subgroup of data. ( first line of values from input file + share of points'  & numofPoints for each process).
	subgroup_ofpointsForEachProcess = dividePointsDataSet(input_numofPoints, allpoints, &numofPointsPerRank, subgroup_ofpointsForEachProcess, processid, numberofprocesses, MPI_POINT_TYPE);

	clustersArray = (Cluster*)calloc(input_numofClusters, sizeof(Cluster));
	checkDynamicAllocation(clustersArray);

	/*Master Initializes first K Points as Cluster Centroids & broadcasts all of these first K Initialized Cluster Centroids to all slave processes (Processes 1..n-1)*/
	clustersArray = initializeClustersDatafields(input_numofClusters, allpoints, clustersArray, MPI_CLUSTER_TYPE, processid);
	// All Processes its own subgroup of data. ( first line of values from input file + share of points'  & numofPoints for each process + K initialized Cluster Centroid Points ).
	clustersArray = kmeansFindFinalClusters(&quality_found, inputQM, &input_numofPoints, &input_numofClusters, t, dt, limit, allpoints, &numofPointsPerRank, subgroup_ofpointsForEachProcess, &processid, &numberofprocesses, MPI_POINT_TYPE, MPI_CLUSTER_TYPE, clustersArray);
	printf("\n ---------------------------------------Main: Return from checkForGoodClusters() ----   Final ClustersArray of Cluster Centres given from process id = %d is returned to Main() to be given as input param to freeAll function--------------------------------------- \n", processid);
	fflush(stdout);
	freeAll(allpoints, clustersArray, subgroup_ofpointsForEachProcess);
	MPI_Finalize();
	return 0;
}
