#include "FunctionPrototypesHeader.h"
#include "GlobalProjectHeader.h"
#include "StructsHeader.h"



MPI_Datatype createMPIDataTypePoint()
{
	Point point;
	MPI_Datatype MPI_Point;
	MPI_Datatype allType[7] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
		MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT };
	int blocklen[7] = { 1, 1, 1 , 1 , 1, 1, 1 };
	MPI_Aint disp[7];

	disp[0] = (char *)&point.x - (char *)&point;
	disp[1] = (char *)&point.y - (char *)&point;
	disp[2] = (char *)&point.z - (char *)&point;
	disp[3] = (char *)&point.Vx - (char *)&point;
	disp[4] = (char *)&point.Vy - (char *)&point;
	disp[5] = (char *)&point.Vz - (char *)&point;
	disp[6] = (char *)&point.idCluster - (char *)&point;
	MPI_Type_create_struct(7, blocklen, disp, allType, &MPI_Point);
	MPI_Type_commit(&MPI_Point);
	return MPI_Point;
}


MPI_Datatype createMPIDataTypeCluster()
{
	Cluster cluster;
	MPI_Datatype typeForMPICluster;
	MPI_Datatype type[6] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE , MPI_INT ,MPI_INT };
	int blocklen[6] = { 1, 1, 1 , 1 , 1, 1 };
	MPI_Aint disp[6];

	disp[0] = (char *)&cluster.x - (char *)&cluster;
	disp[1] = (char *)&cluster.y - (char *)&cluster;
	disp[2] = (char *)&cluster.z - (char *)&cluster;
	disp[3] = (char *)&cluster.diameter - (char *)&cluster;
	disp[4] = (char *)&cluster.id - (char *)&cluster;
	disp[5] = (char *)&cluster.numPointsInCluster - (char *)&cluster;
	MPI_Type_create_struct(6, blocklen, disp, type, &typeForMPICluster);
	MPI_Type_commit(&typeForMPICluster);
	return typeForMPICluster;

}



Cluster* kmeansFindFinalClusters(double* qualityfound, double quality_measure, int* numofPoints, int* num_clusters, double t, double dt, int limit, Point* allpoints, int* numofPointsPerProcess, Point*  pointsPerProcess, int* myid, int* numprocs, MPI_Datatype  MPI_POINT_TYPE, MPI_Datatype MPI_CLUSTER_TYPE, Cluster* allClusters) {
	double quality = 0;
	double elapsedTime = 0;
	int iterationNumber = 0;
	int remainderOffset = 0;


	/*Iterate over different timeframes.
	Per Each timeframe, group of points is moved according to velocity formula to provide different initial point disposition states
	for K-Means' algorithm to run on.*/
	for (double timeframe = 0; timeframe <= t; timeframe += dt)
	{

		/*Input parameters: whole cluster array,size of all clusters, pointsPerProcess, numofPointsPerProcess , limit iterations and dt.
		Here for each process's subset of points we update the points' memberships,
		and as a result update the clusters' centres for each process locallly.
		We perform this parallel k-means implementation until no points' memberships have changed between two iterations
		(clusters stabilize) or until LIMIT num of iterations if clusters don't stabilize.
		*/

		parallel_kmeans(allClusters, *num_clusters, pointsPerProcess, *numofPointsPerProcess, limit, myid);


		collectPoints(*myid, allpoints, *numofPoints, *numofPointsPerProcess, pointsPerProcess, MPI_POINT_TYPE, *numprocs);
		elapsedTime = iterationNumber*dt;

		if (*myid == MASTER)
		{
			// Now Master Calculates quality
			printf("\n\nNow Master Calculates quality iteration Number = %d\n\n", iterationNumber);
			fflush(stdout);

			calculateQuality(*numofPoints, *num_clusters, allpoints, allClusters, &quality);
			// if in this specific dt (logical "external" K-Means Iteration) the quality that has been found is less than or equal to quality_measure given in inputfile.
			if (quality <= quality_measure)
			{
				printf("\n\n -------------------- best quality is found to be = %lf  quality_measure =  %lf ------------------- n\n   ", quality, quality_measure);
				fflush(stdout);
				printGoodClusters(quality, elapsedTime, allClusters, *num_clusters, *myid);
				writeToFile(*num_clusters, allClusters, elapsedTime, quality);
			}
			printf("\n\n------------------quality = %lf------------------ t = %lf     dt = %lf     \n\n", quality, t, dt);
			fflush(stdout);
		}
		MPI_Bcast(&quality, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		if (quality <= quality_measure)
			return allClusters;


		/*Each Process enters Cuda Function with its own set of points*/
		//printf("\n\n Now all processes will perform 	pointsLocationWholeCuda  at  elapsed time  = %lf  sending dt as constant value = %lf \n\n", elapsedTime, dt);
		//fflush(stdout);
		pointsLocationWholeCuda(*numofPointsPerProcess, dt, pointsPerProcess, THREADS_PER_BLOCK, &remainderOffset);

		int numblocks = (*numofPointsPerProcess / THREADS_PER_BLOCK);
		remainderOffset = numblocks* THREADS_PER_BLOCK;
		// if modulo
		if (((*numofPointsPerProcess) % (THREADS_PER_BLOCK)) > 0)
			pointsLocationRemainderOMP(pointsPerProcess, dt, remainderOffset, *numofPointsPerProcess);

		iterationNumber++;

	}
	clusterFindError();
	return allClusters;
}


Cluster* initializeClustersDatafields(int numOfClusters, Point *pointsArray, Cluster* allClusters, MPI_Datatype MPI_CLUSTER_TYPE, int myRank)
{
	int i;

	/*Master initializes all clusters' datafields. First K (numOfClusters) Points in Order of PointsArray will be assigned as initial Clusters. */
	if (myRank == MASTER)
	{
#pragma omp parallel for
		for (i = 0; i < numOfClusters; i++)
		{
			allClusters[i].x = pointsArray[i].x;
			allClusters[i].y = pointsArray[i].y;
			allClusters[i].z = pointsArray[i].z;
			allClusters[i].diameter = 0;
			allClusters[i].id = i;
			allClusters[i].numPointsInCluster = 0;
		}
	}
	//Master broadcasts all first K Initial Points- Clusters' datafields to all salves. (Processes 1..n-1)
	// All Slaves will now contain the first K Initial Cluster Centroids.
	MPI_Bcast(allClusters, numOfClusters, MPI_CLUSTER_TYPE, MASTER, MPI_COMM_WORLD);
	return allClusters;

}

void initialInit(Cluster* arr, int numClusters)
{
	int i;
#pragma omp parallel for
	for (i = 0; i < numClusters; i++)
		arr[i].numPointsInCluster = 0;
}


/*Build PointsPerProcess - chunk of points array for Processes  [0..numOfProcesses-1]
Division of Whole Dataset of points goes as follows:
- Master will contain wholesize + remaindersize points.
- Slaves will contain wholesize points.
*/

Point* dividePointsDataSet(int sizeOfAllPoints, Point *arrPoints, int* numofPointsPerRank, Point *PointsPerProcess, int myRank, int numOfProcesses, MPI_Datatype MPI_POINT_TYPE)
{

	int i, wholesize, remaindersize;
	MPI_Status status;

	wholesize = sizeOfAllPoints / numOfProcesses;
	remaindersize = sizeOfAllPoints % numOfProcesses;

	if (myRank == MASTER)
		*numofPointsPerRank = wholesize + remaindersize;
	else
		*numofPointsPerRank = wholesize;

	PointsPerProcess = (Point*)malloc(sizeof(Point)*(*numofPointsPerRank));
	checkDynamicAllocation(PointsPerProcess);

	/* Master Job */
	if (myRank == MASTER)
	{
		//  Master sends to himself his share of points in Parallel in openMP until (wholesize + remaindersize).
		// In a non-sequential (parallel) order , each thread assigns a point from arrPoints[i] to Master's Point* subgroup of points.
#pragma omp parallel for
		for (i = 0; i < *numofPointsPerRank; i++)
		{
			(PointsPerProcess)[i] = arrPoints[i];
		}
		//  Master sends to all other processes their share of points [ 1.. numOfProcesses-1]
		for (i = 1; i < numOfProcesses; i++)
		{
			MPI_Send(arrPoints + (*numofPointsPerRank + (wholesize*(i - 1))), wholesize, MPI_POINT_TYPE, i, MPI_FLAG, MPI_COMM_WORLD);

		}
	}

	/* Slave Job */
	else {
		// Other Processes recieve their share of points from Master.
		MPI_Recv(PointsPerProcess, wholesize, MPI_POINT_TYPE, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	}
	// After Master finished Scattering all points to all processes, Function will return for each process its own share of points.
	return PointsPerProcess;
}


/* parallel_kmeans.
This function calls sub-operations (implemented in both MPI & OpenMP).
This is the core parallel-kmeans algorithm implementation of the project.
Goal :
Per given clustersArr,numClusters, pointsPerProcess (buffer of points / subgroup of points) for specific process ,
num_of_points for specific process, limit - num of K-Means Iterations & ProcessID.
*/
void parallel_kmeans(Cluster* clustersArr, int numClusters, Point* pointsPerProcess, int num_of_points, int limit, int* myid)
{
	int nearestClusterID, pointIndex, iternal_iterationNumber = 0;
	int unstabilizedClusters_Flag = TRUE;

	/*If not all clusters per given process are stabilized (points' membership have changed)
	and num of iterations is less than LIMIT iterations*/
	while (unstabilizedClusters_Flag && iternal_iterationNumber < limit)
	{
		unstabilizedClusters_Flag = FALSE;
		// Initialize for each Process - its clusters' numberOfPoints to be zero for all clusters
		// in the given process that performs this function.
		initialInit(clustersArr, numClusters);

		/* Per specific subset of points given by any process in MPI executing this function - Iterate over this subset of points.
		This is the Parallel Job, each process does it parallel*/

		for (pointIndex = 0; pointIndex < num_of_points; pointIndex++)
		{
			//For each process - per its given specific points in its buffer of points Array. We will find its nearest cluster.
			nearestClusterID = assign_point_membership(clustersArr, numClusters, myid, pointsPerProcess[pointIndex]);
			// If specific point is now closer to a new found Cluster Centroid.
			if (pointsPerProcess[pointIndex].idCluster != nearestClusterID)
			{
				pointsPerProcess[pointIndex].idCluster = nearestClusterID;
				unstabilizedClusters_Flag = TRUE;
			}
			// Either way - if point's membership is changed or not , per specific cluster 
			// We need to increment the specific process' cluster's num of points. (check why)
			clustersArr[nearestClusterID].numPointsInCluster++;
		}
		// Now per specific process - each cluster in the clustersArr of the local process is updated with its numofpoints 
		// and each point in Process's pointsBuffer is updated with its clusterID membership.

		//Now, for each process' subgroup of points , as a result of at least a single point that has changed its membership 
		// - we will calculate the new total x,y,z of local Cluster Centroids - per given process .
		collectPointsCenters(pointsPerProcess, clustersArr, num_of_points, numClusters);

		// all Processes will contain updated Cluster Centres.
		calculateAverage(numClusters, clustersArr);

		// check if all processes' unstabilizedClusters_Flags are equal to 0. 
		//If so , this function's flag = FALSE (all clusters per given process are stabilized)
		// else - not all clusters per given process are stabilized (points' membership have changed) 
		// - thus continue parallel_kmeans() operation.
		toContinueKMeans(&unstabilizedClusters_Flag);

		iternal_iterationNumber++;
	}
}

// VALID

/* Master Process Broadcasts to all processes : numofpoints,numclusters,t,dt,limit*/
void broadcastData(int* numofpoints, int* numclusters, double* t, double *dt, int *limit, double* qualitymeasure)
{
	MPI_Bcast(numofpoints, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(numclusters, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(t, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(limit, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(dt, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(qualitymeasure, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

}




/*Collect all Data to Master. Now Master has all points & updated final Cluster Centroids.*/
void collectPoints(int myRank, Point *allPoints, int numOfAllPoints, int numPointsPerProcess, Point *PointsForProcess, MPI_Datatype MPI_POINT_TYPE, int numOfProcs)
{
	int i;
	int pointsDeltaOffset = numOfAllPoints / numOfProcs;
	if (myRank == MASTER)
	{

#pragma omp parallel for 
		for (i = 0; i < numPointsPerProcess; i++)
		{
			allPoints[i] = PointsForProcess[i];
		}
		MPI_Status status;

		/* Master recieves points from each other process.*/
#pragma omp parallel for 
		for (i = 1; i < numOfProcs; i++)
		{
			MPI_Recv(allPoints + (numPointsPerProcess + (pointsDeltaOffset*(i - 1))), pointsDeltaOffset, MPI_POINT_TYPE, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		}

	}
	/* Other processes send their share of points to master.*/
	else {
		MPI_Send(PointsForProcess, numPointsPerProcess, MPI_POINT_TYPE, MASTER, 0, MPI_COMM_WORLD);
	}

}



void collectPointsCenters(Point* pointsPerProcess, Cluster* clustersArr, int num_of_points, int numClusters)
{
	int i;

	// Initialize each Cluster's x,y,z 
#pragma omp parallel for
	for (i = 0; i < numClusters; i++) {
		clustersArr[i].x = 0;
		clustersArr[i].y = 0;
		clustersArr[i].z = 0;
	}
	// Calculate Total X,Y,Z for each Cluster in Process - these values are to be referred for future calculations.
	for (i = 0; i < num_of_points; i++) {
		clustersArr[pointsPerProcess[i].idCluster].x += pointsPerProcess[i].x;
		clustersArr[pointsPerProcess[i].idCluster].y += pointsPerProcess[i].y;
		clustersArr[pointsPerProcess[i].idCluster].z += pointsPerProcess[i].z;
	}
}


void calculateQuality(int numberOfPoints, int numberOfClusters, Point* pointsArray, Cluster* clusterArray, double *quality)
{
	int i, j;
	double distance = 0;
	double theQual = 0;
	int denominator = numberOfClusters*(numberOfClusters - 1);

	//The diameter of a cluster is the maximum distance between any two points of the cluster.
	//void giveClustersDiameter(Point* pointsArray, int numberOfPoints, Cluster* clusterArray, int numOfClusters);
	// void giveClustersDiameter(Point* pointsArray, int numberOfPoints, Cluster* clusterArray, int numOfClusters)
	giveClustersDiameter(pointsArray, numberOfPoints, clusterArray, numberOfClusters);

	// Now each Cluster Contains its calculated diameter. Diameters d[0]..d[k-1] are calculated.

#pragma omp parallel for private(distance,j) reduction(+ : theQual)
	for (i = 0; i < numberOfClusters; i++)
	{
		// Each Cluster has its own private distance variable.
		for (j = 0; j < numberOfClusters; j++)
		{
			if (i != j)
			{
				distance = distanceClusters(&clusterArray[i], &clusterArray[j]);
				// add to the Qual the clusterArray[i].diameter/distance(Cluster[i]- Cluster[j] for all)
				theQual += (clusterArray[i].diameter / distance);  
			}
		}
	}
	theQual = theQual / denominator;
	*quality = theQual;

}

// VALID

/* MPI_Allreduce will reduce the values and distribute the results to all processes. */
void calculateAverage(int numClusters, Cluster* arrClusters)
{
	int i, total = 0;
	double x = 0, y = 0, z = 0;

	for (i = 0; i < numClusters; i++)
	{
		/* Sum arrClusters[i].x from each process. Operation Description as follows:
		Process (0) cluster[i]'s total X value of points assigned to cluster[i]
		+   Process (1)  cluster[i]'s total X value of points assigned to cluster[i]
		+... Process (n-1) cluster[i]'s total X value of points assigned to cluster[i]
		*/
		MPI_Allreduce(&(arrClusters[i].x), &x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);



		/* Sum arrClusters[i].y from each process. Operation Description as follows:
		Process (0) cluster[i]'s total Y value of points assigned to cluster[i]
		+   Process (1)  cluster[i]'s total Y value of points assigned to cluster[i]
		+... Process (n-1) cluster[i]'s total Y value of points assigned to cluster[i]
		*/

		MPI_Allreduce(&(arrClusters[i].y), &y, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);



		/* Sum arrClusters[i].z from each process. Operation Description as follows:
		Process (0) cluster[i]'s total Z value of points assigned to cluster[i]
		+   Process (1)  cluster[i]'s total Z value of points assigned to cluster[i]
		+... Process (n-1) cluster[i]'s total Z value of points assigned to cluster[i]
		*/

		MPI_Allreduce(&(arrClusters[i].z), &z, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		// Now all Processes will have [Sigma (x), Sigma (y), Sigma (z)] for each arrClusters[i].
		// Sum arrClusters[i].numPointsInCluster for each arrClusters[i] in each process. Send total numPointsinCluster to all Processes.

		MPI_Allreduce(&(arrClusters[i].numPointsInCluster), &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		/* Each Process Calculate new Clusters' Centre of Gravity. Now all Processes will contain updated Cluster Centres.
		All processes will have average of Sima(x) , Sigma (y), Sigma (z) for each arrCluster[i].
		*/

		//printf("\n\n\n --------------------------------------------     x = %lf   y = %lf     z  =  %lf --------------------------------------------   \n\n\n", x, y, z);
		//fflush(stdout);
		arrClusters[i].x = x / total;
		arrClusters[i].y = y / total;
		arrClusters[i].z = z / total;
	}
}

/*Input unstabilizedClusters_Flag.*/

void toContinueKMeans(int* unstabilizedClusters_Flag)
{
	int sumflags;
	// Send Sum of unstabilizedClusters_Flag of all processes to sumflags variable. 
	MPI_Allreduce(unstabilizedClusters_Flag, &sumflags, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	/*If any of the processes' unstabilizedClusters_Flag is 1 - then sum >=1. This means that at least 
	one of the processes' unstabilizedClusters_Flag == 1
	--> we need to continue parallel_kmeans because for a given process we have found unstabilized clusters 
	(in other words , points' memberships have changed)
	*/
	if (sumflags >= 1)
		*unstabilizedClusters_Flag = TRUE;
	/*Else all processes' unstabilizedClusters_Flag are 0 - then sum < 1. Then we know for sure that we need to stop 
	for each process its parallel_kmeans operation.
	Because only in this case - it means that all processes' new found cluster centroids have been found 
	to be "stabilized" (points' memberships don't change for next iteration).
	*/
	else
		*unstabilizedClusters_Flag = FALSE;
}


void giveClustersDiameter(Point* pointsArray, int numberOfPoints, Cluster* clusterArray, int numOfClusters)
{
	int i;
	// For each cluster in #pragma omp parallel for provide its diameter via its 
#pragma omp parallel for
	for (i = 0; i < numOfClusters; i++)
	{
		clusterArray[i].diameter = clusterDiameterFinder(pointsArray, numberOfPoints, clusterArray[i].id, numOfClusters);
	}
}

double clusterDiameterFinder(Point* pointsArray, int numberOfPoints, int clusterId, int numOfClusters)
{
	int i, j;
	double distance;
	double TheMaxDistance = 0;

	// Find Maximum distance from pointsArray[i] to pointsArray[j]
	for (i = 0; i < numberOfPoints; i++)
	{
		if (pointsArray[i].idCluster == clusterId)
		{
			for (j = (i + 1); j < numberOfPoints; j++)
			{
				if (pointsArray[j].idCluster == clusterId)
				{
					distance = distancePoints(&pointsArray[i], &pointsArray[j]);
					if (distance > TheMaxDistance)
						TheMaxDistance = distance;
				}
			}
		}
	}
	return TheMaxDistance;
}

// VALID

double distanceClusters(Cluster* c1, Cluster* c2)
{
	double dist = 0;
	dist = sqrt(pow((c1->x - c2->x), 2.0) + pow((c1->y - c2->y), 2.0) + pow((c1->z - c2->z), 2.0));
	return dist;
}

// VALID

double distancePoints(Point* point1, Point* point2)
{
	double dist = 0;
	dist = sqrt(pow((point1->x - point2->x), 2.0) + pow((point1->y - point2->y), 2.0) + pow((point1->z - point2->z), 2.0));
	return dist;
}

// VALID

/* Returns the distance between two points - in Square */
double euclid_dist_2(Point p, Cluster c) {
	double euclid_dist_square = 0.0;
	euclid_dist_square = pow((p.x - c.x), 2.0) + pow((p.y - c.y), 2.0) + pow((p.z - c.z), 2.0);
	return euclid_dist_square;
}

// VALID

/* Returns index of nearest cluster centre point to the the given point parameter. This is called by parallel_kmeans() done locally for each process.*/
int assign_point_membership(Cluster* clusters, int numClusters, int *myid, Point point)
{

	int   index = 0, i;
	double dist = 0, min_dist = 0;

	/* find the cluster id that has min distance to object */
	index = clusters[0].id;
	/*Calculate the minimum distance between single point and each Cluster's Middle-Centre Point*/
	min_dist = euclid_dist_2(point, clusters[0]);

	for (i = 1; i < numClusters; i++) {
		dist = euclid_dist_2(point, clusters[i]);
		/* no need square root */
		if (dist < min_dist) { /* find the min and its array index */
			min_dist = dist;
			index = clusters[i].id;
		}
	}

	return(index);
}



void pointsLocationRemainderOMP(Point* pointsArray, double theTime, int remainderOffset, int totalPointsSize)
{
	int i;
#pragma omp parallel for
	for (i = remainderOffset; i < totalPointsSize; i++)
	{
		pointsArray[i].x = pointsArray[i].x + theTime*pointsArray[i].Vx;
		pointsArray[i].y = pointsArray[i].y + theTime*pointsArray[i].Vy;
		pointsArray[i].z = pointsArray[i].z + theTime*pointsArray[i].Vz;
	}
}