#include "FunctionPrototypesHeader.h"
#include "GlobalProjectHeader.h"
#include "StructsHeader.h"

void freeAll(Point* allpoints, Cluster* clustersArray, Point* pointsBufferForEachProcess) {
	free(allpoints);
	free(clustersArray);
	free(pointsBufferForEachProcess);
}


void checkDynamicAllocation(const void* ptr)
{
	if (!ptr)
	{
		printf("Dynamic Memory Allocation Error --> Not Enough Memory !");
		fflush(stdout);
		MPI_Finalize();
		exit(2);
	}
}
