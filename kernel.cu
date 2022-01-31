#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include "StructsHeader.h"


__global__ void movePoints(Point* devPoints, int allPointsSize, unsigned int numofThreadsperBlock, unsigned int numofBlocks, double theTime)
{
	unsigned int index = blockIdx.x * numofThreadsperBlock + threadIdx.x;
	unsigned int threadwork = allPointsSize / ((numofThreadsperBlock)*(numofBlocks));

	for (unsigned int i = index*threadwork; i < (index* threadwork) + threadwork; i++)
	{
		
			devPoints[i].x = devPoints[i].x + theTime*devPoints[i].Vx;
			devPoints[i].y = devPoints[i].y + theTime*devPoints[i].Vy;
			devPoints[i].z = devPoints[i].z + theTime*devPoints[i].Vz;
		
	}
}

//This function is the Bridge between Host Code & GPU.

cudaError_t pointsLocationWholeCuda(int allPointsSize, double theTime, Point* pointsArray, unsigned int threadsperblock, int* remainderOffset)
{

	Point* devPoints = NULL;
	unsigned int numofBlocks = allPointsSize / threadsperblock;
	cudaError_t cudaStatus;

	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		cudaFree(devPoints);
	}
	// Allocated devPoints inside GPU
	cudaStatus = cudaMalloc((void**)&devPoints, allPointsSize * sizeof(Point));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		cudaFree(devPoints);

	}
	// Copy input vector (pointsArray) from host memory (CPU) to GPU's buffer devPoints.
	cudaStatus = cudaMemcpy(devPoints, pointsArray, allPointsSize * sizeof(Point), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		cudaFree(devPoints);
	}

	// Launch a kernel on the GPU
	movePoints << <numofBlocks, threadsperblock >> > (devPoints, allPointsSize, threadsperblock, numofBlocks, theTime);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "movePointsLocation launch failed: %s\n", cudaGetErrorString(cudaStatus));
		cudaFree(devPoints);

	}

	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
		cudaFree(devPoints);
	}

	// Copy output vector from GPU buffer to host memory.
	cudaStatus = cudaMemcpy(pointsArray, devPoints, allPointsSize * sizeof(Point), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		cudaFree(devPoints);
	}

 
	cudaFree(devPoints);
	*remainderOffset = numofBlocks*threadsperblock;
	return cudaStatus;
}


