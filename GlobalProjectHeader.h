#ifndef GLOBALPROJECTHEADER_H_
#define GLOBALPROJECTHEADER_H_


/* --------------------------- Includes ---------------------------*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include <omp.h>
#include "FunctionPrototypesHeader.h"
#include "StructsHeader.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda.h"

/* --------------------------- Constants ---------------------------*/

#define _CRT_SECURE_NO_WARNINGS
#define MASTER 0
#define MPI_FLAG 0
#define KMEANS_TERMINATE_FLAG 1
#define CLOCK_TERMINATE_FLAG 2

#define FILE_INPUT "D:\\ParallelProgramming_FinalProject_315441592\\ParallelProgramming_FinalProject_315441592\\data2100.txt"
#define FILE_OUT "D:\\ParallelProgramming_FinalProject_315441592\\ParallelProgramming_FinalProject_315441592\\outputdata2100.txt"
// 10000 42 5.000000 1.000000 50 1.000000
#define NUMPOINTS 10000
#define NUMCLUSTERS 42
#define T 30.000
#define DT 0.100 
#define LIMIT 200
#define QUALITY_MEASURE  0.100 

#define NUMOFCOORDINATES 3
#define NUM_DIMENSIONS 3
#define TRUE 1
#define FALSE 0
#define THREADS_PER_BLOCK 1024
 

#endif