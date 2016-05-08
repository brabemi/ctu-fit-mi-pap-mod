#define PARALLEL_OPENMP //< define for OpenMP, undefine for sequential
#include <omp.h> // OpenMP library

// if the following variable is NOT defined, program will not use any function nor include the CImg library
// this also means that other libraries are needed
//~ #define CIMG_VISUAL

//~ #define SSE_SQRT

#define LOGGING // all logging output except for the line with "#THREADS #SECONDS" info
#undef LOGGING

#ifdef CIMG_VISUAL
#include "CImg/CImg.h" // lib for visualisation
#endif
#include <cstdio> // printf
#include <cstdlib> // srand
#include <cmath> // sqrt
#include <xmmintrin.h> //SSERsqrt

// include our functions/structs
#include "generator/ioproc.h" // process input file
#include "generator/SimConfig.h" // Struct with simulation config/settings

// include CUDA libs
#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda_runtime_api.h>


#define __CUDA_INTERNAL_COMPILATION__
#include "math_functions.h" // rsqrt
#undef __CUDA_INTERNAL_COMPILATION__

// current mapping of simulation parameters - just to clarify what is what 
#define WIDTH sconf.width // // x-coordinate, width of simulation area
#define HEIGHT sconf.height // y-coordinate, height of simulation area
#define DEPTH sconf.depth // z-coordinate, depth of simulation area
#define MAX_SPEED sconf.max_speed // maximal default speed of particle, eg. 4 => speed in (-4;4)
#define MAX_WEIGHT sconf.max_weight // maximal mass of particle, eg. 4 => speed in (0;4)
// amount of particles
#define AMOUNT sconf.amount // amount of particles
// number of simulation steps is specified in SimConfig struct as simulation_steps
#define SQRT_MAX_WEIGHT (int)sqrt(sconf.max_weight) // constant for calculation particle color {mass/SQRT_MAX_WEIGHT, 255, mass%SQRT_MAX_WEIGHT}
// ------------- END OF MODIFIED MAPPINGs

// parameters for calculations (movement of particles)
#define F_QUOC 0.0005 // compensatory quotient for  MAX_WEIGHT
#define BOUNCE_LOSS 0.8 // conversion rate of velocity on bounce, eg. 0.8 => 80% of speed after bounce

/*
CURRENT STATE:
- finite simulation of NBody
- bouncing, borders
- no particle collisions
- no detection if 

LAST UPDATE (just to clarify what happened - move to CURRENT STATE after review):
- inner loop in simulation is now parallelized with parallel reduction of accelerations
- define maximum steps of simulation - viz SimConfig::simulation_steps
- input (definition of particles) from file/as parameters from comm. line - viz generator.cpp and ioproc.cpp
- input constants as parameters (WIDTH, MAX_SPEED, ...) - viz generator.cpp
- for CImg visualisation, variable "CIMG_VISUAL" must be defined - uncomment #define CIMG_VISUAL at the beginning of this file

- created standalone generator (in directory "generator")
- simulation should be run with: ./simulator input.txt
- valgrind says the program leaks a bit (6 block every time) - with --leak-check=full --show-leak-kinds=all, it seems that libgomp causes that

- reworked makefile (make compile, make clean)

TODO:
...

*/

#ifdef CIMG_VISUAL
using namespace cimg_library; // -> no need to use cimg_library::function()
#else
using namespace std;
#endif

__device__ __host__ bool bounce(float x, float y, float z, const float maxX, const float maxY, const float maxZ) {
	return (x < 0) || (maxX < x) || (y < 0) || (maxY < y) || (z < 0) || (maxZ < z);
}

__device__ __host__ float debounce_vel(float vel, float pos, int min, int max) {
	if((pos < min) || (max < pos)) {
		return -1*BOUNCE_LOSS*vel;
	}
	return BOUNCE_LOSS*vel;
}

__device__ __host__ float debounce_pos (float pos, int min, int max) {
	if(pos < min) {
		return -1 * pos;
	}
	if(pos > max) {
		return max - (pos - max);
	}
	return pos;
}

static void HandleError(cudaError_t err, const char * file, int line){
	if (err != cudaSuccess) {
		printf("%s in %s at line %d\n", cudaGetErrorString(err), file, line);
		exit(EXIT_FAILURE);
	}
}
#define HANDLE_ERROR( err ) (HandleError(err, __FILE__, __LINE__))

__global__ void CopyCoordinatesKernel(float4 * sourceCoords, float4 * newCoords, int n)
{
	int gtid = blockIdx.x * blockDim.x + threadIdx.x;
	if (gtid >= n) return;

	float4 source = sourceCoords[gtid];
	float4 newCoordsVec = newCoords[gtid];

	source.x = newCoordsVec.x;
	source.y = newCoordsVec.y;
	source.z = newCoordsVec.z;
	// source.m is still the same

	sourceCoords[gtid] = source;
}

__global__ void CopyCoordinatesKernelSeq(float4 * sourceCoords, float4 * newCoords, int n)
{
	for (int i = 0; i < n; i++)
	{
		float4 source = sourceCoords[i];
		float4 newCoordsVec = newCoords[i];

		source.x = newCoordsVec.x;
		source.y = newCoordsVec.y;
		source.z = newCoordsVec.z;
		// source.m is still the same

		sourceCoords[i] = source;
	}
}

__device__ float3 perParticleAcceleration(float4 first, float4 second, float3 aXYZ, float eps)
{
	float3 dXYZ;
	dXYZ.x = first.x - second.x;
	dXYZ.y = first.y - second.y;
	dXYZ.z = first.z - second.z;

	float tmp_sum = dXYZ.x * dXYZ.x + dXYZ.y * dXYZ.y + dXYZ.z * dXYZ.z + eps;
	
	//float invr = 1/sqrtf(tmp_sum);
	float invr = rsqrtf(tmp_sum);

	float invr3 = invr*invr*invr;

	float f = F_QUOC * second.w * invr3;

	aXYZ.x += dXYZ.x * f;
	aXYZ.y += dXYZ.y * f;
	aXYZ.z += dXYZ.z * f;
	return aXYZ;
}

#define OBOC // which kind of work distribution is selected

__global__ void OneStepSimulation(float4 * sourceCoords, float4 * newCoords,
	float4 * velocities, float eps, float dt, 
	int n, int offset,
	const float maxX, const float maxY, const float maxZ)
{
	int globalThreadIndex = blockIdx.x * blockDim.x + threadIdx.x;
	if (globalThreadIndex >= n) return;

	// shared memory - used to cache particles from global memory
	extern __shared__ float4 particlesInSM[];
	// pointer to array of particles in global memory 
	float4 * particlesGlobal = sourceCoords;
	// iterators
	int i = 0; // iterator for particles, (for each particle ... (from 1 .. n))
	int blok = 0; // iterator for current block
	float3 aXYZ = { 0.0f, 0.0f, 0.0f }; // acceleration, to be counted
	
	float4 threadVector = particlesGlobal[globalThreadIndex]; // vector of particle assigned to this thread

	// control output
	//if (globalThreadIndex == 1) printf("%f %f %f %f\n", threadVector.x, threadVector.y, threadVector.z, threadVector.w);

	// pre-declaration of variables
	float3 dXYZ; // only assigned to, should be OK
	float tmp_sum; // only assigned to, should be OK
	float invr, invr3; // only assigned to, should be OK
	float f; // only assigned to, should be OK

	int maxParticles = blockDim.x; // amount of particles % threadsPerBlock

#ifdef OBOC // one block per one processor
	if (blockIdx.x == gridDim.x) { maxParticles = n % (blockDim.x+1); /*if (maxParticles == 0) maxParticles == blockDim.x;*/ }
#endif

	// for each particle ...
	for (i = 0; i < n; i += offset) { // offset == threadsPerBlock specified in calling of the kernel

		int indexOfParticle = blok * blockDim.x + threadIdx.x; // index of particle to be stored in SM

		particlesInSM[threadIdx.x] = particlesGlobal[indexOfParticle]; // copy from global to shared memory, "cache"

		__syncthreads(); // wait for every thread so the SM is full

		// count current subblock
		// original condition: j < blockDim.x
		for (int j = 0; j < maxParticles; j++) { // !!! this might cause problems, will work ONLY if every thread in this block copied particle to SM !!!
			
			dXYZ.x = threadVector.x - particlesInSM[j].x;
			dXYZ.y = threadVector.y - particlesInSM[j].y;
			dXYZ.z = threadVector.z - particlesInSM[j].z;

			tmp_sum = dXYZ.x * dXYZ.x + dXYZ.y * dXYZ.y + dXYZ.z * dXYZ.z + eps;

			//float invr = 1/sqrtf(tmp_sum);
			invr = rsqrtf(tmp_sum);

			invr3 = invr*invr*invr;

			f = F_QUOC * particlesInSM[j].w * invr3;

			aXYZ.x += dXYZ.x * f;
			aXYZ.y += dXYZ.y * f;
			aXYZ.z += dXYZ.z * f;
			//aXYZ = perParticleAcceleration(threadVector, particlesInSM[j], aXYZ, eps);
		}
		
		__syncthreads(); // wait for every thread so the SM is free to be edited
		blok++;
	}
	// load velocity from global memory
	float4 vel = velocities[globalThreadIndex];
	// update velocity
	vel.x += dt*aXYZ.x; /* update velocity of particle "i" */
	vel.y += dt*aXYZ.y;
	vel.z += dt*aXYZ.z;
	// no change to .w
	
	// update position of vector
	float4 newVec = threadVector;
	newVec.x = threadVector.x + dt*vel.x + 0.5*dt*dt*aXYZ.x;
	newVec.y = threadVector.y + dt*vel.y + 0.5*dt*dt*aXYZ.y;
	newVec.z = threadVector.z + dt*vel.z + 0.5*dt*dt*aXYZ.z;
	// no change to .w == mass of particle

	// check bouncing
	if (bounce(newVec.x, newVec.y, newVec.z, maxX, maxY, maxZ)) {
		// update of particle velocity, change direction and value (BOUNCE_LOSS)
		vel.x = debounce_vel(vel.x, newVec.x, 0, maxX);
		vel.y = debounce_vel(vel.y, newVec.y, 0, maxY);
		vel.z = debounce_vel(vel.z, newVec.z, 0, maxZ);
		// update of particle position
		newVec.x = debounce_pos(newVec.x, 0, maxX);
		newVec.y = debounce_pos(newVec.y, 0, maxY);
		newVec.z = debounce_pos(newVec.z, 0, maxZ);
	}

	// store new velocities in global memory
	velocities[globalThreadIndex] = vel;

	// store updated position in global memory - array "(x/y/z)new"
	newCoords[globalThreadIndex] = newVec;
}

__global__ void GPUPrintParticles(float4 * particles, int n)
{
	for (int i = 0; i < n; i++)
	{
		float4 myPosition = particles[i];
		printf("%d. %f %f %f %f\n", i, myPosition.x, myPosition.y, myPosition.z, myPosition.w);
	}
}

int main(int argc, char** argv) {
	
	// simulation configuration
	SimConfig sconf;
	int threads;

	#ifdef CIMG_VISUAL
	bool graphics;
	#endif

	// process input
	if(argc == 4)
	{
		threads = atoi(argv[1]);
		#ifdef CIMG_VISUAL
		graphics = atoi(argv[2]) == 1;
		#endif
		#ifdef LOGGING
		printf("Processing input file %s\n", argv[3]);
		#endif
		processInputFile(argv[3], sconf);
	}
	else
	{
		printf("ERROR: Wrong/Missing parameters.\n");
		printf("Expect: \t%s THREADS GRAPHICS[0,1] INPUT_FILE\n", argv[0]);
		return 0;
	}
	#ifdef LOGGING
	printf("Input file processed.\n");	
	printf("Particles: %d\tSteps: %d\n", sconf.amount, sconf.simulation_steps);
	#endif
	
	// iterators
	int i, j;
	
	// pointers to arrays
	float * x = sconf.x;
	float * y = sconf.y;
	float * z = sconf.z;
	
	float * m = sconf.m;
	
	float * vx = sconf.vx;
	float * vy = sconf.vy;
	float * vz = sconf.vz;
#ifdef WITH_CPU
	float * xnew = new float[sconf.amount];
	float * ynew = new float[sconf.amount];
	float * znew = new float[sconf.amount];
#endif
	// create an array that is in form suitable for float4 == 1 particle
	float * hostParticles = new float[sconf.amount * 4]; // x y z m
	float * hostVelocities = new float[sconf.amount * 4]; // vx vy vz 0
	for (i = 0; i < (sconf.amount); i++)
	{
		int base = i * 4;
		hostParticles[base + 0] = x[i];
		hostParticles[base + 1] = y[i];
		hostParticles[base + 2] = z[i];
		hostParticles[base + 3] = m[i];

		hostVelocities[base + 0] = vx[i];
		hostVelocities[base + 1] = vy[i];
		hostVelocities[base + 2] = vz[i];
		hostVelocities[base + 3] = 0.0f;
	}
	// initialize arrays on GPU
	float * devParticles = NULL, * devParticlesNew = NULL;
	float * devVelocities = NULL;
	HANDLE_ERROR(cudaMalloc((void**)&devParticles, (4 * sconf.amount * sizeof(float))));
	HANDLE_ERROR(cudaMalloc((void**)&devParticlesNew, (4 * sconf.amount * sizeof(float))));
	HANDLE_ERROR(cudaMalloc((void**)&devVelocities, (4 * sconf.amount * sizeof(float))));
	// copy content
	HANDLE_ERROR(cudaMemcpy(devParticles, hostParticles, (4 * sconf.amount), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(devVelocities, hostVelocities, (4 * sconf.amount), cudaMemcpyHostToDevice));

#ifdef WITH_CPU
	// variables used in computations
	float ax, ay, az;
	float dx, dy, dz;
	float invr, invr3;
	float f;
#endif
	// definitino of "n" - used in algorithm on site
	int n = sconf.amount;
	
	// constants for computing particle movement 
	float dt = 0.1f; // original value was 0.0001, that was too little for current values of particle parameters 
	float eps = 0.005f;
	
	#ifdef CIMG_VISUAL
	// image ~ "drawing panel"
	CImg<unsigned char> img;
	CImgDisplay main_disp;
	// colours
	const unsigned char red[] = { 255,0,0 }, green[] = { 0,255,0 }, blue[] = { 0,0,255 };

	if(graphics) {
		// image ~ "drawing panel"
		img = CImg<unsigned char> ( sconf.width, sconf.height,1,3);

		// initialization of window
		img.fill(0); //< fill img with black colour
	
		// draw all points
		for(i = 0; i < n; i++) {
			const unsigned char color[] = {(unsigned char) (m[i]/SQRT_MAX_WEIGHT), 255, (unsigned char) (((int) m[i])%SQRT_MAX_WEIGHT)};
			img.draw_point(x[i],y[i],color);
		}
		// create a Window (caption Playground) and fill it with image	
		main_disp = CImgDisplay(img,"Playground");
	}
	#endif

	unsigned steps = 0;

	// count amount of blocks and threads needed
	cudaDeviceProp prop;
	HANDLE_ERROR(cudaGetDeviceProperties(&prop, 0)); // get device properties

	int threadsPerBlock = 0, numOfBlocks = 0;


#ifdef OBOC // one block per one processor
	/* ----------------------------------------------------- */
	/* One block per processor */
	numOfBlocks = prop.multiProcessorCount;
	threadsPerBlock = ((sconf.amount) / numOfBlocks);
	if (threadsPerBlock > prop.maxThreadsPerBlock)
	{
		printf("NOPE NOPE NOPE");
	}
	/* ----------------------------------------------------- */
#endif

#ifdef WARP // one warp per one block - very slow, better to give more warps
	/* ----------------------------------------------------- */
	/* One block per processor */
	threadsPerBlock = prop.warpSize*4; // MUST be power of two
	numOfBlocks = ((sconf.amount - 1) / threadsPerBlock);
	numOfBlocks += 1;
	/* ----------------------------------------------------- */
#endif

#ifdef DEFAULT
	/* ----------------------------------------------------- */
	
	// One possible configuration of blocks and threads
	threadsPerBlock = (prop.maxThreadsPerBlock / 2); // MUST be power of 2
	// why only half of max - SM has limited amount of registers etc.
	// if limits reached -> kernel will not start (HANDLE_ERROR(startkernel<<<>>>()) will detect and print out this error)
	// /2 is just a working guess, may be optimized

	numOfBlocks = ((sconf.amount-1) / threadsPerBlock);
	numOfBlocks += 1; // amount = 1024, maxThreadsPerBlock = 1024 -> 2 blocks, both only half of max threads
	
	/* ----------------------------------------------------- */
#endif
	printf("----------------------------\n");
	printf("Amount of blocks: %d\n", numOfBlocks);
	printf("Amount of threads per block: %d\n", threadsPerBlock);

	// CUDA time measurement with events
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	
	cudaDeviceSynchronize(); // just in case ...

	// Start CUDA recording
	cudaEventRecord(start); 
	// Start OpenMP time measurement
	double t1 = omp_get_wtime();

//#define WITH_CPU

	#ifdef LOGGING
	printf("Starting simulation ...\n");
	#endif
	while (steps < sconf.simulation_steps) {

		// GPU simulation step
		OneStepSimulation << < numOfBlocks, threadsPerBlock, (threadsPerBlock*sizeof(float4)) >> > ((float4 *)devParticles, (float4 *)devParticlesNew,
			(float4 *)devVelocities, eps, dt, n, threadsPerBlock,
			sconf.width, sconf.height, sconf.depth);
		HANDLE_ERROR(cudaPeekAtLastError());
		cudaDeviceSynchronize(); // wait for Kernel to finish

		#ifdef WITH_CPU
		// CPU: compute new coordinates of all particles in parallel
		#ifdef PARALLEL_OPENMP
		#pragma omp parallel for num_threads(threads) private(ax,ay,az,dx,dy,dz,invr,invr3,f)
		#endif
		for (i = 0; i < n; i++) { /* Foreach particle "i" ... */
			ax = 0.0;
			ay = 0.0;
			az = 0.0;

			for (j = 0; j < n; j++) { /* Loop over all particles "j" */
				dx = x[j] - x[i];
				dy = y[j] - y[i];
				dz = z[j] - z[i];


				invr = 1.0 / sqrtf(dx*dx + dy*dy + dz*dz + eps);

				invr3 = invr*invr*invr;
				f = F_QUOC*m[j] * invr3;
				ax += f*dx; /* accumulate the acceleration from gravitational attraction */
				ay += f*dy;
				az += f*dz;
			}

			vx[i] += dt*ax; /* update velocity of particle "i" */
			vy[i] += dt*ay;
			vz[i] += dt*az;

			xnew[i] = x[i] + dt*vx[i] + 0.5*dt*dt*ax; /* update position of particle "i" */
			ynew[i] = y[i] + dt*vy[i] + 0.5*dt*dt*ay;
			znew[i] = z[i] + dt*vz[i] + 0.5*dt*dt*az;

			if (bounce(xnew[i], ynew[i], znew[i], sconf.width, sconf.height, sconf.depth)) {
				// update of particle velocity, change direction and value (BOUNCE_LOSS)
				vx[i] = debounce_vel(vx[i], xnew[i], 0, WIDTH);
				vy[i] = debounce_vel(vy[i], ynew[i], 0, HEIGHT);
				vz[i] = debounce_vel(vz[i], znew[i], 0, DEPTH);
				// update of particle position
				xnew[i] = debounce_pos(xnew[i], 0, WIDTH);
				ynew[i] = debounce_pos(ynew[i], 0, HEIGHT);
				znew[i] = debounce_pos(znew[i], 0, DEPTH);
			}
		}
		#endif




		#ifdef CIMG_VISUAL
		if(graphics) {
			if(main_disp.is_closed()) break;
			
			for(i = 0; i < n; i++) {
				const unsigned char color[] = {(unsigned char) (m[i]/SQRT_MAX_WEIGHT), 255, (unsigned char) (((int) m[i])%SQRT_MAX_WEIGHT)};
				img.draw_circle(x[i],y[i],1,color);
				//~ img.draw_circle(x[i],y[i],2,green);
			}
		}
		#endif

		// GPU copy
		//CopyCoordinatesKernelSeq << <1, 1 >> >((float4 *)devParticles, (float4 *)devParticlesNew, n);
		CopyCoordinatesKernel << < numOfBlocks, threadsPerBlock >> >
			((float4 *)devParticles, (float4 *)devParticlesNew, n);
		HANDLE_ERROR(cudaPeekAtLastError());
		cudaDeviceSynchronize(); // wait for kernel to finish
		
		// CPU copy
		#ifdef WITH_CPU
		#ifdef PARALLEL_OPENMP
		#pragma omp parallel for num_threads(threads)
		#endif
		for(i=0; i<n; i++) { /* copy updated positions back into original arrays */
		if( bounce(xnew[i], ynew[i], znew[i], sconf.width, sconf.height, sconf.depth) ) {
				#ifdef LOGGING
				printf("Particle %d out of borders (x, y, z) = (%0.3f, %0.3f, %0.3f)\n", i, xnew[i], ynew[i], znew[i]);
				#endif
			}
			x[i] = xnew[i];
			y[i] = ynew[i];
			z[i] = znew[i];
		}
		#endif

		#ifdef CIMG_VISUAL
		if(graphics) {
			// redraw the image and show it in the window
			if(steps%500 == 0) {
				img.fill(0); //< black background 
			}
	
			// draw all particles
			for(i = 0; i < n; i++) {
				//~ const unsigned char color[] = {(unsigned char) (m[i]/SQRT_MAX_WEIGHT), 255, (unsigned char) (((int) m[i])%SQRT_MAX_WEIGHT)};
				img.draw_circle(x[i],y[i],1,red);
				//~ img.draw_circle(x[i],y[i],2,green);
			}
			// display the image in window	
			img.display(main_disp);
			
			// wait for some time
			cimg::wait(20); // in milisec
		}
		#endif
		if(100*steps % sconf.simulation_steps == 0) {
			#ifdef LOGGING
			printf("%.1f%% completed\n", 100.0*steps/sconf.simulation_steps);
			#endif
		}
		steps++;
	}
	/*
	// copy end state of simulation from device to host	
	HANDLE_ERROR(cudaMemcpy(hostParticles, devParticles, (4 * sconf.amount), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaPeekAtLastError());

	cudaDeviceSynchronize();
	HANDLE_ERROR(cudaPeekAtLastError());

	// output directly from GPU
	GPUPrintParticles << <1, 1 >> >((float4 *)devParticles, 10);
	// output of values on CPU
	for (int i = 0; i < 10; i++)
	{
		printf("%d. %f %f %f %f\t", i, hostParticles[4*i], hostParticles[4*i + 1], hostParticles[4*i + 2], hostParticles[4*i + 3]);
		printf("%d. %f %f %f %f\n", i, x[i], y[i], z[i], m[i]);
	}*/
	/*
    // Control - differences are quite large in case of many steps (500 and more)
	double diffTol = 100.0;
	for (int i = 0; i < n; i++)
	{
		int base = i * 4;
		if (m[i] != hostParticles[base + 3]) printf("Masses of particle %d are different.\n", i);
		if (abs(x[i] - hostParticles[base + 0]) > diffTol) printf("X-coordinates of particle %d are too(?) different.\n", i);
		if (abs(y[i] - hostParticles[base + 1]) > diffTol) printf("Y-coordinates of particle %d are too(?) different.", i);
		if (abs(z[i] - hostParticles[base + 2]) > diffTol) printf("Z-coordinates of particle %d are too(?) different.", i);

	}*/
	
	double t2 = omp_get_wtime(); // in seconds

	// CUDA time measurement
	cudaEventRecord(stop); // stop recording

	cudaEventSynchronize(stop); // synchronized stop
	float miliseconds = 0; // init time 
	cudaEventElapsedTime(&miliseconds, start, stop); // count time

	#ifdef LOGGING
	printf("Time: %f seconds\n",(t2-t1));
	#endif

	// times should be almost equal (both measures almost the same)
	printf("OpenMP: %d %f\n", threads, (t2-t1));
	printf("CUDA: %f\n", miliseconds/(1000.0));
	
	// free CPU arrays
#ifdef WITH_CPU
	delete [] xnew;
	delete [] ynew;
	delete [] znew;
#endif

	delete[] hostParticles;
	delete[] hostVelocities;

	// free GPU arrays
	cudaFree(devParticles);
	cudaFree(devParticlesNew);
	cudaFree(devVelocities);

	return 0;
}

void recycleBin()
{
#ifdef ABCDEFGH
	// pointers to arrays on GPU
	float * devX, *devY, *devZ;
	float * devXnew, *devYnew, *devZnew;
	float * devM;
	float * devVx, *devVy, *devVz;
	// allocation of arrays for XYZ coordinates
	HANDLE_ERROR(cudaMalloc((void**)&devX, sconf.amount * sizeof(float)));
	HANDLE_ERROR(cudaMalloc((void**)&devY, sconf.amount * sizeof(float)));
	HANDLE_ERROR(cudaMalloc((void**)&devZ, sconf.amount * sizeof(float)));

	HANDLE_ERROR(cudaMalloc((void**)&devM, sconf.amount * sizeof(float)));

	HANDLE_ERROR(cudaMalloc((void**)&devXnew, sconf.amount * sizeof(float)));
	HANDLE_ERROR(cudaMalloc((void**)&devYnew, sconf.amount * sizeof(float)));
	HANDLE_ERROR(cudaMalloc((void**)&devZnew, sconf.amount * sizeof(float)));

	HANDLE_ERROR(cudaMalloc((void**)&devVx, sconf.amount * sizeof(float)));
	HANDLE_ERROR(cudaMalloc((void**)&devVy, sconf.amount * sizeof(float)));
	HANDLE_ERROR(cudaMalloc((void**)&devVz, sconf.amount * sizeof(float)));

	// copy CPU -> GPU
	// cudaMemcpy(to, from, amount, type)
	// coordinates
	HANDLE_ERROR(cudaMemcpy(devX, x, sconf.amount, cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(devY, y, sconf.amount, cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(devZ, z, sconf.amount, cudaMemcpyHostToDevice));
	// masses - maybe to constant/read-only memory?
	HANDLE_ERROR(cudaMemcpy(devM, m, sconf.amount, cudaMemcpyHostToDevice));
	// new coordinates
	HANDLE_ERROR(cudaMemcpy(devXnew, xnew, sconf.amount, cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(devYnew, ynew, sconf.amount, cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(devZnew, znew, sconf.amount, cudaMemcpyHostToDevice));
	// velocities
	HANDLE_ERROR(cudaMemcpy(devVx, vx, sconf.amount, cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(devVy, vy, sconf.amount, cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(devVz, vz, sconf.amount, cudaMemcpyHostToDevice));



	// free GPU arrays
	cudaFree(devX);
	cudaFree(devY);
	cudaFree(devZ);

	cudaFree(devM);

	cudaFree(devXnew);
	cudaFree(devYnew);
	cudaFree(devZnew);

	cudaFree(devVx);
	cudaFree(devVy);
	cudaFree(devVz);
#endif
}


void getDeviceInfo()
{
	int devCount = 0;
	cudaGetDeviceCount(&devCount);

	printf("Amount of devices: %d\n", devCount);

	for (int i = 0; i < devCount; i++) {
		cudaDeviceProp prop;
		HANDLE_ERROR(cudaGetDeviceProperties(&prop, i));
		/*if (cudaStatus != cudaSuccess)
		{
		printf("%s\n", cudaGetErrorString(cudaStatus));
		continue;
		}*/
		printf("Device Number: %d\n", i);
		printf("  Device name: %s\n", prop.name);

		printf("  Multiprocessor count: %d\n", prop.multiProcessorCount);

		printf("  Total global Memory: %d\n", prop.totalGlobalMem);

		printf("  Total const memory: %d\n", prop.totalConstMem);

		printf("  Shared memory per block: %d\n", prop.sharedMemPerBlock);

		printf("  Shared memory per multiprocessor: %d\n", prop.sharedMemPerMultiprocessor);

		printf("  Max threads per block: %d\n", prop.maxThreadsPerBlock);

		printf("  Max threads per multiprocessor: %d\n", prop.maxThreadsPerMultiProcessor);

		size_t size;
		cudaDeviceGetLimit(&size, cudaLimitMallocHeapSize);
		printf("Heap size limit: %d\n", size);

		printf("  Memory Clock Rate (KHz): %d\n",
			prop.memoryClockRate);
		printf("  Memory Bus Width (bits): %d\n",
			prop.memoryBusWidth);
		printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
			2.0*prop.memoryClockRate*(prop.memoryBusWidth / 8) / 1.0e6);
	}
}
