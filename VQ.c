#include "VQ.h"
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>


#define N_CENTROIDS ((int) 5)    //
#define N_CODES		((int) 128)  // Needs to be empircally found
#define MIN_ERROR ((float) 2.0 ) // Needs to be empircally found
#define MAX_REP  ((int) 20)     // Couple with min_error or standalone 
#define N_MFCC ((int)13)

typedef struct {
	float values[N_MFCC];
	float dist_cent[N_CENTROIDS]; //distance from each centeroid
	int assign_cent;
}key_vec;

float euclid_dist2(float* vec1, float* vec2, int length)
{
	int i;
	float error = 0;
	for (i = 0; i < length; i++)
		error += fabsf(powf(vec1[i] - vec2[i], 2));

	return error;
}

/*
float euclid_dist(float* vec1, float* vec2, int length) //No point of using as is slower
{
	return sqrtf(euclid_dist2(vec1, vec2, length));
}
*/

void k_means_plus_init(float* centroids, key_vec* vec_list, int nVectors)
{
	int i, j;
	time_t t;

	int randVec;
	float randVal;
	int idx;

	float* probabilities = (float*)malloc(sizeof(float) * nVectors);
	float total_dist2;

	//Choose one data point randomly for first center

	//srand((unsigned)time(&t));
	randVec = rand() % (nVectors + 1);
	for (i = 0; i < N_MFCC; i++)
		centroids[i] = vec_list[randVec].values[i];

	//Finding the rest of the initial centroids
	for (i = 1; i < N_CENTROIDS; i++) {

		//Calculate distances away from centroid
		total_dist2 = 0;
		for (j = 0; j < nVectors; j++) {
			probabilities[j] = euclid_dist2(centroids + (i - 1) * N_MFCC, vec_list[j].values, N_MFCC);
			total_dist2 += probabilities[j];
		}
		//Make into pdf function based on distances
		for (j = 0; j < nVectors; j++)
			probabilities[j] /= total_dist2;

		//Make it into a cdf function by integration
		for (j = 1; j < nVectors; j++)
			probabilities[j] += probabilities[j - 1];

		//Choose vector index based on probabilites 
		randVal = (float)rand() / (float)RAND_MAX;
		idx = nVectors - 1; //Failsafe in case of rounding over 1 for CDF
		for (j = 0; j < nVectors; j++) {
			if (randVal < probabilities[j]) {
				idx = j;
				break;
			}
		}

		//Assign new centroid
		for (j = 0; j < N_MFCC; j++)
			centroids[i * N_MFCC + j] = vec_list[idx].values[j];
	}

	free(probabilities);
}

void k_means_init(float* centroids, key_vec* vec_list, int nVectors)
{
	int i, j;
	time_t t;

	int rand_index;
	int chosen_vec;
	int lower;
	int upper;

	int* storage = (int*)malloc(sizeof(int) * nVectors);

	//Store indexes to be used in shuffle
	for (i = 0; i < nVectors; i++)
		storage[i] = i;

	//Random, nonrepeat sort by shuffling array
	//srand((unsigned)time(&t));
	for (i = 0; i < N_CENTROIDS - 1; i++) {
		lower = i;
		upper = nVectors - 1;
		rand_index = (rand() % (upper - lower + 1)) + lower;  // (rand() % (upper - lower + 1)) + lower

		chosen_vec = storage[rand_index];
		storage[rand_index] = storage[i];

		for (j = 0; j < N_MFCC; j++)
			centroids[i * N_MFCC + j] = vec_list[chosen_vec].values[j];
	}

	for (j = 0; j < N_MFCC; j++)
		centroids[(N_CENTROIDS - 1) * N_MFCC + j] = vec_list[storage[nVectors - 1]].values[j];


	free(storage);
}


void center(float* centroids, int* count, key_vec* vec_list, int nVectors)
{
	int i, j;

	//Clear data
	for (i = 0; i < N_CENTROIDS; i++)
		for (j = 0; j < N_MFCC; j++)
			centroids[i * N_MFCC + j] = 0;

	//Find the mean
	for (i = 0; i < nVectors; i++)
		for (j = 0; j < N_MFCC; j++)
			centroids[vec_list[i].assign_cent * N_MFCC + j] += vec_list[i].values[j];

	for (i = 0; i < N_CENTROIDS; i++)
		for (j = 0; j < N_MFCC; j++)
			centroids[i * N_MFCC + j] /= (float)count[i];

}

void initialze_list(float* in, key_vec* vec_list, int nVectors)
{
	int i, j;

	for (i = 0; i < nVectors; i++) {
		for (j = 0; j < N_MFCC; j++) {
			vec_list[i].values[j] = in[i * N_MFCC + j];
			vec_list[i].assign_cent = 0;
		}
		for (j = 0; j < N_CENTROIDS; j++)
			vec_list[i].dist_cent[j] = 0;
	}
}


void assign(float* centroids, key_vec* vec_list, int nVectors)
{
	int i, j;
	int prev, min;

	for (i = 0; i < nVectors; i++)
		for (j = 0; j < N_CENTROIDS; j++)
			vec_list[i].dist_cent[j] = euclid_dist2(vec_list[i].values, (centroids + N_MFCC * j), N_MFCC);

	for (i = 0; i < nVectors; i++) {
		vec_list[i].assign_cent = 0;
		min = vec_list[i].dist_cent[0];
		for (j = 1; j < N_CENTROIDS; j++) {
			if (vec_list[i].dist_cent[j] < min) {
				vec_list[i].assign_cent = j;
				min = vec_list[i].dist_cent[j];
			}
		}
	}
}

void count_centroids(int* groupCount, key_vec* vec_list, int nVectors)
{
	int i;
	for (i = 0; i < N_CENTROIDS; i++)
		groupCount[i] = 0;

	for (i = 0; i < nVectors; i++)
		groupCount[vec_list[i].assign_cent]++;
}

float error_calc(float* centroids, key_vec* vec_list, int nVectors)
{
	int i;
	int cent_group;
	float error = 0;

	for (i = 0; i < nVectors; i++) {
		cent_group = vec_list[i].assign_cent;
		error += vec_list[i].dist_cent[cent_group];
	}

	return error;
}
//copies vector 2 into vector 1
void copy(float* vec1, float* vec2, int length1, int length2)
{
	int i, j;
	for (i = 0; i < length1; i++)
		for (j = 0; j < length2; j++)
			vec1[length2 * i + j] = vec2[length2 * i + j];
}
bool zero_element(float* vec, int length)
{
	int i;
	bool zero_found = false;

	for (i = 0; i < length; i++)
		if (vec[i] == 0)
			zero_found = true;

	return zero_found;
}

//nInputs represent number input vectors
float k_means(float* input, int nVectors, float* code)
{
	int rep;
	float error, error_prev, error_best;
	bool change_occurred;
	bool bad_initialize;
	bool copied;

	key_vec* vec_list = (key_vec*)malloc(nVectors * sizeof(key_vec));

	float centroids[N_CENTROIDS * N_MFCC];
	int nGroupCount[N_CENTROIDS];

	//Initialization list of vectors
	initialze_list(input, vec_list, nVectors);

	rep = 0;
	error_best = 0;
	error = 0;
	error_prev = error;
	copied = false;

	while (rep < MAX_REP) {
		//Assign initial centroids
		k_means_init(centroids, vec_list, nVectors);
		change_occurred = true;
		bad_initialize = false;
		while (change_occurred) {
			//Assign vectors to closest centroid, return true if vector reassignment occurred
			assign(centroids, vec_list, nVectors);

			//Count vectors within each group separately 
			count_centroids(nGroupCount, vec_list, nVectors);
			bad_initialize = zero_element(nGroupCount, N_CENTROIDS);
			if (bad_initialize)
				break;

			//Make centroid equal its group mean vector
			center(centroids, nGroupCount, vec_list, nVectors);

			error = error_calc(centroids, vec_list, nVectors);
			change_occurred = error != error_prev;
			error_prev = error;
		}


		if ((error < error_best || !copied) && !bad_initialize) {
			copied = true;
			copy(code, centroids, N_CENTROIDS, N_MFCC);
			error_best = error;
		}
		rep++;
	}

	free(vec_list);

	return error_best;
}
