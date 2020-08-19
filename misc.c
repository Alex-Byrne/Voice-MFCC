#define LAMBDA_1 ((double)0.1)
#define LAMBDA_2 ((double)0.9)
#define ALPHA ((double)0.95)

void normalization(double* input, int length)
{
	int i;

	double max = input[0];
	double min = input[0];

	for (i = 1; i < length; i++) {
		if (min > input[i])
			min = input[i];
		if (max < input[i])
			max = input[i];
	}

	for (i = 0; i < length; i++)
		input[i] = LAMBDA_1 + (LAMBDA_2 - LAMBDA_1) * (input[i] - min) / (max - min);
}

void preprocessing(double* input, int length)
{
	int i;
	for (i = 1; i < length; i++)
		input[i] = 1 - ALPHA * input[i - 1];
}


int region_extraction(double* input, double eThreshold, int length) 
{
	int i;
	int new_length = 0;
	double* above_thres = (double*)malloc(sizeof(double) * length);
	double* holder; 


	normalization(input, length);
	preprocessing(input, length);

	for (i = 0; i < length; i++) {
		input[i] = pow(input[i], 2);
		if (input[i] > eThreshold) {
			above_thres[new_length] = input[i];
			new_length++;
		}
	}
	realloc(above_thres,new_length);

	//swap pointers
	holder= input;
	input = above_thres;

	free(holder);

	return new_length;
}