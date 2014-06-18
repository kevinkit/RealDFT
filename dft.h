#ifndef DFT_H
#define DFT_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#define PI 3.14159


int dft(int size,double* input,double* output)
{
	int i = 0;
	int x = 0;
	double a = 0;
	double sum = 0;
	double sum_2 = 0;
	int computed = 0;
	int b = 0;
	double c = 0;
	int j = 0;	
	double* coef_cos = (double*) malloc(sizeof(double) * (size*size + size + 1));
	double* coef_sin = (double*) malloc(sizeof(double) * (size*size + size + 1));
	
	for(;i < size; i++)
	{
		
		//Kernel Vorbereitung -> später werden diese Werte direkt in den Kernel geschrieben! 
		for(x=j;x < size; x++)
		{
			a = ((double)i*(double)x)/size;
			b = (int) (i*x) / size;
			c = a - (double) b;
			computed = 0;
			sum = 0;
			sum_2 = 0;

			if(c == 0 || c == 1)
			{
				coef_cos[i*size +x] = 1;				
				coef_sin[i*size +x] = 0;
				computed++;
			}
			if(c == 0.25)
			{
				coef_cos[i*size +x] = 0;
				coef_sin[i*size +x] = -1;
				computed++;

			}
			if(c == 0.5)
			{
				coef_cos[i*size +x] = -1;
				coef_sin[i*size +x] = 0;
				computed ++;
			}
			if(c == 0.75)
			{
				coef_cos[i*size +x] = 0;
				coef_sin[i*size +x] = 1;
				computed++;
			}

			if(computed == 0)
			{
				coef_cos[i*size +x] =  cos((2*PI*i*x)/size);
				coef_sin[i*size +x] = -sin((2*PI*i*x)/size);
				computed++;
			}
			if(computed > 1)
			{
				return -1;
			}
				
		}
		j++;
		

		for(x=i; x < size; x++)
		{
			coef_cos[x*size + i] = coef_cos[i*size + (x - i)];
			coef_sin[x*size + i] = coef_sin[i*size + (x - i)];
		}

	}


	for(i=0; i < size; i++)
	{
	
	for(x=0; x < size; x++)
		{
			output[i << 1] += input[i]*coef_cos[i*size +x];
			output[(i << 1) +1] += input[i]*coef_sin[i*size +x];
		}
	}
	return 1;
} 

#endif
