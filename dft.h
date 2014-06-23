#ifndef DFT_H
#define DFT_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#define PI 3.14159



int complex_mult_mat(int size, double* input_1 , double* input_2, double* output)
{
	int i = 0;
	for(; i < size; i++)
	{
		output[i << 1] = input_1[i << 1] * input_2[i<<1] - input_1[(i << 1) +1]*input_2[(i << 1) +1];
		output[(i << 1) +1] = input_1[i << 1] * input_2[(i << 1) +1] + input_1[( i << 1) +1]*input_2[(i<<1)];
		
	}
	return 1;
}

int complex_mult_vec(int size, double* input_1, double* input_2, double* Re, double* Im)
{
	   int i = 0;
        for(; i < size; i++)
        {
                Re[i] = input_1[i << 1] * input_2[i<<1] - input_1[(i << 1) +1]*input_2[(i << 1) +1];
                Im[i] = input_1[i << 1] * input_2[(i << 1) +1] + input_1[( i << 1) +1]*input_2[(i<<1)];
        }
        return 1;
}


//just works for real signals !!!! (input is just a one dimensional vector
int dft(int size,double* input,double* output)
{
	int i = 0;
	int x = 0;
	int computed = 0;
	int b = 0;
	int j = 0;	
	double* coef_cos = (double*) malloc(sizeof(double) * (size*size + size + 1));
	double* coef_sin = (double*) malloc(sizeof(double) * (size*size + size + 1));
	

	double turn  = 0;
	double turn_perc = 0;
	//finding the coeffizients
	for(;i < size; i++)
	{
		
		//Kernel Vorbereitung -> später werden diese Werte direkt in den Kernel geschrieben! 
		for(x=0;x < size; x++)
		{
			computed = 0;
			turn = (((double)i*(double)x/(double)size));
			turn_perc = turn - (int) turn; 

			if(turn_perc == 0 || turn_perc == 1)
			{
				coef_cos[i*size +x] = 1;				
				coef_sin[i*size +x] = 0;
				computed++;
			}
			if(turn_perc == 0.25)
			{
				coef_cos[i*size +x] = 0;
				coef_sin[i*size +x] = -1;
				computed++;

			}
			if(turn_perc == 0.5)
			{
				coef_cos[i*size +x] = -1;
				coef_sin[i*size +x] = 0;
				computed ++;
			}
			if(turn_perc == 0.75)
			{
				coef_cos[i*size +x] = 0;
				coef_sin[i*size +x] = 1;
				computed++;
			}

			if(computed == 0)
			{
				printf("this takes time \n");
				coef_cos[i*size +x] =  cos(turn*2*M_PI);
				coef_sin[i*size +x] = sin(turn*2*M_PI);
				computed++;
			}
			if(computed > 1)
			{
				
				printf("Something went horrible wrong in the dft \n");
				return -1;
			}
				
		}
		j++;
	}

	//The real transformation
	for(i=0; i < size; i++)
	{
	for(x=0; x < size; x++)
		{
			output[i << 1] += input[x]*coef_cos[i*size +x];
			output[(i << 1) +1] += input[x]*(coef_sin[i*size +x]*(-1));

		}
	}
	return 1;
} 

int idft(int size,double* input,double* output)
{
        int i = 0;
        int x = 0;
        int computed = 0;
        int b = 0;
        int j = 0;
        double* coef_cos = (double*) malloc(sizeof(double) * (size*size + size + 1));
        double* coef_sin = (double*) malloc(sizeof(double) * (size*size + size + 1));


        double turn = 0;
        double turn_perc = 0;
        //finding the coeffizients
        for(;i < size; i++)
        {

                //Kernel Vorbereitung -> später werden diese Werte direkt in den Kernel geschrieben!
                for(x=0;x < size; x++)
                {
                        computed = 0;
                        turn = (((double)i*(double)x/(double)size));
                        turn_perc = turn - (int) turn;

                        if(turn_perc == 0 || turn_perc == 1)
                        {
                                coef_cos[i*size +x] = 1;
                                coef_sin[i*size +x] = 0;
                                computed++;
                        }
                        if(turn_perc == 0.25)
                        {
                                coef_cos[i*size +x] = 0;
                                coef_sin[i*size +x] = -1;
                                computed++;

                        }
                        if(turn_perc == 0.5)
                        {
                                coef_cos[i*size +x] = -1;
                                coef_sin[i*size +x] = 0;
                                computed ++;
                        }
                        if(turn_perc == 0.75)
                        {
                                coef_cos[i*size +x] = 0;
                                coef_sin[i*size +x] = 1;
                                computed++;
                        }

                        if(computed == 0)
                        {
                                coef_cos[i*size +x] = cos(turn*2*M_PI);
                                coef_sin[i*size +x] = sin(turn*2*M_PI);
                                computed++;
                        }
                         if(computed > 1)
                        {

                                printf("Something went horrible wrong in the idft \n");
                                return -1;
                        }

                }
                j++;


        }



	printf("size = %d", size);



        //The real transformation
        for(i=0; i < size; i++)
        {
        for(x=0; x < size; x++)
                {

                        output[i << 1] += (input[x << 1]*coef_cos[i*size +x] - input[(x << 1) +1]*coef_sin[i*size +x])/(double)size;
                        output[(i << 1) +1] += (input[x << 1]*(coef_sin[i*size +x]) + input[(x << 1) + 1] * coef_cos[i*size +x])/(double)size;


                }
        }
        return 1;
}




//just works when both vectors have the same length! 
int convolute(int size, double* input,double* input_2, double* output)
{
	int i = 0;
        int x = 0;
        int computed = 0;
        int b = 0;
        int j = 0;
	
	//fill up vectors
	double* input_filled = (double*) malloc(sizeof(double)*((4*size) +1));
	double* input_2_filled = (double*) malloc(sizeof(double)*((4*size) +1));
	
	double turn = 0;
        double turn_perc = 0;

	for(i=0; i < size; i++)
	{
		printf("input[%i] = %f \n",i, input[i]);
	}


	//making a new vector
	for(i=0; i < size; i++)
	{
		input_filled[i] = input[i];
		input_2_filled[i] = input_2[i];
	}

	//Zero-Adding to get rid off zyclic
	for(i=0; i < 3*size; i ++)
	{
		input_filled[i + size] = 0;
                input_2_filled[i + size] = 0;
	}
	size = 3*size +1;
	for(i = 0; i < size; i++)
	{
		printf("input_1 now[%i] = %f \t input_2 now[%i] = %f \n", i, input_filled[i], i , input_2_filled[i]);
	}

	double* coef_cos = (double*) malloc(sizeof(double) * ((size*size + size) +2));
	double* coef_sin = (double*) malloc(sizeof(double) * ((size*size + size) +2));

	double* trans_one = (double*) malloc(sizeof(double) * ((size << 2) + 2));
        double* trans_two = (double*) malloc(sizeof(double)*((size << 1) +2));

        //convolution result
        double* conv = (double*) malloc(sizeof(double)*((size << 2) +2));

	//finding the coeffizients
        for(i=0;i < size; i++)
        {

                //Kernel Vorbereitung -> später werden diese Werte direkt in den Kernel geschrieben!
                for(x=0;x < size; x++)
                {
                        computed = 0;
                        turn = (((double)i*(double)x/(double)size));
                        turn_perc = turn - (int) turn;

                        if(turn_perc == 0 || turn_perc == 1)
                        {
                                coef_cos[i*size +x] = 1;
                                coef_sin[i*size +x] = 0;
                                computed++;
			}
                        if(turn_perc == 0.25)
                        {
                                coef_cos[i*size +x] = 0;
                                coef_sin[i*size +x] = -1;
                                computed++;

                        }
                        if(turn_perc == 0.5)
                        {
                                coef_cos[i*size +x] = -1;
                                coef_sin[i*size +x] = 0;
                                computed ++;
                        }
                        if(turn_perc == 0.75)
                        {
                                coef_cos[i*size +x] = 0;
                                coef_sin[i*size +x] = 1;
                                computed++;
                        }

                        if(computed == 0)
                        {
                                coef_cos[i*size +x] = cos(turn*2*M_PI);
                                coef_sin[i*size +x] = sin(turn*2*M_PI);
                                computed++;
                        }
                         if(computed > 1)
                        {

                                printf("Something went horrible wrong in the idft \n");
                                return -1;
                        }

                }
                j++;


        }


		//Transformation		
	        for(i=0; i < size; i++)
        	{
        		for(x=0; x < size; x++)
                	{
                        	trans_one[i << 1] += input_filled[x]*coef_cos[i*size +x];
                        	trans_one[(i << 1) +1] += input_filled[x]*(coef_sin[i*size +x]*(-1));
                		trans_two[i << 1] += input_2_filled[x]*coef_cos[i*size +x];
				trans_two[(i << 1) +1] += input_2_filled[x]*(coef_sin[i*size +x]*(-1));

			}
        	}


		//The multiplikation itself -> Convoltuon
		for(i=0; i < size; i++)
		{
			complex_mult_mat(size, trans_one, trans_two, conv);
			printf("Re[%i] = %f \t Im[%i] = %f \t Re[%i] = %f  \t Im[%i] = %f \n", i, trans_one[i << 1] , i , trans_one[( i << 1) +1], i , trans_two[i << 1], i , trans_two[(i << 1) +1]);  			


		}

		idft(size, conv, output);


}



#endif

