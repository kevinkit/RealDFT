#ifndef DFT_H
#define DFT_H


#include <stdlib.h>
#include <stdio.h>

#include <math.h>

#define PI 3.14159




//this function could be used to forbid memory corruption but I am not quite sure if works.........
double* memallec(double* array, int* size, int index,int size_index)
{
        int i = 0;

        FILE* file;

        file = fopen("errorlog.txt", "a+");


        if(index > size[size_index] - 1)
        {
		        fprintf(file, "Memory Corruption please check your boundaries, the size %d is too small \n in the index %d \n Tried to access %d ", size[size_index], size_index, index);

                printf("Something gone wrong \n");
                printf("I did that for you but you should check your memory allocation \n");
                printf("Wanted to access %d but size is only %d", size[size_index], index);
                printf("Changed the size to %d", 2*size[size_index] + index);

                double* swap = (double*) malloc(sizeof(double) * (size[size_index] +1));

                for(i=0; i < size[size_index]; i++)
                {
                        swap[i] = array[i];
                }

                double* newsth= (double*) malloc(sizeof(double) *( 2*size[size_index] +index + 1));



                for(i=0; i < size[size_index]; i++)
                {
                        newsth[i] = array[i];
                }
                free(array);
                size[0] = 2*size[size_index] + index;
                fclose(file);
                return newsth;

                }
        else
        {
                fclose(file);
                return array;
        }


}



//comlex multiplication of two vectors with the same length -> vector[i << 1] is Re and vector[(i << 1) + 1] is the Im
int complex_mult_mat(int size, double* input_1 , double* input_2, double* output)
{
	int i = 0;
	for(; i < size; i++)
	{
		output[i << 1] = input_1[i << 1] * input_2[i<<1] - input_1[(i << 1) +1]*input_2[(i << 1) +1];
		output[(i << 1) +1] = input_1[i << 1] * input_2[(i << 1) +1] + input_1[( i << 1) +1]*input_2[(i<<1)];
		if(((i<<1)+1) >= size*2 +1)
		{
			printf("!! SEGMENTATION FAULT IN complex_mult_mat !! \n");
			printf("wanted to access: %d \t But size is %d \n", (i << 1) +1, size);
			return -1;	

		}
	}
	return 1;
}

//complex multiplikation , the input is still a vector with vector[i << 1] = Re and vector[(i << 1) +1] is the Im but The output Re is the Real part and Im is the imaginary Output
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


//Used for adding zeros at the end of the vector
//BE CAREFUL THIS FUNCTION GOES OVER THE SIZE OF YOUR "NORMAL" ARRAY ! WHEN YOU ALLOCATE YOUR ARRAY MAKE SURE THAT YOU HAVE DOUBLE THE SIZE YOU THINK IT SHOULD HAVE
void zeroadding(int size, double* input, double* output)
{
	int i = 0;
	int size_array[1] = size;
	for(; i < size; i++)
	{
		output[i] = input[i];
		output[i + size] = 0;	
	}	
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
	free(coef_cos);
	free(coef_sin);
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

	}

	//The real transformation
	for(i=0; i < size; i++)
	{
		for(x=0; x < size; x++)
		{

			output[i << 1] += (input[x << 1]*coef_cos[i*size +x] - input[(x << 1) +1]*coef_sin[i*size +x])/(double)size;
			//		output[(i << 1) +1] += (input[x << 1]*(coef_sin[i*size +x]) + input[(x << 1) + 1] * coef_cos[i*size +x])/(double)size;

		}
	}
	

	free(coef_cos);
	free(coef_sin);

	return 1;

	

}




//just works when both vectors have the same length! 
int convolute(int size, double* input,double* input_2, double* output)
{

	printf("inside \n");
	int i = 0;
	int x = 0;
	int computed = 0;
	int b = 0;
	int j = 0;

	//fill up vectors
	double* input_filled = (double*) malloc(sizeof(double)*((2*size) +2));
	double* input_2_filled = (double*) malloc(sizeof(double)*((2*size) +2));
	printf("array \n");
	double turn = 0;
	double turn_perc = 0;

	zeroadding(size, input, input_filled); //does not leak;
	printf("zerroadding \n");
	zeroadding(size, input_2, input_2_filled); //does not leak;
	printf("zeroadding 2\n");


	size = 2*size +1;
	printf("size refreshed \n");

	
	int coef_size = (size)*(size) + size + 1;
	int trans_size = 2*(size) + 1;
	int conv_size = 2*(size) + 1;
	printf("after some sit \n");	
	int size_array[3];
	
	size_array[0] = coef_size;
	size_array[1] = trans_size;
	size_array[2] = conv_size;

	double* coef_cos = (double*) malloc(sizeof(double) * (coef_size));
	double* coef_sin = (double*) malloc(sizeof(double) * (coef_size));

	double* trans_one = (double*) malloc(sizeof(double) *(trans_size));
	double* trans_two = (double*) malloc(sizeof(double)* (trans_size));

	//convolution result
	double* conv = (double*) malloc(sizeof(double)*conv_size);

	//finding the coeffizients
	for(i=0;i < size; i++)
	{

		//Kernel Vorbereitung -> später werden diese Werte direkt in den Kernel geschrieben!
		for(x=0;x < size; x++)
		{
			printf("still on the run %d  \n", i);
			computed = 0;
			turn = (((double)i*(double)x/(double)size));
			turn_perc = turn - (int) turn;

			//			printf("size = %d \t allokiert = %d \tZugriff auf =  %d\n",size, size*size + size +2, i*size + x);


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
			if(computed > 1 || (i*size +x) >= coef_size)
			{

				printf("Something went horrible wrong in the idft \n");
				printf("Wanted to acces %d but size is only %d", i*size +x, size*size + size +2);
				return -1;
			}

		//	         coef_cos = memallec(coef_cos, size_array, i*size + x, 0);
                //		 coef_sin = memallec(coef_sin, size_array, i*size +x, 0);



		}

	}


	//Transformation		
	for(i=0; i < size; i++)
	{
		printf("on the run again %d\n", i);
		for(x=0; x < size; x++)
		{
			trans_one[i << 1] += input_filled[x]*coef_cos[i*size +x];
			trans_one[(i << 1) +1] += input_filled[x]*(coef_sin[i*size +x]*(-1));
			trans_two[i << 1] += input_2_filled[x]*coef_cos[i*size +x];
			trans_two[(i << 1) +1] += input_2_filled[x]*(coef_sin[i*size +x]*(-1));
			if((i << 1) +1 >= trans_size-1 || i*size +x >= coef_size-1)
			{
				printf("SEGMENATION FAULT IN CONVOLOTION \n");
				printf("wanted to access %d but size is only %d", (i << 1) +1, trans_size);
				printf("wanted to accces %d but size is only %d", i*size +x, coef_size);
			}

		
	
		}

	}

	//The multiplikation itself -> Convoltuon
	complex_mult_mat(size, trans_one, trans_two, conv);

	idft(size, conv, output);


	free(coef_cos);
        free(coef_sin);
        free(trans_one);
        free(trans_two);
        free(conv);
        free(input_filled);
        free(input_2_filled);

}

int convolute_td(int size, double* input , double* input_2, double* convolute)
{
	int i = 0;
	int x = 0;

	double* input_fill = (double*) malloc(sizeof(double) * (3*size +1));
	double* input_fill_2 = (double*) malloc(sizeof(double) * (3*size +1));

	zeroadding(size, input, input_fill);
	zeroadding(size, input_2, input_fill_2);

	size = 2*size +1;
	for(; i < size; i++)
	{
		for(x=0; x < i; x++)
		{
			convolute[i] += input_fill[i-x] * input_fill_2[size - x];		
		}	
	}
	free(input_fill);
	free(input_fill_2);
	return 1;

}

int findcofs(int size, double* coef_cos, double* coef_sin)
{
	int x = 0;
	int i = 0;
	int computed = 0;
	int b = 0;

	double turn = 0;
	double turn_perc = 0;

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


	}
	return 1;
}

int conv_cofs(int size,double* input, double* input_2,  double* coef_cos, double* coef_sin, double* output)
{

	double* trans_one = (double*) malloc(sizeof(double) * ((size << 1) +1));
	double* trans_two = (double*) malloc(sizeof(double) * ((size << 1) +1));

	double* input_filled = (double*) malloc(sizeof(double) * (3*size +1));
	double* input_2_filled = (double*) malloc(sizeof(double) * (3*size +1));

	zeroadding(size, input, input_filled);
	zeroadding(size, input_2, input_2_filled);
	size = 2*size +1;

	int i = 0;
	int x = 0;



	double* conv = (double*) malloc(sizeof(double) * (size << 1) +2);


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
	}

	for(i=0; i < size; i++)
	{
		for(x=0; x < size; x++)
		{

			output[i << 1] += (conv[x << 1]*coef_cos[i*size +x] - conv[(x << 1) +1]*coef_sin[i*size +x])/(double)size;
			//		output[(i << 1) +1] += (conv[x << 1]*(coef_sin[i*size +x]) + conv[(x << 1) + 1] * coef_cos[i*size +x])/(double)size;


		}
	}



	free(trans_one);
	free(trans_two);
	free(input_filled);
	free(input_2_filled);


}


void fastdft(int size, double input, double* coef_cos, double* coef_sin, double* Re, double* Im)
{
        int i = 0;
        int j = 0;
        int x = 0;
        int g = 0;

        double* input_filled = (double*) malloc(sizeof(double) * (2*size +2));
        zeroadding(size, input, input_filled);
        size = 2*size +1;
        double* Re_matrix = (double*) malloc(sizeof(double)*(size*size + size + 2)); 
        double* Im_matrox = (double*) malloc(sizeof(double)*(size*size + size + 2));
        for(j=0; j < size; j++)
        {
                for(i=j; i < size; i++)
                {       
                        Re_matrix[size*j + i] = coef_cos[i*size + j] * input_filled[i];
                        Im_matrix[size*j + i] = coef_sin[i*size + j] * input_filled[i];
                
                        Re_matrix[size*i +j] = Re_matrix[size*j +i];
                        Im_matrix[size*i +j] = Im_matrix[size*j +i];
                        
                }
        
        }

        for(j = 0; j  < size; j++)
        {
                for(i=0; i < size; i++)
                {
                        Re[j] += Re_matrix[size*i +j];
                        Im[j] += Im_matrix[size*i +j];

                }

        }

}











#endif

