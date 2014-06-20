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
                        printf("\n \n \n \t \t %f", turn_perc);
                        printf("\n %f", turn);

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
                printf("\n ----------------------\n ");
        for(x=0; x < size; x++)
                {
                        output[i << 1] += input[x]*coef_cos[i*size +x];
                        output[(i << 1) +1] += input[x]*(coef_sin[i*size +x]*(-1));


                        printf("\t %f \t \t %f \t %f \t %f \t %f \t", input[x], output[i << 1], output[(i << 1) +1], cos(2*M_PI*i*x/size), sin(2*M_PI*i*x/size));

                        printf("%f \t %f \t \n", coef_cos[i*size +x], coef_sin[i*size +x]);



                }
        }
        return 1;
}





#endif
