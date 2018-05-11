#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <assert.h>
#include <time.h>
#include <string.h>

#include "super.h"
#include "randomlib.h"

//runs the simulation with several noise values with the given parameters, calculating the average velocity with the respective error
void noise_simulation(int N, double L, double v0, double tstep, int tmax, int tmin, double range, double etamax, double etastep, int nsamples)
{
	//printf("kurr");
	double vave;
	double va_cum;
	double va_2_cum;
	double va_4_cum;
	double binder;
	double error;
	
	//write the input into the same directory as the results, just to be sure
	
	//printf("kurr\n");
	char filunimi3[100]; 
	sprintf(filunimi3,"/scratch/networks/pniskala/noise_simulations/N = %d L = %d/input.txt", N, (int)L);
	FILE* input = fopen(filunimi3,"w");
	fprintf(input,"N = %d\n",N);
	fprintf(input,"L = %.1f\n",L);
	fprintf(input,"v0 = %.2f\n",v0);
	fprintf(input,"range = %.1f\n",range);
	fprintf(input,"tstep = %.1f\n",tstep);
	fprintf(input,"etamax = %.1f\n",etamax);
	fprintf(input,"etastep = %.1f\n",etastep);
	fprintf(input,"tmax = %d\n",tmax);
	fprintf(input,"nsamples = %d\n",nsamples);
	fclose(input);
	
	//define the filename and path for the result and raw data files.
	char filename1[100];
	char filename2[100];
	
	sprintf(filename1,"/scratch/networks/pniskala/noise_simulations/N = %d L = %d/results.txt", N,(int)L);
	sprintf(filename2,"/scratch/networks/pniskala/noise_simulations/N = %d L = %d/raw_data.txt", N,(int)L);
	FILE* filu;//= fopen(filename1,"w");
	FILE* raw_data;// = fopen(filename2,"w");
	printf("Noise simulation is starting...\n");

	//this loop handles the different eta values
	for (double eta = 0.0; eta <= etamax; eta += etastep)
	{
		va_cum = 0;
		va_2_cum = 0;
		va_4_cum = 0;
		error = 0;
		binder = 0;
		
		//this loop runs nsamples number of simulations and saves the data
		for (int i = 0; i < nsamples; i++)
		{
			Simulation* simu = init_simulation(N,L,v0,tstep,tmax,tmin,range,eta);
			vave = run(simu,0);
			raw_data = fopen(filename2,"a"); //open the raw data file
			fprintf(raw_data,"%f %f\n", eta, vave); //write the sample data
			fclose(raw_data); //close the raw data file
			va_cum += vave;
			va_2_cum += vave*vave;
			va_4_cum += vave*vave*vave*vave;
			free_simu(simu); //free the allocated memory
		}
		//calculate the averages
		double va = va_cum / (double)nsamples;
		double va_2 = va_2_cum / (double)nsamples;
		double va_4 = va_4_cum / (double)nsamples;
		
		//calculate the moments / cumulants
		error = va_2-va*va;
		binder = 1.0 - va_4/(3.0*va_2*va_2);
		
		//write the data to files
		filu = fopen(filename1,"a");
		printf("%f %f %f %f\n",eta,va,error,binder);
		fprintf(filu,"%f %f %f %f\n",eta,va,error,binder);
		fclose(filu);
	}
	

}

//runs the simulation with several noise values with the given parameters, calculating the average velocity with the respective error
void density_simulation(double eta, double L, double v0, double tstep, int tmax, int tmin, double range, double Nmax, int nsamples)
{
	//Simulation* simu = NULL;
	double vave;
	double va_cum;
	double va_2_cum;
	double va_4_cum;
	double binder;
	double error;
	
	char filename1[100];
	char filename2[100];
	
	char filunimi3[100]; 
	
	sprintf(filunimi3,"/scratch/networks/pniskala/eta = %.1f L = %d/input.txt", eta, (int)L);
	FILE* input = fopen(filunimi3,"w");
	fprintf(input,"eta = %.1f\n",eta);
	fprintf(input,"L = %.1f\n",L);
	fprintf(input,"v0 = %.2f\n",v0);
	fprintf(input,"range = %.1f\n",range);
	fprintf(input,"tstep = %.1f\n",tstep);
	fprintf(input,"Nmax = %d\n",(int)Nmax);
	//fprintf(input,"etastep = %.1f\n",etastep);
	fprintf(input,"tmax = %d\n",tmax);
	fprintf(input,"nsamples = %d\n",nsamples);
	fclose(input);
	
	//create the right paths for the files
	sprintf(filename1,"/scratch/networks/pniskala/eta = %.1f L = %d/results.txt", eta,(int)L);
	sprintf(filename2,"/scratch/networks/pniskala/eta = %.1f L = %d/raw_data.txt", eta,(int)L);
	
	//
	FILE* filu;// = fopen(filename1,"w");
	FILE* raw_data;// = fopen(filename2,"w");
	printf("Density simulation is starting...\n");
	//printf("%f",etamax);
	for (double N = 80; N <= Nmax; N += 80)
	{
		va_cum = 0;
		va_2_cum = 0;
		va_4_cum = 0;
		error = 0;
		binder = 0;
		for (int i = 0; i < nsamples; i++)
		{
			Simulation* simu = init_simulation(N,L,v0,tstep,tmax,tmin,range,eta);
			vave = run(simu,0);
			raw_data = fopen(filename2,"a");
			fprintf(raw_data,"%f %f\n", ((double)N)/(L*L), vave);
			fclose(raw_data);
			va_cum += vave;
			va_2_cum += vave*vave;
			va_4_cum += vave*vave*vave*vave;
			//vaves+=vavelist[i];
			free_simu(simu);
		}
		//calculate the averages
		double va = va_cum / (double)nsamples;
		double va_2 = va_2_cum / (double)nsamples;
		double va_4 = va_4_cum / (double)nsamples;
		
		//calculate the moments / cumulants
		error = va_2-va*va;
		binder = 1.0 - va_4/(3.0*va_2*va_2);
		
		//write the data to files
		filu = fopen(filename1,"a");
		printf("%f %f %f %f\n",((double)N)/(L*L),va,error,binder);
		fprintf(filu,"%f %f %f %f\n",((double)N)/(L*L),va,error,binder);
		fclose(filu);
	}
}


//for running several simulations with same parameters and investigating the time evolution
void time_average(int N, double L, double v0, double tstep, double eta, int tmax, int tmin, double range, int nsamples)
{

	printf("Starting time evolution thingy\n");
	FILE* data;
	//FILE* raw_data = fopen("/scratch/networks/pniskala/newest_run/raw_data.txt","w");
	
	Simulation** simu_list = malloc(nsamples*sizeof(Simulation*));
	double* vave_list = malloc(tmax*sizeof(double));
	
	char filename1[110];
	
	sprintf(filename1,"./time_ave_N%d_L%d_eta%.1f_nsamples%d.txt",N,(int)L,eta,nsamples);
	//printf("%s\n",filename1);
	data = fopen(filename1,"w");
	fclose(data);
	
	for (int i = 0; i < nsamples; i++)
	{
		simu_list[i] = init_simulation(N, L, v0, tstep, tmax, tmin, range, eta);
	}
	
	for (int t = 0; t < tmax; t++)
	{
		vave_list[t] = 0;	
	}
	
	for (int t = 0; t < tmax; t++)
	{
		double average = 0;
		for (int i = 0; i < nsamples; i++)
		{
			average += run_one_step(simu_list[i]);
		}
		vave_list[t] = average/nsamples;
		data = fopen(filename1,"a");
		fprintf(data, "%d %f\n", t, (average/((double)nsamples)));
		printf("%d %f\n", t, (average/((double)nsamples)));
		fclose(data);
	}
	
	for (int i = 0; i < nsamples; i++)
	{
		free_simu(simu_list[i]);
	}
	free(simu_list);
	free(vave_list);
}


void read_input()
{
	printf("reading input...\n");
	FILE* file = fopen("input.txt", "r");
	if (file != NULL)
	{
		char line[128];
		char* value;
		int mode, N, tmax, tmin, nsamples;
		double L, v0, range, tstep, eta, etamax, Nmax, etastep;
		while (fgets(line, 128, file) != NULL)
		{
			char* par = strtok(line," =");
			if (strcmp(par,"mode")==0)
			{
				value = strtok(NULL," =");
				printf("mode = %s\n",value);
				mode = atoi(value);
			} else if (strcmp(par,"N")==0)
			{
				value = strtok(NULL," =");
				printf("N = %s\n",value);
				N = atoi(value);
			} else if (strcmp(par,"L")==0)
			{
				value = strtok(NULL," =");
				printf("L = %s\n",value);
				L = atof(value);
			} else if (strcmp(par,"v0")==0)
			{
				value = strtok(NULL," =");
				printf("v0 = %s\n",value);
				v0 = atof(value); 
			} else if (strcmp(par,"tstep")==0)
			{
				value = strtok(NULL," =");
				printf("tstep = %s\n",value);
				tstep = atof(value);
			} else if (strcmp(par,"nsamples")==0)
			{
				value = strtok(NULL," =");
				printf("nsamples = %s\n",value);
				nsamples = atoi(value);
			} else if (strcmp(par,"tmax")==0)
			{
				value = strtok(NULL," =");
				printf("tmax = %s\n",value);
				tmax = atoi(value);
			} else if (strcmp(par,"range")==0)
			{
				value = strtok(NULL," =");
				printf("range = %s\n",value);
				range = atof(value);
			} else if (strcmp(par,"eta")==0)
			{
				value = strtok(NULL," =");
				printf("eta = %s\n",value);
				eta = atof(value);
			} else if (strcmp(par,"etamax")==0)
			{
				value = strtok(NULL," =");
				printf("etamax = %s\n",value);
				etamax = atof(value);
			} else if (strcmp(par,"Nmax")==0)
			{
				value = strtok(NULL," =");
				printf("Nmax = %s\n",value);
				Nmax = atof(value);
			} else if (strcmp(par,"etastep")==0)
			{
				value = strtok(NULL," =");
				printf("etastep = %s\n",value);
				etastep = atof(value);
			} else if (strcmp(par,"tmin")==0)
			{
				value = strtok(NULL," =");
				printf("tmin = %s\n",value);
				tmin = atoi(value);
			}
		}
		fclose(file);
		if (mode == 1) //runs a single simulation run and writes the bird positions to a file
		{
			printf("single simulation chosen\n");
			Simulation* simu = init_simulation(N,L,v0,tstep,tmax,tmin,range,eta);
			run(simu,1);
			free_simu(simu);
		} else if (mode == 2) //runs the simulation nsamples times per noise value upto etamax (max noise)
		{
			printf("noise simulation chosen\n");
			noise_simulation(N,L,v0,tstep,tmax,tmin,range,etamax,etastep,nsamples);
		} else if (mode == 3) //runs the simulation nsamples times per density value upto Nmax (max particle number)
		{
			printf("density simulation chosen\n");
			density_simulation(eta,L,v0,tstep,tmax,tmin,range,Nmax,nsamples);
		}
		else if (mode == 4) //runs several simulations with the same parameters
		{
			printf("time averaging simulation chosen\n");
			time_average(N,L,v0,tstep,eta,tmax,tmin,range,nsamples);
		} else //mode selection invalid or no mode selected in the input file
		{
			printf("No mode selected!");
		}
		
	} else
	{
		printf("File didn't open! Check if the file exists or call Batman.\n");
	}
}

void test_random()
{
	FILE* rfilu = fopen("random_test.txt","w");
	for (int i = 0; i < 726800; i++)
	{
		double randomi = get_noise(3.0);
		fprintf(rfilu, "%f\n",randomi);
	}
	for (int i = 0; i < 20000; i++)
	{
		double randomi = get_noise(3.0);
		fprintf(rfilu, "%f\n",randomi);
	}
	for (int j = 0; j < 6; j++)
	{
		int random2 = RandomUniform()*(4096*4096);
		printf("%d\n",random2);
	}
	fclose(rfilu);
}

int main(void)
{
	int ti = time(NULL);

	unsigned long int sec = time(NULL);
	unsigned int IJ = sec % 30000;
	unsigned int KL = sec % 4241;
	
	//printf("%d %d\n",IJ,KL);
	//IJ = 1802;
	//KL = 9373;
	
	RandomInitialise(IJ,KL);
	
	read_input();
	//test_random();

	unsigned long int te = time(NULL);
	printf("It took %ld seconds to run the simulation. Also, I'm batman.\n",te-ti);	
}


