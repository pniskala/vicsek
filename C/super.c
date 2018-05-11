#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <assert.h>

#include "randomlib.h"
#include "super.h"

//define the global constant of pi
const double PI = 3.14159265358979323846;
const double TWOPI = 2.0*3.14159265358979323846;


//define the Simulation struct (and type), containing all the information about the simulation
struct structSimulation
{
	Bird** birds;
	Cell*** cells;
	int N;
	double L;
	double v0;
	double noise;
	double tstep;
	double range;
	int time;
	int tmax;
	int tmin;
	double mva_list[200]; //moving average array
	Bird* bfirst;
};

//The simulation area is divided in to cells. A cell knows the list of birds inside the cell.
struct Cell_s
{
	int nbirds; //how many birds are in the cell
	//Bird** birds;
	int cx;
	int cy;
	Bird* bfirst;
	Bird* blast;
	//int* bidx; //bird index list
};


//Birds live in the cells. They know the cell they are in, and the neighboring cells.
struct Bird_s
{
	double a; //angle
	double na; //new angle
	double x;
	double y;
	int cx; //cell coordinate
	int cy; //cell coordinate
	double vx;
	double vy;
	double v0;
	Cell* ncells[9];
	Bird* snext;

	Bird* cnext;
	Bird* cprev;
	double L;
	//Simulation* simu;
	int idx;
	Simulation* simu;
	double cum_vx; //for calculating the other average angle
	double cum_vy; //for calculating the other average angle
};


//Initializes the simulation
Simulation* init_simulation(int N, double L, double v0, double tstep, int tmax, int tmin, double range, double noise)
{
	Simulation* simu = malloc(sizeof(Simulation));
	simu->N = N;
	simu->L = L;
	simu->v0 = v0;
	simu->tstep = tstep;
	simu->range = range;
	simu->time = 0;
	simu->tmax = tmax;
	simu->tmin = tmin;
	simu->noise = noise;
	simu->birds = malloc(N*sizeof(Bird*));
	simu->cells = malloc(((int)L)*sizeof(Cell**));
	
	//create the cells
	for (int i = 0; i < (int)simu->L; i++)
	{
		simu->cells[i] = malloc(((int)L)*sizeof(Cell*));
		for (int j = 0; j < (int)simu->L; j++)
		{
			simu->cells[i][j] = create_cell(i,j);
		}
	}
	
	simu->birds[0] = create_bird(simu,L,v0,0);
	simu->bfirst = simu->birds[0];
	assign_cc(simu->birds[0]);
	assign_ncells(simu->birds[0],L);
	add_bird(simu->birds[0]->ncells[0], simu->birds[0]);
	//create the birds and add them to the cells
	for (int i = 1; i < simu->N; i++)
	{
		simu->birds[i] = create_bird(simu,L,v0,i);
		simu->birds[i-1]->snext = simu->birds[i];
		assign_cc(simu->birds[i]);
		assign_ncells(simu->birds[i],L);
		add_bird(simu->birds[i]->ncells[0], simu->birds[i]);
	}
	
	for (int i = 0; i < 200; i++)
	{
		simu->mva_list[i] = 0;
	} 
	
	
	return simu;
}

Bird** get_birds(Simulation* simu)
{
	assert(simu->birds);
	return simu->birds;
}

Cell*** get_cells(Simulation* simu)
{
	assert(simu);
	return simu->cells;
}


//moves all the birds
void move_birds(Simulation* simu)
{
	assert(simu);
	for (int i = 0; i < simu->N; i++)
	{
		move(simu->birds[i],simu->L, simu->tstep);
	}
}


//updates all the birds angles and moves them
double update_birds(Simulation* simu)
{
	double vx = 0.0;
	double vy = 0.0;
	Bird* next = simu->bfirst;
	while (next != NULL)
	{
		move(next,simu->L, simu->tstep);
		update(next);
		vx += next->vx;
		vy += next->vy;
		next = next->snext;
	}
	return sqrt(vx*vx + vy*vy)/(simu->v0*(double)simu->N);
}

//calculates the (absolute) average velocity for dem birds
double calc_vave(Simulation* simu)
{
	double vx = 0;
	double vy = 0;
	Bird* next = simu->bfirst;
	while (next != NULL)
	{
		vx += next->vx;
		vy += next->vy;
		next = next->snext;
	}
	return sqrt(vx*vx + vy*vy)/(simu->v0*(double)simu->N);
}


//this function calculates the new angles for the birds
void calc_angles(Simulation* simu)
{
	//assert(simu);
	Bird* snext = simu->bfirst;
	while (snext != NULL)
	{
		double vx = 0;
		double vy = 0;
		for (int j = 0; j < 9; j++)
		{
			Bird* next = snext->ncells[j]->bfirst;
			while (next != NULL)
			{
				if (distance(snext, next, simu->L) < simu->range)
				{
					vx += next->vx;
					vy += next->vy;
				}
				next = next->cnext;
			}
		}
		snext->na = atan2f(vy,vx) + get_noise(simu->noise);
		snext = snext->snext;
	}
}


//runs the simulation
double run(Simulation* simu, int write)
{
	assert(simu);
	FILE* loc_file;
	FILE* time_vave_file;
	FILE* time_cumave;
	double cnew = 0, cold = 1;
	double mva = 0;
	if (write==1)
	{
		loc_file = fopen("./birdloc.txt","w");
		time_vave_file = fopen("./time_vave.txt","w");
		time_cumave = fopen("./time_cumave.txt","w");
	}
	double vave = calc_vave(simu);
	cnew = vave;
	for (int i = 0; i < simu->tmax; i++)
	{
		cold = cnew;
		if (write==1)
		{	
			print_birds(simu,loc_file);
			fprintf(time_vave_file,"%d %f %f\n",simu->time,vave,cnew);
			printf("%d %f %f %f\n",simu->time,vave,cnew,mva);
		}
		calc_angles(simu);
		vave = update_birds(simu);
		simu->time++;
		cnew = (cold*((double)simu->time-1.0)+vave)/((double)(simu->time));
		mva = update_mva_list(simu,cnew);
		if ((fabs(mva-simu->mva_list[160]) < 0.00001) && (simu->time > simu->tmin))
		{
			break;
		}
		if (vave > 0.9999)
		{
			vave = 1.0;
			break;
		}
	}
	if(write==1)
	{
		fclose(time_vave_file);
		fclose(loc_file);
		fclose(time_cumave);
	}
	printf("%d %f\n",simu->time,vave);
	return vave;
}


//runs the simulation for one time step
double run_one_step(Simulation* simu)
{
	calc_angles(simu);
	double vave = update_birds(simu);
	simu->time++;
	return vave;
}

//updates the moving average array by putting the given value to the last place
double update_mva_list(Simulation* simu, double vave)
{
	double mva = vave;
	for (int i = 1; i < 200; i++)
	{
		simu->mva_list[i-1] = simu->mva_list[i];
		mva += simu->mva_list[i-1];
	}
	simu->mva_list[199] = vave;
	return (mva /= 200.0);
}

//calculates the distance between two birds and applies the boundary conditions in the process
double distance(Bird* bird1, Bird* bird2, double L)
{
	double Lp2 = L/2.0;
	double dx = fabs(bird1->x-bird2->x);
	double dy = fabs(bird1->y-bird2->y);
	if (dx > Lp2)
	{
		dx = L - dx;	
	} else if (dx < -Lp2)
	{
		dx = L + dx;
	}
	if (dy > Lp2)
	{
		dy = L - dy;
	}
	return dx*dx + dy*dy;
}

//applies the boundary conditions
double get_coord(double coord, double L)
{
	if (coord > L)
	{
		coord = coord - L;
	} else if (coord < 0)
	{
		coord = coord + L;
	}
	return coord;
}

int get_cc(int cc, double L)
{
	if (cc < 0)
	{
		cc = cc + (int)L;
	} else if (cc > (int)L - 1)
	{
		cc = cc - (int)L;
	}	
	return cc;
}

double get_noise(double eta)
{
	return eta*(RandomUniform()-0.5);
}

void free_birds(Simulation* simu)
{
	//assert(simu);
	for (int i = 0; i < simu->N; i++)
	{
		free(simu->birds[i]);
	}
	free(simu->birds);
}

void free_cells(Simulation* simu)
{
	assert(simu);
	for (int i = 0; i < (int)simu->L; i++)
	{
		for (int j = 0; j < (int)simu->L; j++)
		{
			//free(simu->cells[i][j]->birds);
			free(simu->cells[i][j]);
		}
		free(simu->cells[i]);
	}
	free(simu->cells);
}

void free_simu(Simulation* simu)
{
	free_birds(simu);
	free_cells(simu);	
	free(simu);
}

void print_birds(Simulation* simu, FILE* file)
{
	assert(simu);
	for (int i = 0; i < simu->N; i++)
	{
		fprintf(file,"%d %f %f %f\n",simu->time, simu->birds[i]->x,simu->birds[i]->y,simu->birds[i]->a);
	}
}

//CELLS

Cell* create_cell(int cx, int cy)
{
	Cell* cell = malloc(sizeof(Cell));
	cell->cx = cx;
	cell->cy = cy;
	cell->nbirds = 0;
	cell->bfirst = NULL;
	cell->blast = NULL;
	return cell;
}

int get_cell_cx(Cell* cell)
{
	return cell->cx;
}

void add_bird(Cell* cell, Bird* bird)
{
	if (cell->nbirds == 0)
	{
		cell->bfirst = bird;
		cell->blast = bird;
	} else
	{
		bird->cprev = cell->blast;
		cell->blast->cnext = bird;
		cell->blast = bird;
	}
	cell->nbirds++;
}

void remove_bird(Cell* cell, Bird* bird)
{
	if ((bird->cprev == NULL) && (bird->cnext == NULL))
	{
		cell->bfirst = NULL;
		cell->blast = NULL;
	} else if ((bird->cprev != NULL) && (bird->cnext == NULL))
	{
		cell->blast = bird->cprev;
		bird->cprev->cnext = NULL;
		bird->cprev = NULL;
	} else if ((bird->cprev == NULL) && (bird->cnext != NULL))
	{
		cell->bfirst = bird->cnext;
		bird->cnext->cprev = NULL;
		bird->cnext = NULL;
	} else
	{
		bird->cprev->cnext = bird->cnext;
		bird->cnext->cprev = bird->cprev;
		bird->cprev = NULL;
		bird->cnext = NULL;
	}
	cell->nbirds--;
}

//BIRDS

Bird* create_bird(Simulation* simu, double L, double v0, int idx)
{
	Bird* newBird = malloc(sizeof(Bird));
	newBird->simu = simu;
	newBird->v0 = v0;
	newBird->L = L;
	newBird->x = (double)L*RandomUniform();
	newBird->y = (double)L*RandomUniform();
	newBird->cx = (int)newBird->x;
	newBird->cy = (int)newBird->y;
	newBird->a = (double)TWOPI*RandomUniform();
	newBird->na = newBird->a;
	newBird->vx = v0*cosf(newBird->a);
	newBird->vy = v0*sinf(newBird->a);
	newBird->idx = idx;
	newBird->cprev = NULL;
	newBird->cnext = NULL;
	newBird->snext = NULL;
	newBird->cum_vx = 0;
	newBird->cum_vy = 0;
	return newBird;
}

//Moves a bird and, if necessary, updates the list of neighboring cells.
void move(Bird* bird, double L, double tstep)
{
	bird->x = get_coord(bird->x + bird->vx*tstep, L);
	bird->y = get_coord(bird->y + bird->vy*tstep, L);
	int ocx = bird->cx;
	int ocy = bird->cy;
	assign_cc(bird);
	if ((ocx != bird->cx) || (ocy != bird->cy))
	{
		remove_bird(bird->ncells[0], bird);
		assign_ncells(bird,L);
		add_bird(bird->ncells[0], bird);
	}
}

void update(Bird* bird)
{
	bird->a = bird->na;
	if (bird->na > TWOPI)
	{
		bird->a = bird->a - TWOPI;
	} else if (bird->na < 0)
	{
		bird->a = bird->a + TWOPI;
	}
	bird->vx = bird->v0*cosf(bird->a);
	bird->vy = bird->v0*sinf(bird->a);
}

//sets the new angle for the bird and limits it to interval from 0 to 2pi
void set_new_angle(Bird* bird, double na)
{
	bird->na = na;	
}

//gives the bird its cell coordinates according to its location
void assign_cc(Bird* bird)
{
	bird->cx = (int)bird->x;
	bird->cy = (int)bird->y;	
}

void assign_ncells(Bird* bird, double L)
{
	int cx = bird->cx;
	int cy = bird->cy;
	int xright = get_cc(cx+1, L);
	int xleft = get_cc(cx-1, L);
	int yup = get_cc(cy+1, L);
	int ydown = get_cc(cy-1,L);
	
	bird->ncells[0] = bird->simu->cells[cx][cy];
	bird->ncells[1] = bird->simu->cells[xleft][cy];
	bird->ncells[2] = bird->simu->cells[xleft][yup];
	bird->ncells[3] = bird->simu->cells[cx][yup];
	bird->ncells[4] = bird->simu->cells[xright][yup];
	bird->ncells[5] = bird->simu->cells[xright][cy];
	bird->ncells[6] = bird->simu->cells[xright][ydown];
	bird->ncells[7] = bird->simu->cells[cx][ydown];
	bird->ncells[8] = bird->simu->cells[xleft][ydown];
}

Cell** get_ncells(Bird* bird)
{
	return bird->ncells;
}

