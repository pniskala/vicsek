#ifndef SUPER_H
#define SUPER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include "randomlib.h"


typedef struct structSimulation Simulation;

typedef struct Cell_s Cell;

typedef struct Bird_s Bird;


Simulation* init_simulation(int N, double L, double v0, double tstep, int tmax, int tmin, double range, double noise);

Bird** get_birds(Simulation* simu);

Cell*** get_cells(Simulation* simu);

void move_birds(Simulation* simu);

double update_birds(Simulation* simu);

double calc_vave(Simulation* simu);

void calc_angles(Simulation* simu);

double run(Simulation* run, int write);

double run_one_step(Simulation* simu);

double update_mva_list(Simulation* simulation, double vave);

double distance(Bird* bird1, Bird* bird2, double L);

double get_coord(double coord, double L);

int get_cc(int cc, double L);

double get_noise(double eka);

long int getseed();

void print_birds(Simulation* simu, FILE* file);

void free_birds(Simulation* simu);

void free_simu(Simulation* simu);

void free_cells(Simulation* simu);

Cell* create_cell(int cx, int cy);

int get_cell_cx(Cell* cell);

void add_bird(Cell* cell, Bird* bird);

void remove_bird(Cell* cell, Bird* bird);

Bird* create_bird(Simulation* simu, double L, double v0, int idx);

void move(Bird* bird, double L, double tstep);

void update(Bird* bird);

void set_new_angle(Bird* bird, double na);

void assign_cc(Bird* bird);

void assign_ncells(Bird* bird, double L);

Cell** get_ncells(Bird* bird);

#endif
