/* *************************
 * Author: xbx1992
 * Created Time:  Tue 11  July 2015
 * LastModified:  Tue 11July 2015 04:00:27 PM CST
 * C File Name:  configclass.h
 * ************************/
#ifndef CONFIGCLASS_H
#  define CONFIGCLASS_H

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <vector>

#include "xg_datetime.h"
using namespace std;

struct daconfig{
	/* ensemble size*/
	size_t K;

	/*grid information*/
	double grid_xres;
	double grid_yres;
	int XSIZE;
	int YSIZE;
	
	/*start time,lag time and  end time*/
	datetime start_date;	
	datetime end_date;

	/*time step*/
	int DA_length;
	
	int DA_itr;
	
	
	/*number of lag*/
	int nlag;

	/*inflation factor*/
	double rho;

	/*number of region's classification*/
	size_t Class;
	
	/*size of state vector*/
	size_t M;
};

struct fileconfig{
	/*DA diractory*/
	const char * DA_dir;
	
	/*input diractory*/
	const char * in_dir;

	/*out diractory*/
	const char * out_dir;

	const char * lai_dir;	
	
	/*landcover file name*/
	const char * landcover_file;

	/*BPLUT file name*/
	const char * BPLUT_file;

	/*observation file name*/
	const char * obs_f;
};

struct observation{
    double lat, lon;
    int geosI,geosJ,alt;
    double tau;
    double value,mdm;
    vector<double> modeled_v;
};


struct mod_data{
    double geosI,geosJ;
    double tau;
    double value,mdm;
};

struct mpi_config{
    int size, rank, name_len;
    char* host_name;
};

typedef vector<observation> OBS_DATA;

typedef map<int, map<int, int > > REGIONS_MAP;

typedef map<pair<int, int>, vector<pair<double, int> > > OBS_INDEX;
#endif
