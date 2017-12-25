/* *************************
 * Author: xbx1992
 * Created Time:  
 * LastModified:  Tue 11 Dec 2014 04:00:27 PM CST
 * C File Name: 
 * ************************/
#include <iostream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <set>
#include "bpch.h"
#include "common.h"
#include "xg_datetime.h"
#include "state.h"
#include "init.h"
#include "configclass.h"
#include "observationdata.h"
#include "observeoperation.h"
#include "sample.h"
#include "ensrf.h"
#include <Eigen/Dense>
#include "mpi.h"
#include "xg_math_vector.h"
#include "xg_code.h"
//#ifdef _OPENMP
//#include "omp.h"
//#endif
using namespace Eigen;
using namespace std;





 typedef map<int, map<int, int > > REGIONS_MAP;
REGIONS_MAP regions;
map<int, vector<pair<int, int> > > regions_index;

void read_Regions(const char *land_f)
{
    	printf("read lc:%s\n",land_f);
	    GDALDataset *data = (GDALDataset*) GDALOpen(land_f, GA_ReadOnly);
    	assert(data != NULL);
    	GDALRasterBand *band = data -> GetRasterBand(1);
    	assert(band != NULL);

	int NaNs_land = band -> GetNoDataValue();
    	int YSize = band->GetYSize(), XSize = band->GetXSize();
		for(int line = 0; line < YSize; ++line){
        	vector<double> lines;
            	lines = ReadLine(band, line);
        	for(int i = 0; i < XSize; ++i){
			regions[line+1][i+1] = lines[i];
			regions_index[lines[i]].push_back(make_pair<int, int>(line+1,i+1));	
		}
	}
    	GDALClose(data);
	
}



int main(int argc, char *argv[])
{
    /**********************************
     * Init MPI environment
     * *******************************/
    int mpi_size, mpi_rank, mpi_name_len;
    char mpi_host_name[MPI_MAX_PROCESSOR_NAME];
    mpi_config mpic;
    MPI_Init(&argc, &argv); // initialize MPI
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Get_processor_name(mpi_host_name, &mpi_name_len);
    mpic = SAVE_MPI_Config(mpi_size, mpi_rank, mpi_name_len, mpi_host_name);
    /************************************
     * read configuration file 
     ************************************/

    DAConfig da_config;
    daconfig dc;
    fileconfig fc;
    da_config.Initialize(argv[1],dc,fc);
    
    /**********************************
     * Register GDAL
     **********************************/
    GDALAllRegister();
    /* load  REGIONS*/
    

    read_Regions(fc.landcover_file); //read landcover file
	/*****************************************************
     * Set up the initial state vector x (SCALING FACTOR)
     *****************************************************/
    /*
     * generate the first ensemble of state vector x 
     *
     *  DONE by PROCESS-0, then broadcast to MPI_COMM_WORLD
     */
    statevector<double> x_b_lagMembers(dc.nlag , dc.Class, dc.K);
    OBS_DATA obs;
    /************************************************
     *
     * START first year
     *
     ************************************************/
	x_b_lagMembers = generate_x(mpi_rank,dc,regions_index);
     debug("finish generate_x,start assimilation cycle");
    /************************************************
     *
     * START ASSIMILATION CYCLE
     *
     ************************************************/
    clock_t time_first = clock(),time_second;
     
    
    for(datetime DA_i = dc.start_date; DA_i < dc.end_date; DA_i += timespan(0,0,0,0,dc.DA_length)){

        /* set up time */
 	datetime DA_start = DA_i;
        datetime DA_end   = DA_i + timespan(0,0,0,0, dc.DA_length);
        double   tau_start= DA_start.tau();
        double   tau_end  = DA_end  .tau();
        printf("DA_length = %d, START ASSIMILATION CYCLE %s(%lf) - %s(%lf)\n", dc.DA_length, DA_start.str().c_str(), tau_start, DA_end.str().c_str(), tau_end);
        /************************
         * read observation data
         *
         * which is y_o
         ************************/
	 obs_c DA_Obs;
 	 obs = DA_Obs.Initialize(DA_start, DA_end, dc, fc);
	 map<int, double> R = DA_Obs.Get_R(obs); 
        /* build up index to speed up searching */
	OBS_INDEX obs_index = build_index(obs);
	/****************************
	 * statevector to grid
	 *
	 * *****************************/
	 char mk_str[128];
            
    	if(mpi_rank == 0) {
       	    	sprintf(mk_str, "mkdir %s/diagnose", fc.DA_dir);
            	cmd(mk_str);
            	
		sprintf(mk_str, "statevector.%s.nc", DA_i.str("YYYYMMDDhhmmss").c_str());
    	    	bool flag = x_b_lagMembers.save_apri_vector(mk_str);	
	    }   
	 
	    /***************************
	     * run ensemble
	     *
	     *************************/
 	
	    obs_operation OBS_Operation;	
	
	    OBS_Operation.Initialize(dc, fc,mpic, mpi_rank, x_b_lagMembers, obs, DA_start, DA_end);
		
   	 
	    MPI_Barrier(MPI_COMM_WORLD);
	 
    
	    double *x_b_p = (double*)malloc( dc.nlag * dc.Class * dc.K * sizeof(double));
	
	if(mpi_rank == 0){
          
		/*Sample*/
		
		printf("begin sample!"); 
    		time_second = Cal_time(time_first); 
		
		Sample Sample_obs;
		
		Sample_obs.Initialize(obs, DA_start, DA_end,dc, fc);
		
		printf("finsih sample!"); 
    		time_first = Cal_time(time_second); 
            /*
             *start assimilate GRID POINT by GRID POINT 
             */
            	debug("START ASSIMILATING REGION");

		optimize OP;
		OP.Initialize(x_b_lagMembers, obs, dc);
		MatrixXf x_b_bar(dc.M, 1); 
		x_b_lagMembers = OP.Run(dc, obs, R, x_b_bar, regions, regions_index,tau_start,fc,DA_start,DA_end);	
		     /*****************
 		     *
                     * STEP 5 generate all_scaling
                     * ***************/
		    debug("generate all scaling");
		    OBS_Operation.Run_posterior(dc, fc, x_b_bar,obs,DA_start,DA_end);
		    
		    /******************
 		     *STEP 6 Propagate
 		     *****************/
                    debug("save post");
		    sprintf(mk_str, "statevector.%s.nc", DA_i.str("YYYYMMDDhhmmss").c_str());
    	    	    bool flag = x_b_lagMembers.save_post_vector(mk_str);	
            	    sprintf(mk_str, "mv statevector.%s.nc %s/diagnose/",DA_i.str("YYYYMMDDhhmmss").c_str(), fc.DA_dir);
            	    cmd(mk_str);
		    
		    debug("propagate");

		    bool rt = x_b_lagMembers.Propagate(dc,x_b_bar);
    		   		    
		    x_b_p = x_b_lagMembers.expand();
		    //debug(x_b_p);
		    MPI_Bcast(x_b_p, dc.M * dc.K, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	            /****************
 		    *
 		    * print P_a
 		    * ***************/ 
		    debug("P_a");
		      
		}	
   		else{
   		//debug(x_b_p); 
		MPI_Bcast(x_b_p, dc.M * dc.K, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		reverse_expand(x_b_lagMembers,x_b_p, dc.nlag, dc.Class,dc.K);
	 
		}
   	 	MPI_Barrier(MPI_COMM_WORLD);
		free(x_b_p);
	    }// END OF DA CYCLE
	 	
	    MPI_Finalize();
	    return 0;
	}

			    


