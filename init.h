/* *************************
 * Author: xbx1992
 * Created Time:  Tue 11  July 2015
 * LastModified:  Tue 11July 2015 04:00:27 PM CST
 * C File Name:  init.h
 * ************************/
#ifndef INIT_H
#  define INIT_H

#include <cstring>
#include <map>
#include "common.h"
#include "xg_datetime.h"
#include "state.h"
#include "configclass.h"

class DAConfig{
	private:
	map<string, string> DA_config;

	public:
	DAConfig(){}
	void Initialize(const char * conf, daconfig& dc, fileconfig& fc)
	{
    		/************************************
     		 * read configuration file 
     		 ************************************/
		ReadConfig(conf, DA_config);
	 		
		/*initialization ds's file*/
	
		fc.DA_dir = DA_config["DA_dir"].c_str();
		
		fc.in_dir = DA_config["input_dir"].c_str();
	
		fc.out_dir = DA_config["out_dir"].c_str();
		
		fc.lai_dir = DA_config["lai_dir"].c_str();
		
		fc.obs_f = DA_config["obs"].c_str();
		
		fc.landcover_file = DA_config["landcover"].c_str(); 		
		fc.BPLUT_file = DA_config["BPLUT_file"].c_str(); 		
			
		/*initialization ds's base configuration*/
		dc.K = atoi(DA_config["ensemble"].c_str());

		dc.grid_xres = atof(DA_config["res_lon"].c_str());	
		dc.grid_yres = atof(DA_config["res_lat"].c_str());	
    		dc.XSIZE     = int(360 / dc.grid_xres + 1e-5);
    		dc.YSIZE     = int(180 / dc.grid_yres + 1e-5);

		datetime start_date(DA_config["start_time"], DA_config["date_format"]);
		dc.start_date = start_date;
		datetime end_date(DA_config["end_time"], DA_config["date_format"]);
		dc.end_date = end_date;

		dc.DA_length = atoi(DA_config["DA_length"].c_str()); // assimilation time step     [minute]
		dc.DA_itr = atoi(DA_config["DA_itr"].c_str()); // assimilation step length   [minute]

		dc.nlag = atoi(DA_config["nlag"].c_str());

		dc.rho = atof(DA_config["inflation"].c_str());

		dc.Class = atoi(DA_config["Class"].c_str());
		dc.M = dc.Class * dc.nlag;

	}



};	

	
statevector<double> generate_x(int mpi_rank, daconfig dc,map<int,vector<pair<int,int> > >& regions_index)
{
    		double *x_b_b = (double*)malloc(dc.nlag * dc.Class * dc.K * sizeof(double));
    		statevector<double> x_b_lagMembers(dc.nlag, dc.Class, dc.K);
		if(mpi_rank == 0) {
			x_b_lagMembers.Initialization(dc,regions_index);
			//x_b_lagMembers = readvector("/wps/home/linxf/DA/GPP_DA_Result/diagnose/statevector.20060129000000.nc",dc); 
			x_b_b = x_b_lagMembers.expand();
			MPI_Bcast(x_b_b, dc.nlag * dc.Class * dc.K, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   		}
   		else{
    			MPI_Bcast(x_b_b, dc.nlag * dc.Class * dc.K, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   			reverse_expand(x_b_lagMembers,x_b_b, dc.nlag, dc.Class, dc.K);
		}
		//free(x_b_b);
		//x_b_p = NULL;
		return x_b_lagMembers;	
}	

mpi_config SAVE_MPI_Config(int mpi_size, int mpi_rank, int mpi_name_len, char * mpi_host_name)
{
	mpi_config mpic;
	mpic.size = mpi_size;
	mpic.rank = mpi_rank;
	mpic.name_len = mpi_name_len;
	//strcpy(mpic.host_name, mpi_host_name);

	return mpic;
}
#endif
