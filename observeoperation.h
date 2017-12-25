/* *************************
 * Author: xbx1992
 * Created Time:  Tue 11  July 2015
 * LastModified:   08 Oct 2015 20:09 PM CST
 * C File Name:  observeoperation.h
 * FUNCTION:
 *      Calculate GPP from ECMWF data
 * USAGE:
 *      ./<this bin> <SSRD tif> <d2m tif> <Tmp tif> <LAI tif> <landcover tif> <BPLUT table> <output tif> <month> <day> <hour>
 * ************************/
#ifndef OBERVEOPERATION_H
#  define OBERVEOPERATION_H

#include <iostream>
#include <cstring>
#include <map>
#include <vector>
#include "common.h"
#include "configclass.h"
#include "xg_code.h"
using namespace std;

map<int, vector<double> > BPLUT;

int hour,day,month;


void init_BPLUT(const char *filename){
    FILE *fp = fopen(filename, "r");
    int c;
    vector<double> v(10);
    while(fscanf(fp, "%d", &c) != EOF){
        for(int i = 0; i < 10; ++i){
            fscanf(fp, "%lf", &v[i]);
        }
        BPLUT[c] = v;
    }
	fclose(fp);
}


void show_BPLUT( ){
    for(map<int, vector<double> > ::iterator it  = BPLUT.begin(); it != BPLUT.end(); ++it){
        printf("%d", it->first);
        for (size_t i = 0; i < it->second.size(); ++i){
            printf("\t%f", (it->second)[i]);
        }
        printf("\n");
    }
        printf("\n");
}

map<int,double> LUE_x;
double CosZs;

void cal_CosZs(const double lat,const double lon,const int hour,const int day,const int month,const double clumping){
        int days;
	double LAI_su,LAI_sh,Delta,Hsolar,hr,Lat_arc,LAIsu,LAIsh;
	// CosZs��ȡ����Сʱ���м�ʱ�̴������, h: unit/hour
	for(size_t m = 1;m < month;++m){
		days += monthlength[m];
	}
	days = days + day;
			// Caculating Solar Zenith
			Delta = 0.006918 - 0.399912*cos(2.0*3.1415926*days/365.0) + 0.070257*sin(2.0*3.1415926*days/365.0)- 0.006758*cos(4.0*3.1415926*days/365.0) + 0.000907*sin(4.0*3.1415926*days/365.0);	
		    // World time to local time
			hr = hour + lon/15;
			if(hr > 24)
				hr = hr - 24;
			else if(hr < 0)
				hr = 24 + hr;
			/*-------------------------*/
			Lat_arc = 3.1415926 * lat / 180.0; //unit:rad
			Hsolar = 2.0 * 3.1415926 * (hr -12.0) / 24.0; // if do not need convert time, you should use h instead of hr
			CosZs=cos(Lat_arc)*cos(Delta)*cos(Hsolar) + sin(Lat_arc)*sin(Delta);
}
		

double GPP(const vector<double> & _a, const vector<double> & Na,int y,int x){
    vector<double> a(_a);//ÕâÀï·µ»ØµÄ_aÊÇ¶à¸öÍ¼²ã£¬ÊÇÒ»¸övector,¶ÔÃ¿¸öÏñÔª¶à²ã²Ù×÷
    assert(a.size() == 5);
    double lat=90-y,lon=x-180;
    for(int i = 0; i < 5; ++i){
        if(abs(a[i] - Na[i]) < 1e-13){//no-data value exists
                return -1;
        }
    }
        if(abs(a[4] - 0) < 1e-13)//no-data value exists
                return -1;
    double SSRD = a[0]; // [ w/m2 ]
    double d2m = a[1] - 273.15; // how to calculate?
    double TMP = a[2] - 273.15; 
    double lai_o = a[3]*0.1; 
    int LC =int( a[4] + 0.1 );

    if(BPLUT.find(LC) ==  BPLUT.end()){//no vegetated landcover type
        return -1;
    }
 
    double t_max = BPLUT[LC][0]; // [ C ]
    double t_min = BPLUT[LC][1];
    double vpd_max = BPLUT[LC][2]; // [ KPa ]
    double vpd_min = BPLUT[LC][3]; 
    double clumping= BPLUT[LC][4]; 
    double stem_o  = BPLUT[LC][5]; 
    double albedo  = BPLUT[LC][6]; 
    double LUEmax = 0;
    LUEmax = LUE_x[LC] * LUE_def[LC];  
	//cal_CosZs(lat,lon,hour,day,month,clumping);  
    double laio=lai_o+stem_o;
    
    /*zhanghf cal*/
    
    double VPD=0.61078*(exp((17.27*TMP)/(237.3+TMP))-exp((17.27*d2m)/(237.3+d2m)));//kPa
    
    /*begin caculating fpar*/
    double fpar = 1-exp(-0.5*lai_o);
    //double fpar = 1 - albedo - exp(-0.5 * clumping * laio /CosZs);
    
    if(fpar <0)
    	return 0;
    
    double tmp_scalar = TMP > t_max ? 1 : ( TMP < t_min ? 0 :   (TMP - t_min) / (t_max - t_min) );
    double vpd_scalar = VPD > vpd_max ? 0 : ( VPD < vpd_min ? 1 :   (VPD- vpd_max) / (vpd_min - vpd_max) );
    double APAR = fpar * 0.5 * SSRD;
    
    if(APAR < 0)
	return 0;
    
    return (LUEmax * APAR * tmp_scalar * vpd_scalar);
    
    /*end*/
}

class obs_operation
{
	public:
	obs_operation(){}
	void Initialize(daconfig dc,fileconfig fc, mpi_config mpic, int mpi_rank,statevector<double> x_b_lagMembers,const vector<observation>& obs, const datetime & DA_start, const datetime & DA_end)
	{
        	for(size_t i = 0; i < dc.K; ++i)
		{
	   	if(size_t(mpi_rank % mpic.size) == (i % mpic.size)){
			//printf("PROCESS %d start handling ensemble-%.3lu on %s...\n", mpic.rank, i, mpic.host_name);
			printf("PROCESS %d start handling ensemble-%.3lu...\n", mpi_rank, i);

               		char cmd_str[128];//string buffer of command 
                	/*
                 	 * go into the running directory of ensemble-i
                 	 */
		
    		 	/*save LUEmax*/
        		for(size_t j = 0; j < dc.Class; ++j)
				LUE_x[j + 1] = x_b_lagMembers[0][j][i];
              		
			sprintf(cmd_str, "mkdir %s/ensemble-%.3lu", fc.DA_dir, i );
                	cmd(cmd_str);
			
    			for(datetime DA_itr = DA_start; DA_itr < DA_end; DA_itr += timespan(0,0,0,0,dc.DA_itr)){

     			vector<string> infs; 

		 	/*write ssrd file*/	
               		char ssrd_str[128];//string buffer of command 
               		sprintf(ssrd_str, "%s/ssrd/ec-ei-fc012up2tr3-sfc-glb100x100-ssrd_%d%.2lu%.2lu_00p03.hdf.%d.tif.trans.tif", fc.in_dir, DA_itr.year, DA_itr.month, DA_itr.day, DA_itr.hour);
			infs.push_back(ssrd_str);
			
		 	/*write d2m file*/	
               		char d2m_str[128];//string buffer of command 
               		sprintf(d2m_str, "%s/d2m/ec-ei-fc012up2tr3-sfc-glb100x100-d2m_%d%.2lu%.2lu_00p03.hdf.%d.tif.trans.tif", fc.in_dir, DA_itr.year, DA_itr.month, DA_itr.day, DA_itr.hour);
			infs.push_back(d2m_str);
			
		 	/*write Tmp file*/	
               		char Tmp_str[128];//string buffer of command 
               		sprintf(Tmp_str, "%s/t2m/ec-ei-fc012up2tr3-sfc-glb100x100-t2m_%d%.2lu%.2lu_00p03.hdf.%d.tif.trans.tif", fc.in_dir, DA_itr.year, DA_itr.month, DA_itr.day, DA_itr.hour);
			infs.push_back(Tmp_str);

		 	/*write LAI file*/	
               		char LAI_str[128];//string buffer of command 
               		sprintf(LAI_str, "%s/offset.global_30s_%d_%d_trans.tif.1x1.tif", fc.lai_dir, DA_itr.year, DA_itr.month);
			infs.push_back(LAI_str);
			
		 	/*write landcover file*/	
               		char landcover_str[128];//string buffer of command 
               		sprintf(landcover_str, "%s",fc.landcover_file);
			infs.push_back(landcover_str);

			
		 	/*write BPLUT file*/	
               		char BPLUT_str[128];//string buffer of command 
               		sprintf(BPLUT_str, "%s",fc.BPLUT_file);
			init_BPLUT(BPLUT_str);	
			

		 	/*write output file*/	
               		char output_str[128];//string buffer of command 
               		sprintf(output_str, "%s/ensemble-%.3lu/GPP_ecmwf_improvd_LAI_MOD12_SLM_bylxfandxubx.%d.%d.%d.%d.tif",fc.out_dir, i, DA_itr.year, DA_itr.month, DA_itr.day, DA_itr.hour);
			string o_s = output_str; 
			
			//assert(infs.size() == 4);

			GDALAllRegister();
    			
			vector<string> inp;
    			
			printf("Reducing (\n");
    			for(int f = 0; f < infs.size(); ++f){
        			printf("%s\n",infs[f].c_str());
        			inp.push_back(infs[f]);
    			}
    			printf(") into %s\n",o_s.c_str() );
    
                	/* 
                 	 * run  model 
                 	 */
       	       		cout<<"start begining run ensemble-"<<i<<endl;
               		
			sprintf(cmd_str, "(cd %s/ensemble-%.3lu;)",fc.DA_dir, i);
                	cmd(cmd_str);
			
			Reduce(GPP, inp, o_s, GDT_Float64, -1);

			sprintf(cmd_str, "(cd ..;)");
                	cmd(cmd_str);
			
			}
			printf("PROCESS %d FINISHED...\n", i);
			break;
			}
			
        	}// END OF GENERATING y_b for each ensemble member
		//printf("PROCESS %d on %s FINISHED...\n", mpi_rank, mpic.host_name);
	}


	void Run_posterior(daconfig dc, fileconfig fc, MatrixXf x_b_bar,OBS_DATA obs, datetime DA_start, datetime DA_end)
	{
			
		    	/**********************
		     	 *  SAVE scal_x_a into file
		     	 **********************/
		    	char buf_str[128];//string buffer of command 
		    	sprintf(buf_str, "%s/posterior_scaling_factor.%s.geos.2x25", fc.DA_dir, DA_start.str("YYYYMMDDhhmmss").c_str());
		     	/*************************************************
		     	*
		     	* START POSTERIOR RUNNING
		     	*
		     	* ***********************************************/
		    	printf("START POSTERIOR RUNNING...\n");
		    	char cmd_str[128];//string buffer of command 
		    
		    	/*
		     	* go into the running directory of ensemble-i
		     	*/
		    	sprintf(cmd_str, "mkdir %s/posterior", fc.DA_dir );
		    	cmd(cmd_str);

    		  	/*save LUEmax*/
        		for(size_t j = 0; j < dc.Class; ++j)
				LUE_x[j + 1] = x_b_bar(j,0);
			
    			for(datetime DA_itr = DA_start; DA_itr < DA_end; DA_itr += timespan(0,0,0,0,dc.DA_itr)){

     			vector<string> infs; 

		 	/*write ssrd file*/	
               		char ssrd_str[128];//string buffer of command 
               		sprintf(ssrd_str, "%s/ssrd/ec-ei-fc012up2tr3-sfc-glb100x100-ssrd_%d%.2lu%.2lu_00p03.hdf.%d.tif.trans.tif", fc.in_dir, DA_itr.year, DA_itr.month, DA_itr.day, DA_itr.hour);
			infs.push_back(ssrd_str);
			
		 	/*write d2m file*/	
               		char d2m_str[128];//string buffer of command 
               		sprintf(d2m_str, "%s/d2m/ec-ei-fc012up2tr3-sfc-glb100x100-d2m_%d%.2lu%.2lu_00p03.hdf.%d.tif.trans.tif", fc.in_dir, DA_itr.year, DA_itr.month, DA_itr.day, DA_itr.hour);
			infs.push_back(d2m_str);
			
		 	/*write Tmp file*/	
               		char Tmp_str[128];//string buffer of command 
               		sprintf(Tmp_str, "%s/t2m/ec-ei-fc012up2tr3-sfc-glb100x100-t2m_%d%.2lu%.2lu_00p03.hdf.%d.tif.trans.tif", fc.in_dir, DA_itr.year, DA_itr.month, DA_itr.day, DA_itr.hour);
			infs.push_back(Tmp_str);

		 	/*write LAI file*/	
               		char LAI_str[128];//string buffer of command 
               		sprintf(LAI_str, "%s/offset.global_30s_%d_%d_trans.tif.1x1.tif", fc.lai_dir, DA_itr.year, DA_itr.month);
			infs.push_back(LAI_str);
			
		 	/*write landcover file*/	
               		char landcover_str[128];//string buffer of command 
               		sprintf(landcover_str, "%s",fc.landcover_file);
			infs.push_back(landcover_str);

			
		 	/*write BPLUT file*/	
               		char BPLUT_str[128];//string buffer of command 
               		sprintf(BPLUT_str, "BPLUT_UMD.csv");
			init_BPLUT(BPLUT_str);	
		 	

		 	/*write output file*/	
               		char output_str[128];//string buffer of command 
               		sprintf(output_str, "%s/posterior/GPP_ecmwf_improvd_LAI_MOD12_SLM_bylxfandxubx.%d.%d.%d.%d.tif",fc.out_dir, DA_itr.year, DA_itr.month, DA_itr.day, DA_itr.hour);
		      
			string o_s = output_str; 
			
			assert(infs.size() == 5);

			GDALAllRegister();
    			
			vector<string> inp;
    			
			printf("Reducing (\n");
    			for(int f = 0; f < infs.size(); ++f){
        			printf("%s\n",infs[f].c_str());
        			inp.push_back(infs[f]);
    			}
    			printf(") into %s\n",o_s.c_str() );
    
                	/* 
                 	 * run  model 
                 	 */
               		
			sprintf(cmd_str, "(cd %s/posterior;)",fc.DA_dir);
                	cmd(cmd_str);
			
			Reduce(GPP, inp, o_s, GDT_Float64, -1);

			sprintf(cmd_str, "(cd ..;)");
                	cmd(cmd_str);
		}
	}
};

#endif
