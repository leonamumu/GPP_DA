/* *************************
 * Author: xbx1992
 * Created Time:  Tue 11  July 2015
 * LastModified:  Tue 11July 2015 04:00:27 PM CST
 * C File Name:  observationdata.h
 * ************************/
#ifndef OBERVATIONDATA_H
#  define OBERVATIONDATA_H

#include <iostream>
#include <cstring>
#include <map>
#include <vector>
#include "xg_datetime.h"
#include "configclass.h"
#include "gdal_priv.h"
#include "assert.h"
#include "cpl_conv.h"
#include "cpl_string.h"
#include "xg_code.h"
using namespace std;


class obs_c
{
	public:
	obs_c(){}
	OBS_DATA Initialize(const datetime & DA_start, const datetime & DA_end,daconfig dc, fileconfig fc)
	{
		OBS_DATA obs;
    		GDALAllRegister();
    		
		for(datetime DA_itr = DA_start; DA_itr < DA_end; DA_itr += timespan(0,0,0,0,dc.DA_itr)){
		
			observation buf_o;
        		double dx = dc.grid_xres;
        		double dy = dc.grid_yres;
		    
                	char GPP_str[128];//string buffer of command 
			sprintf(GPP_str, "%s/GPP_ecmwf_improvd_LAI_MOD12_SLM_xuglc_bylxf.%d.%d.%d.%d.tif", fc.obs_f, DA_itr.year, DA_itr.month, DA_itr.day, DA_itr.hour);

			printf("reading stations output :%s\n", GPP_str);
			datetime obs_time(DA_itr.year, DA_itr.month, DA_itr.day, DA_itr.hour);
			int ii = 0;
                        vector<double> rands(2701); //观测值空间分辨率为1°x1°，根据间隔5°采样一次观测值：构建随机向量73*37
			vector_randn_boost(rands, 2701, 0, 1, -1, 1);
    			
			GDALDataset *data = (GDALDataset*) GDALOpen(GPP_str, GA_ReadOnly);
    			assert(data != NULL);
    			GDALRasterBand *band = data -> GetRasterBand(1);
    			assert(band != NULL);

			int NaNs = band -> GetNoDataValue();
    			int YSize = band->GetYSize(), XSize = band->GetXSize();
    			for(int line = 0; line < YSize; line=line + 5){
        			vector<double> lines;
            			lines = ReadLine(band, line);
        			for(int i = 0; i < XSize; i =i + 5){
				    	buf_o.geosI = i + 1;
					buf_o.geosJ = line + 1;
					buf_o.alt   = 1;
					buf_o.tau = obs_time.tau();
					buf_o.mdm = 1;
					if(lines[i] == NaNs || lines[i] != lines[i] || (abs(lines[i]) < 0.000001))//观测为0不进入同化
						continue;
					buf_o.value = lines[i] + rands[ii];
					obs.push_back(buf_o);
					++ii;
				}
			}
    			GDALClose(data);
		}
			return obs;
	}	
	map<int, double> Get_R(OBS_DATA& obs)
	{
		   /********************
 		    * error covariance
 		    * *****************/ 
		    map<int, double> R_test;
		    for(size_t i = 0; i < obs.size(); ++i){
			    R_test[i] = pow(obs[i].mdm ,2);

                    }
		    return R_test;	
	}
};


OBS_INDEX build_index(const vector<observation>& obs){
    OBS_INDEX ans;
    for(size_t i = 0; i < obs.size(); ++i){
        ans[pair<int,int>(obs[i].geosI, obs[i].geosJ)].push_back(pair<double, int>(obs[i].tau, i));
    }
    for(OBS_INDEX::iterator it = ans.begin(); it!=ans.end(); ++it){
        sort(it->second.begin(), it->second.end(), cmp_1st<double, int>);
    }
    return ans;
}
#endif
