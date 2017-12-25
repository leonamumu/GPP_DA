/* *************************
 * Author: xbx1992
 * Created Time:  Tue 11  July 2015
 * LastModified:  Tue 11July 2015 04:00:27 PM CST
 * C File Name:  sample.h
 * ************************/
#ifndef SAMPLE_H
#  define SAMPLE_H

#include "bpch.h"
#include "observationdata.h"
#include "xg_datetime.h"
#include "configclass.h"

int lookup_obs(const OBS_INDEX& index, int x, int y, double tau, double max_diff = 0.25,int i=0){
    
    printf("look for tau = %lf...", tau);
    OBS_INDEX::const_iterator it = index.find(pair<int,int>(x,y));
    if(it != index.end()){
        int a = 0, b = it->second.size() - 1;
        while(a < b - 1){
            int mid = (a + b) / 2;
            if( it->second[mid].first > tau){
                b = mid;
            }
            else 
                a = mid;
        }
        double da = fabs(tau - it->second[a].first);
        double db = fabs(tau - it->second[b].first);
	//if(i==0)
	//	printf("x:%d,y:%d,a:%d,b:%d,da:%lf,db:%lf,tau:%lf\n",x,y,a,b,da,db,tau);
        if(da <= db && da <= max_diff){
            printf(" index = %d\n",it->second[a].second );
            return it->second[a].second;
        }
        else if(da > db && db <= max_diff){
            printf(" index = %d\n",it->second[b].second );
            return it->second[b].second;
        }
	else {
            printf(" index = -1\n");
            return -1;
        }
    }
    else {
        printf(" index = -1\n");
        return -1;
    }
}

double GetTifPixelValue(const char* tif, double x, double y);

class Sample{
	public:
	Sample(){}
	void Initialize(OBS_DATA& obs, datetime DA_start,datetime DA_end, daconfig dc, fileconfig fc)
	{
	    	vector<mod_data> mod;
	    	mod_data buf_mod;
	    	vector<size_t> index_v(obs.size(),-1); 
    		
    		GDALAllRegister();
		clock_t time_f = clock(),time_s;
            	
		for(size_t i = 0; i < dc.K; ++i) {
                	char cmd_str[128];//string buffer of command 
                	/*
                 	 * read modeled observation 
                 	 */
			
		for(datetime DA_itr = DA_start; DA_itr < DA_end; DA_itr += timespan(0,0,0,0,dc.DA_itr)){
			sprintf(cmd_str, "%s/ensemble-%.3lu/GPP_ecmwf_improvd_LAI_MOD12_SLM_bylxfandxubx.%d.%d.%d.%d.tif", fc.DA_dir, i, DA_itr.year, DA_itr.month, DA_itr.day, DA_itr.hour);
                    	printf("reading stations output :%s\n", cmd_str);
			datetime obs_time(DA_itr.year, DA_itr.month, DA_itr.day, DA_itr.hour);
			
			GDALDataset *data = (GDALDataset*) GDALOpen(cmd_str, GA_ReadOnly);
    			assert(data != NULL);
    			GDALRasterBand *band = data -> GetRasterBand(1);
    			assert(band != NULL);

			int NaNs = band -> GetNoDataValue();
    			int YSize = band->GetYSize(), XSize = band->GetXSize();
    			for(int line = 0; line < YSize; line=line + 10){
        			vector<double> lines;
            			lines = ReadLine(band, line);
        			for(int ii = 0; ii < XSize; ii =ii + 10){
				    	buf_mod.geosI = ii + 1;
					buf_mod.geosJ = line + 1;
					buf_mod.tau = obs_time.tau();
					buf_mod.value = lines[ii];
					if(buf_mod.value == NaNs || buf_mod.value != buf_mod.value)
						continue;
					//if(i==0 )
					//	cout<<"sample-"<<i<<","<<buf_mod.geosI<<","<<buf_mod.geosJ<<","<<buf_mod.value<<endl;
					mod.push_back(buf_mod);
				}
			}
    			GDALClose(data);
		}	
			printf("finish read!"); 
    		
			time_s = Cal_time(time_f); 
			
			OBS_INDEX ans;
			for(size_t j = 0;j < mod.size();++j){
				ans[pair<int,int>(mod[j].geosI,mod[j].geosJ)].push_back(pair<double,int>(mod[j].tau,j));
			}
			

			for(OBS_INDEX::iterator it = ans.begin();it!=ans.end();++it){
				sort(it->second.begin(), it->second.end(), cmp_1st<double, int>);	
			}
			
			printf("finish short!"); 
    			time_f = Cal_time(time_s);
			
			for(size_t j =0; j < obs.size(); ++j){
				if(obs[j].tau <= (DA_start + timespan(0,0,0,0,dc.DA_length)).tau())
				{
				if(i == 0)
					index_v[j] = lookup_obs(ans,obs[j].geosI,obs[j].geosJ,obs[j].tau, 0.25, 0);
				if(index_v[j] >= 0 && (mod[index_v[j]].value < -0.00001 || mod[index_v[j]].value == 0 || mod[index_v[j]].value > 0.00001))
					if(abs(mod[index_v[j]].value) < 1e-6)
						 obs[j].modeled_v.push_back(0);
					else
					     obs[j].modeled_v.push_back(mod[index_v[j]].value);
				}
				else
					break;	
			}
			mod.clear();
			
			vector<mod_data>(mod).swap(mod);
			
			printf("finish sample!"); 
    			time_s = Cal_time(time_f);
            	}

            	/*
             	 * DEBUG OUTPUT
             	 *
             	 * content of obs
             	 */
            	for(size_t i = 0; i < obs.size(); ++i){
			printf("obs-%lu:\n(%d, %d, %d, [%lf])-> %.10lf\n", i, obs[i].geosI, obs[i].geosJ, obs[i].alt,obs[i].tau, obs[i].value);
                	cout << obs[i].modeled_v << endl;
            	}
	
	}	
};



#endif
