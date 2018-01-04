/* *************************
 * Author: xbx1992
 * Created Time:  Tue 12  July 2015
 * LastModified:  Tue 11July 2015 04:00:27 PM CST
 * C File Name:  ensrf.h
 * ************************/
#ifndef ENSRF_H
#  define ENSRF_H

#include <algorithm>
#include <Eigen/Dense>
#include <map>
#include "xg_debug.h"
#include "xg_datetime.h"
#include "observationdata.h"
#include "observeoperation.h"
#include "configclass.h"
using namespace Eigen;
class optimize
{
	private:
	    MatrixXf P_b;
            MatrixXf Y_b, y_o, y_b, y_b_bar, y_b_bar_met;
	    MatrixXf x_b_inflation,PHt,HPHt,K_gain;
	    MatrixXf Coef;
	    map<int, double> HPHR,alpha;
	    int index_x;
	    double L;
	    vector<int> y_b_nan;
	public:
            	MatrixXf X_b,x_b_bar;
		optimize(){}
		void Initialize(statevector<double> x_b_lagMembers, OBS_DATA obs,daconfig dc)
		{
		     	/*******************
 		      	 *	resize matrix
 		      	 ******************/
			 x_b_bar.resize(dc.M,1),x_b_inflation.resize(dc.M,dc.K),Coef.resize(dc.M,1),y_o.resize(obs.size(),1),y_b.resize(obs.size(),dc.K),y_b_bar.resize(obs.size(),1),y_b_bar_met.resize(obs.size(),1);
	     		/*
             		 * CALCULATE X_b
             		 */
	    		X_b = x_b_lagMembers.transfer(); 
	    		
			//debug(X_b);
	    		x_b_bar = X_b.rowwise().mean();			    
	    		for(size_t i = 0; i < X_b.cols(); ++i)
            		{
				
				for(size_t j = 0; j < x_b_bar.rows();++j){
					X_b(j,i) = X_b(j,i) - x_b_bar(j,0); 
				} 
	    		}  
	    		//debug(X_b);
			/*prior converience*/
	    		P_b =  (X_b * X_b.transpose()) / (dc.K - 1);
            		/*
             		 * BUILD INDEX of obs with sufficient y_b
             		 */
				for(size_t i = 0; i < obs.size(); ++i){
                	
				y_o(i,0) = obs[i].value;
				if(obs[i].modeled_v.size() == dc.K){
		     			for(size_t j=0;j < dc.K;++j)
						y_b(i,j) = obs[i].modeled_v[j];
               	}
				else
					y_b_nan.push_back(i);
	    		}	
	    		y_b_bar_met = y_b_bar = y_b.rowwise().mean();
	    		debug(y_b_bar);
	   		Y_b = y_b; 
			for(size_t i = 0;i < Y_b.rows(); ++i)
            			for(size_t j = 0; j < Y_b.cols(); ++j)
					Y_b(i,j) = Y_b(i,j) - y_b_bar(i,0); 

			
		    	K_gain = MatrixXf::Zero(dc.M, obs.size());

		}
		
		statevector<double> Run(daconfig dc, OBS_DATA obs,map<int, double> R,MatrixXf& x_b_bar_temp, REGIONS_MAP regions,map<int, vector<pair<int, int> > >& regions_index, double tau_start,fileconfig fc,const datetime & DA_start, const datetime & DA_end)
		{
		    statevector<double> x_b_lagMembers(dc.nlag,dc.Class,dc.K);
				
		     /******************
 		     * Pre for obs
 		     * *****************/ 		
		    x_b_bar_temp = x_b_bar;
		     for(size_t n = 0;n < obs.size(); ++n){
		     	   /******************
                     	    * STEP 1 check Y_b
                            * ***************/
				 vector<int>::iterator result = find( y_b_nan.begin( ), y_b_nan.end( ), n ); //查找n
				 if(result!= y_b_nan.end())
				 {
					y_b_nan.erase(result); 
					continue;
				 }
				if( Y_b(n,0) == -1 || (abs(y_b_bar_met(n,0)) > 1e6) || obs[n].value != obs[n].value)
				{
					continue;
				}
		     	   /******************
                     	    * STEP 2 calculate PHt
                            * ***************/
			    PHt = X_b * (Y_b.row(n).transpose());
			   for(size_t i = 0;i < PHt.cols();++i)
				for(size_t j = 0; j < PHt.rows();++j)
					PHt(j,i) = PHt(j,i) / (dc.K - 1);
    		   	   /******************************
     		            * generate Coef
     		            ****************************/
		            debug("before PHt");
			    int index_r = regions[obs[n].geosJ][obs[n].geosI];
			    debug(index_r);
			    if(index_r > 12 || index_r == 11 || index_r == 0)
			    	continue;	    
			    int a_x,a_y,b_x,b_y;
			    int index_v = regions_index[index_r].size() / 2;
			    debug(index_v);
			    a_x =  regions_index[index_r][index_v].first; 
			    a_y =  regions_index[index_r][index_v].second;
			    debug(obs[n].geosJ);
			    debug(obs[n].geosI);
				for(size_t i = 0; i < dc.M; ++i) 
			    {
					Coef(i,0) = 0;
			    		
					index_v = regions_index[i % dc.Class].size() / 2;
					b_x =  regions_index[i % dc.Class][index_v].first; 
					b_y = regions_index[i % dc.Class][index_v].second;
					if(abs(b_x - a_x) > 5 || abs(b_y - a_y) > 5)
			   			continue;
					L = sqrt(pow(a_x - b_x,2) +  pow(a_y - b_y,2));
					Coef(i,0) = exp(-L/7);
			   }
			   Coef(index_r,0) = 1;
			   debug(index_v);
			   //debug(Coef);
			    
		     	   /******************
                            * STEP 3 calculate HPHR
                            * ***************/
			    MatrixXf a(Y_b.row(n));
				
				HPHt = a * a.transpose();
				HPHR[n] = HPHt(0,0) / (dc.K - 1) + R[n];
			    debug(HPHt);
			    debug(HPHR[n]);
		     	   /******************
                            * STEP 4 calculate K_gain and post_x_b
			    *
			    * change plus and minus
                            * ***************/
			   for(size_t i = 0;i < K_gain.rows();++i)
			   {
					
		           	//debug(x_b_bar_temp(i,0));
		           	double tt = x_b_bar_temp(i,0);
				K_gain(i,n) = PHt(i,0) / HPHR[n];
				double temp = x_b_bar_temp(i,0);
				debug(temp);
				//x_b_bar_temp(i,0) = x_b_bar_temp(i,0) + K_gain(i,n) * (y_o(n,0) - y_b_bar(n,0));
				x_b_bar_temp(i,0) = x_b_bar_temp(i,0) +  Coef(i,0) * K_gain(i,n) * (y_o(n,0) - y_b_bar(n,0));
				if(index_r == i)
						printf("opt_obs-%d:%d,%lf,%lf,%lf,%lf,%lf\n",n,i,temp,x_b_bar_temp(i,0),K_gain(i,n) ,y_o(n,0),y_b_bar(n,0));
				debug(K_gain(i,n));
				debug(Coef(i,0));
		        debug(x_b_bar_temp(i,0));
			   }
			   alpha[index_r] = 1.0 / (1.0 + sqrt(R[n] / HPHR[n])); 
			   //debug(alpha);
			   /*****************
 			   * update the deviations from the mean state vector
 			   * ***************/ 
			   debug(alpha[index_r]);
			   for(size_t k = 0;k < dc.K;++k)
                           {
				for(size_t i = 0;i < dc.M;++i)
				{
				X_b(i,k) = X_b(i,k) - Coef(i,0) * alpha[index_r]  * K_gain(i,n) *(Y_b(n,k));
				//debug(K_gain(i,n));
				}
				//debug(Y_b(n,k));
			   } 
				   //debug(X_b);
				   /***************
				   * updat the ensemble of sampled CO2 concentrations
				   * **************/  
				   int length_L = 0;
				   for(size_t m = n+1;m < obs.size();++m)
				   {
					L = sqrt(pow(obs[n].geosJ - obs[m].geosJ,2) +  pow(obs[n].geosI - obs[m].geosI,2));
					(regions[obs[m].geosJ-1][obs[m].geosI-1] < 12 && regions[obs[m].geosJ-1][obs[m].geosI-1] != 0)?:length_L=7;length_L=10;
					MatrixXf a2(Y_b.row(n)),b2(Y_b.row(m)),Y_b_term(1,1);
					Y_b_term = a2 * b2.transpose();
					double fac =(Y_b_term(0,0) / HPHR[n]) / ( dc.K - 1);
					debug(fac);
					//if(index_r == m)
					{
					    //y_b_bar(m,0) = y_b_bar(m,0) + fac * (y_o(n,0) - y_b_bar(n,0));
					    y_b_bar(m,0) = y_b_bar(m,0) + exp(-L/length_L) * fac * (y_o(n,0) - y_b_bar(n,0));
					
					for(size_t k = 0;k < dc.K;++k)
						//Y_b(m,k) = Y_b(m,k) - alpha[index_r] * fac *(Y_b(n,k));
						Y_b(m,k) = Y_b(m,k) - exp(-L/length_L) * alpha[index_r] * fac *(Y_b(n,k));
				   
					}
				   }
				   //debug(y_b_bar);
				   //debug(Y_b);
				   for(size_t m = 0;m < n + 1;++m)
				   {
					L = sqrt(pow(obs[n].geosJ - obs[m].geosJ,2) +  pow(obs[n].geosI - obs[m].geosI,2));
					(regions[obs[m].geosJ-1][obs[m].geosI-1] < 12 && regions[obs[m].geosJ-1][obs[m].geosI-1] != 0)?:length_L=7;length_L=10;
					MatrixXf a2(Y_b.row(n)),b2(Y_b.row(m)),Y_b_term(1,1);
					 Y_b_term = a2 * b2.transpose();
					double fac = (Y_b_term(0,0) / HPHR[n]) / (dc.K - 1);
					//if(index_r == m)
					{
					//y_b_bar(m,0) = y_b_bar(m,0) + fac * (y_o(n,0) - y_b_bar(n,0));
					y_b_bar(m,0) = y_b_bar(m,0) + exp(-L/length_L) * fac * (y_o(n,0) - y_b_bar(n,0));
					
					for(size_t k = 0;k < dc.K;++k)
						//Y_b(m,k) = Y_b(m,k) - alpha[index_r] * fac *(Y_b(n,k));
						Y_b(m,k) = Y_b(m,k) - exp(-L/length_L) * alpha[index_r] * fac *(Y_b(n,k));
					}
				   }
				   debug(n);
			    }
			    /*****************
			     * STEP 4 x_b
			     ****************/
			    /*****************
			     * STEP 4 inflation x_b
			     ****************/
			    for(size_t i = 0;i < dc.M; ++i)
				    for(size_t j = 0;j < dc.K; ++j)
				//X_b(i,j) = x_b_bar_temp(i,0) + dc.rho * X_b(i,j); 
				X_b(i,j) = x_b_bar_temp(i,0) + (1 + Coef(i,0) * (dc.rho - 1)) * X_b(i,j);
				//debug(x_b);
				debug(x_b_bar_temp);
				debug(X_b);
				
		    /*****************
 		     *
                     * STEP 5 print Y_b_bar and obs
                     * ***************/
		    for(size_t i = 0;i < y_b_bar.rows(); ++i)
		    	for(size_t j = 0;j < dc.K; ++j)
			    Y_b(i,j) = y_b_bar(i,0) + Y_b(i,j); 
		    y_b_bar = Y_b.rowwise().mean();	    
		
		    for(size_t i = 0;i < y_b_bar.rows(); ++i)
			    printf("%d %d %lf %.10lf %.10lf %.10lf\n",obs[i].geosI,obs[i].geosJ,obs[i].tau,obs[i].value,y_b_bar(i,0),y_b_bar_met(i,0));
		    x_b_lagMembers = reverse_transfer(X_b,dc.nlag,dc.Class);	    	 
	    
	
		    return x_b_lagMembers;   
		}
};
	
#endif
