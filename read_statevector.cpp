/*Netcdf文件转换为TIF格式的文件*/
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <set>
#include <assert.h>
#include <map>
#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "ogrsf_frmts.h"
#include <netcdfcpp.h>
using namespace std;

int main(int argc, char *argv[])
{
		
		NcFile s_f(argv[1], NcFile::ReadOnly);
		if(!s_f.is_valid())
		{
			cout<<"couldn't open file"<<endl;
		//	return 0;
		}
		double temp[1 * 13 * 60];
	    s_f.get_var("apri_state")->get(temp, 1 ,13, 60);
		map< int, map<int, map<int, double> > > r_tmp;
		for(size_t i =0;i<1 * 13 * 60;++i)
			r_tmp[i / (60 * 13)][(i % (60 * 13)) / (60)][(i % (60 * 13)) % (60)] = temp[i];

		for(size_t i = 0;i < 1;++i)
		{
		  cout<<i<<" week"<<endl;

		  for(size_t j = 0;j < 13;++j)
		  {
		      vector<double> ensemble;
			  double sum = 0;
			  for(size_t z= 0;z< 60;z++)
			  {
			       ensemble.push_back(r_tmp[i][j][z]);
				   sum += r_tmp[i][j][z];
			  }
			  double mean = sum / 60;
			  printf("%lf\n",mean);
		  }
		}
	    return -1;
}
