
INS_IBM=-I/share/soft/rhel6.2/netcdf-4.1.1Parallel/include
LD_IBM=-L/share/soft/rhel6.2/netcdf-4.1.1Parallel/lib

INS_HDF5=-I/share/soft/rhel6.2/hdf5-1.8.8Parallel/include/
LD_HDF5=-L/share/soft/rhel6.2/hdf5-1.8.8Parallel/lib

INS_GDAL=-I/data1/xubx/develop_tools_az/gdal/include/
LD_GDAL=-L/data1/xubx/develop_tools_az/gdal/lib/

#INS_EIGEN=-I/data1/xubx/develop_tools_az/eigen-eigen-36fd1ba04c12/eigen-eigen-36fd1ba04c12
INS_EIGEN=-I/wps/home/linxf/DA/DA_tools/eigen-eigen-5a0156e40feb

INS_BOOST=-I/wps/home/linxf/programFiles/boost/boost_az/include

GDAL = -I/wps/home/linxf/programFiles/gdal-2.2.0/gdal_az/include/ -L/wps/home/linxf/programFiles/gdal-2.2.0/gdal_az/lib/ -lgdal

NC_LIB=-lnetcdf_c++ -lnetcdf
HDF_LIB=-lhdf5_hl -lhdf5
FLAG=-std=c++0x -w
CXX=mpic++
CXX_IBM=mpic++ 


all:main


run:main
	bsub -I -n 20 -R "span[ptile=1]" mpirun -np 50 ./main
    
debug:main_d
	mpirun -np 2 ./main_d

main:main.cpp *.h
	${CXX_IBM} ${FLAG} -o $@ $< ${INS_EIGEN} ${INS_IBM} ${INS_HDF5}  ${LD_IBM} ${NC_LIB} ${LD_HDF5} ${HDF_LIB} ${GDAL}  ${INS_BOOST}

main_d:main.cpp *.h
	${CXX} -lnetcdf ${FLAG} -DDEBUG -g ${INS_local} ${LD_local} ${LIB} -o $@ $<


