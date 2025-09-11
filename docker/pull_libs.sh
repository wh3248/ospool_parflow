#######
# Pull parflow and libs from verde to use for docker build
#######
mkdir -p parflow
mkdir -p parflow/bin
mkdir -p parflow/lib
mkdir -p parflow/other_lib
mkdir -p parflow/openmpi
mkdir -p parflow/pmix
#scp hmei-hydro@verde.princeton.edu:/home/SHARED/software/parflow/7c8e7f0/bin/* parflow/bin
#scp hmei-hydro@verde.princeton.edu:/home/SHARED/software/parflow/7c8e7f0/lib/* parflow/lib
#scp hmei-hydro@verde.princeton.edu:/usr/local/netcdf/gcc/hdf5-1.10.6/openmpi-4.1.0/4.7.4/lib64/libnetcdf.so.18* parflow/other_lib
#scp hmei-hydro@verde.princeton.edu:/usr/local/hdf5/gcc/openmpi-4.1.0/1.10.6/lib64/libhdf5_hl.so.100* parflow/other_lib
#scp hmei-hydro@verde.princeton.edu:/usr/local/hdf5/gcc/openmpi-4.1.0/1.10.6/lib64/libhdf5.so.103* parflow/other_lib
#scp hmei-hydro@verde.princeton.edu:/usr/local/openmpi/4.1.6/gcc/lib64/libmpi_usempif08.so* parflow/other_lib
#scp hmei-hydro@verde.princeton.edu:/usr/local/openmpi/4.1.6/gcc/lib64/libmpi_usempi_ignore_tkr.so* parflow/other_lib
#scp hmei-hydro@verde.princeton.edu:/usr/local/openmpi/4.1.0/gcc/lib64/libmpi_mpifh.so.40* parflow/other_lib
#scp hmei-hydro@verde.princeton.edu:/usr/lib64/libgfortran.so.5* parflow/other_lib
#scp hmei-hydro@verde.princeton.edu:/usr/lib64/libquadmath.so.0* parflow/other_lib
#scp hmei-hydro@verde.princeton.edu:/usr/local/openmpi/4.1.0/gcc/lib64/libmpi.so.40* parflow/other_lib
#scp hmei-hydro@verde.princeton.edu:/usr/lib64/libsz.so.2* parflow/other_lib
#scp hmei-hydro@verde.princeton.edu:/usr/local/openmpi/4.1.0/gcc/lib64/libopen-rte.so.40* parflow/other_lib
#scp hmei-hydro@verde.princeton.edu:/usr/local/openmpi/4.1.0/gcc/lib64/libopen-pal.so.40* parflow/other_lib
#scp hmei-hydro@verde.princeton.edu:/usr/lib64/libhwloc.so.15* parflow/other_lib
#scp hmei-hydro@verde.princeton.edu:/usr/lib64/libevent_core-2.1.so.6* parflow/other_lib
#scp hmei-hydro@verde.princeton.edu:/usr/lib64/libevent_pthreads-2.1.so.6* parflow/other_lib
#scp hmei-hydro@verde.princeton.edu:/usr/lib64/libaec.so.0* parflow/other_lib
#scp hmei-hydro@verde.princeton.edu:/usr/lib64/libcrypto.so* parflow/other_lib
#scp hmei-hydro@verde.princeton.edu:/usr/lib64/libpmix.so.2* parflow/other_lib
#scp -r hmei-hydro@verde.princeton.edu:/usr/local/openmpi/4.1.0 parflow/openmpi
#scp -r hmei-hydro@verde.princeton.edu:/usr/share/pmix/* parflow/pmix
