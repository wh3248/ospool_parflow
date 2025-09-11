# -----------------------------------------------------------------------------
# Build options (* = default)
# -----------------------------------------------------------------------------
#
#  BASE_IMAGE
#    *  Red Hat Version 9
#
#  BUILD_TYPE
#    *  Release
#
# -----------------------------------------------------------------------------

FROM docker.io/rockylinux/rockylinux:9


RUN dnf install -y && gcc-gfortran && \
    dnf clean all


#RUN dnf install -y netcdf && \
#    dnf clean all

WORKDIR /parflow
COPY parflow/bin /parflow/bin
COPY parflow/lib /parflow/lib
COPY parflow/netcdf/* /lib64
COPY parflow/hdf5/* /lib64
COPY parflow/other/* /lib64
COPY parflow/openmpi/* /lib64
ENV LD_LIBRARY_PATH=/parflow/lib:/usr/lib64:${LIB_LIBRARY_PATH}
ENV PARFLOW_DIR=/parflow


USER root
