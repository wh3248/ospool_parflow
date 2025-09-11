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

RUN dnf install -y openssh-clients && dnf clean all

USER root
WORKDIR /parflow
COPY parflow/bin /parflow/bin
COPY parflow/lib /parflow/lib
COPY parflow/other_lib /parflow/other_list
COPY parflow/openmpi/4.1.0/gcc /usr/local/openmpi/4.1.0/gcc
COPY parflow/pmix /usr/share/pmix
ENV LD_LIBRARY_PATH=/parflow/lib:/usr/lib64:/parflow/other_list:${LIB_LIBRARY_PATH}
ENV PARFLOW_DIR=/parflow
