DOCKER=apptainer
DOCKER_REPOSITORY=hydroframeml

IMAGE_VERSON=0.2
IMAGE_NAME=ospool
${DOCKER} build short.sif short.def
if [ $? -ne 0 ]; then
    echo apptainer image build failed
    exit -1
fi
