DOCKER=podman
DOCKER_REPOSITORY=hydroframeml

IMAGE_VERSON=0.1
IMAGE_NAME=ospool
${DOCKER} image rm --f ${IMAGE_NAME}:${IMAGE_VERSON}
dangling_images=$(${DOCKER} images -f "dangling=true" -q)
if [ -z "$dangling_images" ]
then
    echo "No dangling images"
else
    echo "Remove dangling images"
    ${DOCKER} rmi --force $(${DOCKER} images -f "dangling=true" -q)
fi
${DOCKER} build -t ${IMAGE_NAME}:${IMAGE_VERSON} -f Dockerfile .
if [ $? -ne 0 ]; then
    echo Docker image build failed
    exit -1
fi
${DOCKER} tag ${IMAGE_NAME}:${IMAGE_VERSON} docker.io/${DOCKER_REPOSITORY}/${IMAGE_NAME}:${IMAGE_VERSON}
if [ $? -ne 0 ]; then
    echo Docker tag creation failed
    exit -1
fi
${DOCKER} push docker.io/${DOCKER_REPOSITORY}/${IMAGE_NAME}:${IMAGE_VERSON}
if [ $? -ne 0 ]; then
    echo Docker tag push failed
    exit -1
fi
