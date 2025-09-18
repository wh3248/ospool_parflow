rm -rf ./parflow_build small.sif
apptainer build --fakeroot small.sif small.def
if [ $? -ne 0 ]; then
    echo apptainer image build failed
    exit -1
else
    echo apptainer image buidl succeeded!
fi
