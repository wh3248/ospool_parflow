rm -rf ./parflow_build ub.sif
apptainer build ub.sif ub.def
if [ $? -ne 0 ]; then
    echo apptainer image build failed
    exit -1
else
    echo apptainer image build succeeded!
fi
