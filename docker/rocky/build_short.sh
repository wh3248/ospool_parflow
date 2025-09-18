rm -rf ./parflow_build short.sif
apptainer build short.sif short.def
if [ $? -ne 0 ]; then
    echo apptainer image build failed
    exit -1
else
    echo apptainer image build succeeded!
fi
