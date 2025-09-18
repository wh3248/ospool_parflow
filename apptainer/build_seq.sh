rm -rf parflow_seq.sif
apptainer build parflow_seq.sif parflow_seq.def
if [ $? -ne 0 ]; then
    echo apptainer parflow_seq.sif build failed
    exit -1
else
    echo apptainer parflow_seq.sif build succeeded!
fi
