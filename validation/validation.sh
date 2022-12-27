# User defined variables
CONDA_ENV_NAME="CompChem"

# Locate the conda script and source it
CONDA_INFO=`conda info | egrep "base environment :"`
CONDA_INFO=(${CONDA_INFO// / })

source "${CONDA_INFO[3]}/etc/profile.d/conda.sh"
conda activate $CONDA_ENV_NAME

# Copy the content of the source folder to a temporary production folder and enter it
mkdir tmp
cp -r ./source/* ./tmp
cd tmp

# Download parameters for DFTB+
wget https://dftb.org/fileadmin/DFTB/public/slako/3ob/3ob-3-1.tar.xz
unxz 3ob-3-1.tar.xz
tar -xvf 3ob-3-1.tar
mv 3ob-3-1 ./dftb
rm 3ob-3-1.tar

# Run DFTB+ calculation
cd ./dftb/water_opt

export OMP_NUM_THREADS=1
mpirun -np 1 dftb+ > water_opt.out

cd ../../

# Move the files to their final destination
# TO BE IMPLEMENTED

# Remove the temporary directory
cd ../
rm -rf tmp

