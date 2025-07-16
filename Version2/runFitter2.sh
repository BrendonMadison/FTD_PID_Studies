#!/bin/bash
#SBATCH --job-name=FTDFIT             # Job name
#SBATCH --partition=bigjay            # Partition Name (Required)
#SBATCH --mail-type=ALL               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=b047m507@ku.edu   # Where to send mail	
#SBATCH --mail-type=FAIL              # Only send mail if the job fails
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=8gb                    # Job memory request
#SBATCH --time=0-24:00:00             # Time limit days-hrs:min:sec
#SBATCH --output=FTDFIT_%j.log        # Standard output and error log

#To compile the executable:
#g++ -std=c++17 -O2 FitterExecutable_2.cpp -o FitterExecutable_2 `root-config --cflags --libs` -lMinuit

#To run:
#sbatch runFitter2.sh Output/gam_zpole_20k.root FitterGamTest.root

module list

module load root/6.32.08

INAME=$1
ONAME=$2

echo 'Input file name     '$INAME
echo 'Output file name    '$ONAME

echo 'Running LumiCal_MakePixel job as '
echo 'User '$USER

pwd
hostname
date

echo 'PATH '
echo $PATH
 
echo 'LD_LIBRARY_PATH'
echo $LD_LIBRARY_PATH

echo 'Run HGCAL_testbeam'

#MWORK=/panfs/pfs.local/work/wilson/b047m507/GEANT4/Install/HGCAL-build
MWORK=/kuhpc/work/wilson/b047m507/GEANT4/Install/HGCAL-build

cd ${MWORK}

./FitterExecutable_2 ${INAME} hits ${ONAME} ExTree

echo 'Finished running fitter!'

date

exit
