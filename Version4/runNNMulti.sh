#!/bin/bash
#SBATCH --job-name=FTD_NN             # Job name
#SBATCH --partition=bigjay            # Partition Name (Required)
#SBATCH --mail-type=ALL               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=b047m507@ku.edu   # Where to send mail	
#SBATCH --mail-type=FAIL              # Only send mail if the job fails
#SBATCH --ntasks=1                   # Run on a single CPU
#SBATCH --mem=16gb                     # Job memory request
#SBATCH --time=0-24:00:00             # Time limit days-hrs:min:sec
#SBATCH --output=FTD_NN_%j.log        # Standard output and error log

module list

module load root/6.32.08
#module load geant4

#sbatch runNN.sh


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

#root -l 'TMVA_MultiClass.C'

#root -l 'TMVA_MultiClass_v2.C'

root -l 'TMVA_MultiClass_v3.C'

#hadd TOTg100.root OUTg100*

echo 'Finished running FTD NN Test'

date

exit
