#!/bin/sh --login

#SBATCH --account=coastal
#SBATCH --job-name=a70_CHI_ATM_WAV2OCNv2.1
#SBATCH -q batch
#SBATCH --time=00:30:00
#SBATCH --ntasks=13
#SBATCH --mail-user=saeed.moghimi@noaa.gov
#SBATCH --mail-type=ALL
#SBATCH --output=a70_CHI_ATM_WAV2OCNv2.1.out.log
#SBATCH --error=a70_CHI_ATM_WAV2OCNv2.1.err.log


# -- load ENV variables
source   /scratch2/COASTAL/coastal/save/Saeed.Moghimi/models/NEMS/tests/new_nems_app/temp/ADC-WW3-NWM-NEMS/modulefiles/hera/ESMF_NUOPC

srun ./NEMS.x 


