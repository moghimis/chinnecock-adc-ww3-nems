rm -f field_*.nc
rm PET*.ESMF_LogFile
cp -fv /scratch4/COASTAL/coastal/save/Saeed.Moghimi/models/NEMS/tests/NSEModel_try_err_branches/HWRF2ADC_test02/NEMS/exe/NEMS.x .
source /scratch4/COASTAL/coastal/save/Saeed.Moghimi/models/NEMS/tests/NSEModel_try_err_branches/HWRF2ADC_test02/modulefiles/theia/fv3-saeed
ln -sfv   NEMS.x  esm1
mpirun -np 5 ./esm1




