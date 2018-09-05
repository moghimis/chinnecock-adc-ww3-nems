
ncf1 = '/disks/NASARCHIVE/saeed_moghimi/codes/adc_cap/07-adcirc_esm_v10/v6_wave_file_ike/chinnecock_ww3_file_esm1_no_wave/fort.63.nc'
ncf2 = '/disks/NASARCHIVE/saeed_moghimi/codes/adc_cap/07-adcirc_esm_v10/v6_wave_file_ike/chinnecock_ww3_file_esm1_wave/fort.63.nc'


import netCDF4

nc1 = netCDF4.Dataset(ncf1)
ncv1 = nc1.variables
z1 = ncv1['zeta'][:]
#
nc2 = netCDF4.Dataset(ncf2)
ncv2 = nc2.variables
z2 = ncv2['zeta'][:]
#

dz = z2-z1
plot(dz[-1,:])




