# Script for plotting unstructured WW3 output files in netCDF format, created by ww3_ounf
# Andre van der Westhuysen, 08/21/13
#
import matplotlib
matplotlib.use('Agg',warn=False)

import netCDF4
import numpy as np
import numpy.ma as ma
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import scipy.io as sio
import datetime
from scipy.io.matlab.mio5 import varmats_from_mat

# Set site ID
siteid = 'EC2001_WW3'

#Get date
now = datetime.datetime.now()
nowstr = str(now)
nowstr = nowstr[0:4]+nowstr[5:7]+nowstr[8:10]

#Single file netCDF reading
#ncf='ww3.20130816_hs.nc'
#nco=netCDF4.Dataset(ncf)

#Multiple file netCDF reading
ncf='*_hs.nc'
nco=netCDF4.MFDataset(ncf)

#Get fields to plot
lon=nco.variables['longitude'][:]
lat=nco.variables['latitude'][:]
timeindays=nco.variables['time'][:]
hs=nco.variables['hs'][:]

#Set up a Mercator projection basemap
plt.figure()
m=Basemap(projection='merc',llcrnrlon=lon.min(),urcrnrlon=lon.max(),
llcrnrlat=lat.min(),urcrnrlat=lat.max(),resolution='h')
x,y=m(lon,lat)

#Loop through each time step and plot results
for ind in range(0, len(timeindays)): 
   dt = datetime.datetime.combine(datetime.date(1990, 1, 1), datetime.time(0, 0)) + datetime.timedelta(days=timeindays[ind])
   dstr = datetime.date.strftime(dt,'%Y%m%d%H:%M:%S')
   dstr = dstr[0:8]+' '+dstr[8:17]
   print 'Plotting '+dstr

   #Deal with exception values using a mask
   par=np.double(hs[ind,:])
   par=ma.array(par,mask=np.isnan(par),fill_value=np.nan)

   ramp=np.arange(0,5.0,0.01)

   plt.clf()
   m.contourf(x,y,par,levels=ramp,tri=True,cmap=plt.cm.jet)
   m.colorbar(location='right')
   m.fillcontinents()
   m.drawcoastlines()
   m.drawmapboundary()
   m.drawmeridians(np.arange(-100,-60,5),labels=[0,0,0,5])
   m.drawparallels(np.arange(5,50,5),labels=[5,0,0,0])

   figtitle = siteid+': '+dstr+' (m)'
   plt.title(figtitle)

   dtlabel = datetime.date.strftime(dt,'%Y%m%d%H%M%S')
   dtlabel = dtlabel[0:8]+'_'+dtlabel[8:14]
   filenm = nowstr+'_'+siteid+'_'+dtlabel+'.png'
   plt.savefig(filenm,dpi=150,bbox_inches='tight',pad_inches=0.5)

   del(par)




