# Script for plotting unstructured WW3 output files in netCDF format, created by ww3_ounf
# Andre van der Westhuysen, 08/21/13
#
import matplotlib
matplotlib.use('Agg',warn=False)

import datetime
from matplotlib.tri import Triangulation, TriAnalyzer, LinearTriInterpolator
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from mpl_toolkits.basemap import Basemap

# Set site ID
siteid = 'Shinnecock_WW3'

#Get date
now = datetime.datetime.now()
nowstr = str(now)
nowstr = nowstr[0:4]+nowstr[5:7]+nowstr[8:10]

#Single file netCDF reading
#ncf='ww3.20130816_hs.nc'
#nco=netCDF4.Dataset(ncf)

#Multiple file netCDF reading
ncf='ww3.Constant.20151214_hs.nc'
nco=netCDF4.Dataset(ncf)

#Get fields to plot
lon=nco.variables['longitude'][:]
lat=nco.variables['latitude'][:]
timeindays=nco.variables['time'][:]
hs=nco.variables['hs'][:]

reflon=np.linspace(lon.min(),lon.max(),1000)
reflat=np.linspace(lat.min(),lat.max(),1000)
reflon,reflat=np.meshgrid(reflon,reflat)

plt.figure()

#Loop through each time step and plot results
for ind in range(0, len(timeindays)): 
   dt = datetime.datetime.combine(datetime.date(1990, 1, 1), datetime.time(0, 0)) + datetime.timedelta(days=timeindays[ind])
   dstr = datetime.date.strftime(dt,'%Y%m%d%H:%M:%S')
   dstr = dstr[0:8]+' '+dstr[8:17]
   print 'Plotting '+dstr

   par=np.double(hs[ind,:])

   flatness=0.2  # flatness is from 0-.5 .5 is equilateral triangle
   tri=Triangulation(lon,lat)
   mask = TriAnalyzer(tri).get_flat_tri_mask(flatness)
   tri.set_mask(mask)

   tli=LinearTriInterpolator(tri,par)
   par_interp=tli(reflon,reflat)

   #Set up a Mercator projection basemap
   plt.clf()
   m=Basemap(projection='merc',llcrnrlon=reflon.min(),urcrnrlon=reflon.max(),\
       llcrnrlat=reflat.min(),urcrnrlat=reflat.max(),resolution='h')

   x,y=m(reflon,reflat)
   m.pcolormesh(x,y,par_interp,vmax=1.4,shading='flat',cmap=plt.cm.jet)
   m.colorbar(location='right')
   m.drawcoastlines()
   m.fillcontinents()

   m.drawmeridians(np.arange(lon.min(),lon.max(),0.2),labels=[0,0,0,0.2])
   m.drawparallels(np.arange(lat.min(),lat.max(),0.2),labels=[0.2,0,0,0])
   m.drawmapboundary()

   #m.drawmeridians(np.arange(-100,-60,5),labels=[0,0,0,5])
   #m.drawparallels(np.arange(5,50,5),labels=[5,0,0,0])

   figtitle = siteid+': '+dstr+' (m)'
   plt.title(figtitle)

   dtlabel = datetime.date.strftime(dt,'%Y%m%d%H%M%S')
   dtlabel = dtlabel[0:8]+'_'+dtlabel[8:14]
   filenm = nowstr+'_'+siteid+'_'+dtlabel+'_tri.png'
   plt.savefig(filenm,dpi=150,bbox_inches='tight',pad_inches=0.5)

   del(par)




