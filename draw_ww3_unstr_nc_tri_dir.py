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
siteid = 'Erie_WW3'

#Get date
now = datetime.datetime.now()
nowstr = str(now)
nowstr = nowstr[0:4]+nowstr[5:7]+nowstr[8:10]

#Single file netCDF reading
#ncf='ww3.20130816_hs.nc'
#nco=netCDF4.Dataset(ncf)

#Multiple file netCDF reading
ncf='ww3.????????_hs.nc'
nco=netCDF4.MFDataset(ncf)
ncf2='ww3.????????_dir.nc'
nco2=netCDF4.MFDataset(ncf2)

#Get fields to plot
lon=nco.variables['longitude'][:]
lat=nco.variables['latitude'][:]
timeindays=nco.variables['time'][:]
hs=nco.variables['hs'][:]

timeindays2=nco2.variables['time'][:]
dir=nco2.variables['dir'][:]

reflon=np.linspace(lon.min(),lon.max(),1000)
reflat=np.linspace(lat.min(),lat.max(),1000)
reflon,reflat=np.meshgrid(reflon,reflat)

#Read obs point locations
data=np.loadtxt('erie_ndbc.loc', comments = '$')
buoylon=data[:,0]
buoylat=data[:,1]

plt.figure()

#Loop through each time step and plot results
for ind in range(0, len(timeindays)): 
   dt = datetime.datetime.combine(datetime.date(1990, 1, 1), datetime.time(0, 0)) + datetime.timedelta(days=timeindays[ind])
   dstr = datetime.date.strftime(dt,'%Y%m%d%H:%M:%S')
   dstr = dstr[0:8]+' '+dstr[8:17]
   print 'Plotting '+dstr

   par=np.double(hs[ind,:])
   par2=np.double(dir[ind,:])

   flatness=0.10  # flatness is from 0-.5 .5 is equilateral triangle
   tri=Triangulation(lon,lat)
   mask = TriAnalyzer(tri).get_flat_tri_mask(flatness)
   tri.set_mask(mask)

   tli=LinearTriInterpolator(tri,par)
   par_interp=tli(reflon,reflat)
   tli2=LinearTriInterpolator(tri,par2)
   par2_interp=tli2(reflon,reflat)

   #Set up a Mercator projection basemap
   plt.clf()
   m=Basemap(projection='merc',llcrnrlon=reflon.min(),urcrnrlon=reflon.max(),\
       llcrnrlat=reflat.min(),urcrnrlat=reflat.max(),resolution='h')
   x,y=m(reflon,reflat)

   u=np.cos(np.pi/180*(270-par2_interp))
   v=np.sin(np.pi/180*(270-par2_interp))

   m.pcolormesh(x,y,par_interp,vmax=4.0,shading='flat',cmap=plt.cm.jet)
   m.colorbar(location='right')
   rowskip=np.floor(par2_interp.shape[0]/25)
   colskip=np.floor(par2_interp.shape[1]/25)
   m.quiver(x[0::rowskip,0::colskip],y[0::rowskip,0::colskip],\
         u[0::rowskip,0::colskip],v[0::rowskip,0::colskip], \
         color='black',pivot='middle',units='xy',alpha=0.7)
   #m.drawcoastlines()
   #m.fillcontinents()

   m.drawmeridians(np.arange(int(lon.min()),int(lon.max())+1,1),labels=[0,0,0,1])
   m.drawparallels(np.arange(int(lat.min()),int(lat.max())+1,0.5),labels=[0.5,0,0,0])
   m.drawmapboundary()

   #Plot NDBC locations
   #From Matplotlib Basemap manual, plot(): If latlon keyword is set to True, x,y are 
   # intrepreted as longitude and latitude in degrees. ***Data and longitudes are automatically 
   # shifted to match map projection region for cylindrical and pseudocylindrical projections, 
   # and x,y are transformed to map projection coordinates.***
   m.plot(buoylon,buoylat,'w.',latlon=True,markerfacecolor='w',markeredgecolor='k',markersize=8)

   #m.drawmeridians(np.arange(-100,-60,5),labels=[0,0,0,5])
   #m.drawparallels(np.arange(5,50,5),labels=[5,0,0,0])

   figtitle = siteid+' (12 km): Hsig: '+dstr+' (m)'
   plt.title(figtitle)

   dtlabel = datetime.date.strftime(dt,'%Y%m%d%H%M%S')
   dtlabel = dtlabel[0:8]+'_'+dtlabel[8:14]
   filenm = nowstr+'_'+siteid+'_'+dtlabel+'.png'
   plt.savefig(filenm,dpi=150,bbox_inches='tight',pad_inches=0.5)

   del(par)




