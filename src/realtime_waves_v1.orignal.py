def read_data(filename,latmax=None) :
### read in u,v,z data from file could be replaced by other routine
    import numpy as np
    import netCDF4 as nc

    if latmax == None :
        latmax=24

    fields=nc.Dataset('uvz_6hr_moan83day_mofc7day_20190101.nc','r')
    latrange=np.where(abs(fields.variables['latitude'][:]) <= latmax)

    u=fields.variables['u'][:,:,latrange[0],:]
    v=fields.variables['v'][:,:,latrange[0],:]
    z=fields.variables['ht'][:,:,latrange[0],:]
    lons=fields.variables['longitude'][:]
    lats=fields.variables['latitude'][latrange[0]]
    press=fields.variables['p'][:]
    times=fields.variables['t'][:]

    return u,v,z, lons,lats,press,times

def uz_to_qr(u,z,g_on_c) :
# transform u,z to q, r using q=z*(g/c) + u; r=z*(g/c) - u 

    q=z*g_on_c+u
    r=z*g_on_c-u

    return q,r

def filt_project (qf,rf,vf,lats,y0,waves,pmin,pmax,kmin,kmax,c_on_g)  :
    
    import numpy as np

# find size of arrays 
    nf=qf.shape[0]  
    nz=qf.shape[1]
    nlats=lats.size
    nk=qf.shape[3]
# Find frequencies and wavenumbers corresponding to pmin,pmax and kmin,kmax in coeff matrices
    f=np.fft.fftfreq(nf,0.25)
    fmin=np.where((f >= 1./pmax))[0][0]
    fmax=(np.where((f > 1./pmin))[0][0])-1
    f1p=fmin
    f2p=fmax+1
    f1n=nf-fmax 
    f2n=nf-fmin+1
    k1p=kmin
    k2p=kmax+1
    k1n=nk-kmax
    k2n=nk-kmin+1
 
### note that need to adjust for the fact that array referencing doesn't include last point!!!!!!

# Define the parobolic cylinder functions
    spi2=np.sqrt(2*np.pi)
    dsq=np.array([spi2,spi2,2*spi2,6*spi2]) # Normalization for the 1st 4 paroblic CF
    d=np.zeros([dsq.size,nlats])
    y=lats[:]/y0
    ysq=y**2
    d[0,:]=np.exp(-ysq/4.0)
    d[1,:]=y*d[0,:]
    d[2,:]=(ysq-1.0)*d[0,:]
    d[3,:]=y*(ysq-3.0)*d[0,:]
    dlat=np.abs(lats[0]-lats[1])
    
    qf_Kel=np.zeros([nf,nz,nk],dtype=complex)
    qf_mode=np.zeros([dsq.size,nf,nz,nk],dtype=complex)
    rf_mode=np.zeros([dsq.size,nf,nz,nk],dtype=complex)
    vf_mode=np.zeros([dsq.size,nf,nz,nk],dtype=complex)

# reorder the spectral coefficents to make the latitudes the last dimension
    qf=np.transpose(qf,[0,1,3,2])
    rf=np.transpose(rf,[0,1,3,2])
    vf=np.transpose(vf,[0,1,3,2])

    for m in np.arange(dsq.size) :
        if m == 0:   # For eastward moving Kelvin Waves
            qf_Kel[f1n:f2n,:,k1p:k2p]=np.sum(qf[f1n:f2n,:,k1p:k2p,:]\
                   *np.squeeze(d[m,:])*dlat,axis=-1)/(dsq[m]*y0)
            qf_Kel[f1p:f2p,:,k1n:k2n]=np.sum(qf[f1p:f2p,:,k1n:k2n,:]\
                   *np.squeeze(d[m,:])*dlat,axis=-1)/(dsq[m]*y0)
        # For westward moving waves
        qf_mode[m,f1n:f2n,:,k1n:k2n]=np.sum(qf[f1n:f2n,:,k1n:k2n,:]\
                   *np.squeeze(d[m,:])*dlat,axis=-1)/(dsq[m]*y0)
        qf_mode[m,f1p:f2p,:,k1p:k2p]=np.sum(qf[f1p:f2p,:,k1p:k2p,:]\
                   *np.squeeze(d[m,:])*dlat,axis=-1)/(dsq[m]*y0)
        rf_mode[m,f1n:f2n,:,k1n:k2n]=np.sum(rf[f1n:f2n,:,k1n:k2n,:]\
                   *np.squeeze(d[m,:])*dlat,axis=-1)/(dsq[m]*y0)
        rf_mode[m,f1p:f2p,:,k1p:k2p]=np.sum(rf[f1p:f2p,:,k1p:k2p,:]\
                   *np.squeeze(d[m,:])*dlat,axis=-1)/(dsq[m]*y0)
        vf_mode[m,f1n:f2n,:,k1n:k2n]=np.sum(vf[f1n:f2n,:,k1n:k2n,:]\
                   *np.squeeze(d[m,:])*dlat,axis=-1)/(dsq[m]*y0)
        vf_mode[m,f1p:f2p,:,k1p:k2p]=np.sum(vf[f1p:f2p,:,k1p:k2p,:]\
                   *np.squeeze(d[m,:])*dlat,axis=-1)/(dsq[m]*y0)
    
    uf_wave=np.zeros([waves.size,nf,nz,nlats,nk],dtype=complex)
    zf_wave=np.zeros([waves.size,nf,nz,nlats,nk],dtype=complex)
    vf_wave=np.zeros([waves.size,nf,nz,nlats,nk],dtype=complex)

    for w in np.arange(waves.size) :
        if waves[w] == 'Kelvin' :
            for j in np.arange(nlats) :
                uf_wave[w,:,:,j,:]=\
                    0.5*qf_Kel[:,:,:]*d[0,j]
                zf_wave[w,:,:,j,:]=\
                    0.5*qf_Kel[:,:,:]*d[0,j]*c_on_g

        if waves[w] == 'WMRG' :
            for j in np.arange(nlats) :
                uf_wave[w,:,:,j,:]=\
                    0.5*qf_mode[1,:,:,:]*d[1,j]
                zf_wave[w,:,:,j,:]=\
                    0.5*qf_mode[1,:,:,:]*d[1,j]*c_on_g
                vf_wave[w,:,:,j,:]=\
                    vf_mode[0,:,:,:]*d[0,j]

        if waves[w] == 'R1' :
            for j in np.arange(nlats) :
                uf_wave[w,:,:,j,:]=\
                    0.5*(qf_mode[2,:,:,:]*d[2,j]-rf_mode[0,:,:,:]*d[0,j])
                zf_wave[w,:,:,j,:]=\
                    0.5*(qf_mode[2,:,:,:]*d[2,j]+rf_mode[0,:,:,:]*d[0,j])*c_on_g
                vf_wave[w,:,:,j,:]=\
                    vf_mode[1,:,:,:]*d[1,j]
        if waves[w] == 'R2' :
            for j in np.arange(nlats) :
                uf_wave[w,:,:,j,:]=\
                    0.5*(qf_mode[3,:,:,:]*d[3,j]-rf_mode[1,:,:,:]*d[1,j])
                zf_wave[w,:,:,j,:]=\
                    0.5*(qf_mode[3,:,:,:]*d[3,j]+rf_mode[1,:,:,:]*d[1,j])*c_on_g
                vf_wave[w,:,:,j,:]=\
                    vf_mode[2,:,:,:]*d[2,j]

    return uf_wave,zf_wave,vf_wave

def write_data(u_wave,z_wave,v_wave,lons,lats,press,times,waves) :

    import numpy as np
    import netCDF4 as nc

    outfile='./waves_6hr_moan83day_mofc7day_20190101.nc'
    ncout=nc.Dataset(outfile,'w')
    ncout.description='Wave field anomalies using realtime methodoloy of Yang (2019)'
    ncout.history='Created by realtime waves.py'

    ncout.createDimension('longitude',lons.size)
    varunit=ncout.createVariable('longitude',lons.dtype,('longitude'))
    varunit.setncattr('units','degrees_east')
    varunit.setncattr('long_name','longitude')
    ncout.variables['longitude'][:]=lons

    ncout.createDimension('latitude',lats.size)
    varunit=ncout.createVariable('latitude',lats.dtype,('latitude'))
    varunit.setncattr('units','degrees_east')
    varunit.setncattr('long_name','latitude')
    ncout.variables['latitude'][:]=lats

    ncout.createDimension('pressure',press.size)
    varunit=ncout.createVariable('pressure',press.dtype,('pressure'))
    varunit.setncattr('units','degrees_east',)
    varunit.setncattr('long_name','pressure')
    ncout.variables['pressure'][:]=press

    ncout.createDimension('time',times.size)
    varunit=ncout.createVariable('time',times.dtype,('time'))
    varunit.setncattr('units','days since 2019-01-01 00:00:00')
    varunit.setncattr('long_name','time')
    ncout.variables['time'][:]=times

    ncout.createDimension('wave_type',waves.size)
    varunit=ncout.createVariable('wave_type',waves.dtype,('wave_type'))
    varunit.setncattr('long_name','wave_type')
    ncout.variables['wave_type'][:]=waves

    u_unit=ncout.createVariable('u_wave',np.float32,('wave_type','time','pressure','latitude','longitude'))
    u_unit.setncattr('units','m s-1')
    u_unit.setncattr('longname','Wave Zonal Wind')
    u_unit.setncattr('standard_name','eastward_wind')

    v_unit=ncout.createVariable('v_wave',np.float32,('wave_type','time',\
           'pressure','latitude','longitude'))
    v_unit.setncattr('units','m s-1')
    v_unit.setncattr('longname','Wave Meridional Wind')
    v_unit.setncattr('standard_name','northward_wind')
    
    z_unit=ncout.createVariable('z_wave',np.float32,('wave_type','time',\
           'pressure','latitude','longitude'))
    z_unit.setncattr('units','m')
    z_unit.setncattr('longname','Wave Geopotential Height')
    z_unit.setncattr('standard_name','geopotential_height')

    ncout.variables['u_wave'][:]=u_wave
    ncout.variables['v_wave'][:]=v_wave
    ncout.variables['z_wave'][:]=z_wave

    return

    
#---------------------------------------------------------
# start of main routine
    
import numpy as np

# define some phyiscal parameters
g=9.8
beta=2.3e-11
radea=6.371e6
spd=86400.
ww=2*np.pi/spd


# Define some parameters spefic to the methodology
latmax=24.  #   +/- latitude range over which to process data.
kmin=2      # minimum zonal wavenumber
kmax=40     # maximum zonal wavenumber
pmin=2.0      # minimum period (in days)
pmax=30.0   # maximum period (in days)
y0=6.0      # meridional trapping scale (degrees)
waves=np.array(['Kelvin','WMRG','R1','R2']) # List of wave types to output
#waves=np.array(['Kelvin','WMRG']) # List of wave types to output

y0real= 2*np.pi*radea*y0/360.0   # convert trapping scale to metres
ce=2*y0real**2*beta
g_on_c=g/ce
c_on_g=ce/g


### read in 90 days of u,v,z data at 6 hourly time resolution and equatorial grid point
###  Methodology test and applied with 1 degree resolution data, but not dependent on it
### For Met Office operational forecasts the data would be 83 days of analysis and 7 days of forecast

#read data
u,v,z, lons,lats,press,times = read_data('uvz_6hr_moan83day_mofc7day_20190101.nc',latmax=24)


#convert u,z to q,r
q,r = uz_to_qr(u,z,g_on_c)


# Fourier transform in time and longitude

qf=np.fft.fft2 (q,axes=(0,-1))
rf=np.fft.fft2 (r,axes=(0,-1))
vf=np.fft.fft2 (v,axes=(0,-1))


# Project onto individual wave modes
uf_wave,zf_wave,vf_wave=filt_project(qf,rf,vf,lats,y0,waves,pmin,pmax,kmin,kmax,c_on_g)

# Inverse Fourier transform in time and longitude

u_wave=np.real(np.fft.ifft2 (uf_wave,axes=(1,-1)))
z_wave=np.real(np.fft.ifft2 (zf_wave,axes=(1,-1)))
v_wave=np.real(np.fft.ifft2 (vf_wave,axes=(1,-1)))

print(u_wave[:,0,0,0,0])
print(z_wave[:,0,0,0,0])
print(v_wave[:,0,0,0,0])


write_data(u_wave,z_wave,v_wave,lons,lats,press,times,waves)


# End of main routine
###-----------------------------------------------------------------
