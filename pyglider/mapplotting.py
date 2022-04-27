## copy paste ##
import xarray as xr
import logging
import math
import numpy as np
import matplotlib.pyplot as plt
import os
import yaml
import matplotlib.dates as mdates
import matplotlib.units as munits
import gsw
from datetime import datetime
#test1
import cartopy.crs as ccrs
#
# test2
import matplotlib
import glob
import csv
import simplekml
import seawater
import subprocess
import xarray
import pandas as pd
import pyglider.utils as pygutils
from bs4 import BeautifulSoup

def test(dir,rawnc,figfile):

    converter = mdates.ConciseDateConverter()
    munits.registry[np.datetime64] = converter

    lon0 = -129
    lat0 =  51.5


    Waypoints={}
    Waypoints['QCS02'] = {'lon': -128 - 17/60 - 34 /3600, 'lat': 51 + 41/60 + 39/3600}
    Waypoints['SS5'] = {'lon': -128 - 30/60, 'lat': 51 + 28/60}
    Waypoints['MPA1'] = {'lon': -128 - 43/60 -44 / 3600, 'lat': 51 + 23/60}
    Waypoints['MPA2'] = {'lon': -129 - 55 / 3600, 'lat': 51 + 19/60 + 14 / 3600}
    Waypoints['Shelf'] = {'lon': -129 - 49 /60 - 51 / 3600, 'lat': 51 + 5/60 + 13 / 3600}

    wps = np.zeros((len(Waypoints.keys()), 2))
    for nn, k in enumerate(Waypoints.keys()):
        print('k', k)
        wp = Waypoints[k]
        x, y  = get_xy_lonlat(wp['lon'], wp['lat'], lon0, lat0)
        wps[nn, 0] = x
        wps[nn, 1] = y
    print(wps)



    # In[15]:


    proj = ccrs.Mercator(central_longitude=-129.25, min_latitude=50.5, max_latitude=52)


    fig = plt.figure(figsize=(8, 9), constrained_layout=True)

    gs = fig.add_gridspec(1, 1)
    ax = fig.add_subplot( gs[0, :], projection = proj, facecolor='0.5')

    with xr.open_dataset('/Users/juliaputko/processing/british_columbia_3_msl_2013.nc') as topo0: #jpnote: ? ../../../bathy/british_columbia_3_msl_2013.nc
        topo0 = topo0.isel(lon=slice(6500, 13000, 3), lat=slice(2000, 5500, 3)) 
        #ax = axs['map']
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                            linewidth=0.3, color='gray')
        ax.pcolormesh(topo0.lon, topo0.lat, topo0.Band1, transform=ccrs.PlateCarree(),
                    rasterized=True, alpha=0.4, vmin=-500, vmax=0, zorder=0)
       # ax.set_xlim(-131, -127.5)
        ax.set_extent((-130.5, -128, 50.8, 52))
        
        
   
       # with xr.open_dataset('/Users/juliaputko/processing/dfo-bb046/dfo-bb046-20210511/L0-gridfiles/dfo-bb046-20210511_grid.nc') as ds:  #jpnote: dfo-bb046-20210511_grid.nc ? change to ? was 'L0-gridfiles/dfo-bb046-20210423_grid.nc
        with xr.open_dataset(dir) as ds: 
            ds = get_xy(ds, lat0, lon0)

            ax.plot(ds.longitude, ds.latitude, 'm.', markersize=2,          transform=ccrs.PlateCarree(), zorder=5) 
            ds['alongx'], ind, distline = get_alongx(wps, ds.x, ds.y)

            # fit alongx to get a speed...
            ind = np.where(ds.time > ds.time[-1] - np.timedelta64(36,'h') )[0][0]
            dt = np.float(ds.time[-1] - ds.time[ind]) / 1e9 / 24/ 3600
            dx = ds.alongx[-1] - ds.alongx[ind]
            tt = (ds.time[ind:] - ds.time[ind]).astype(float) / 24/3600 /1e9
            print(tt)
            p = np.polyfit(tt, ds.alongx[ind:].values, 1)
            print('p', p)

        for nn, k in enumerate(Waypoints.keys()):
            wp = Waypoints[k]
            x, y  = get_xy_lonlat(wp['lon'], wp['lat'], lon0, lat0)
            wps[nn, 0] = x
            wps[nn, 1] = y
            print(wp['lon'], wp['lat'], x, y)
            ax.plot(wp['lon'], wp['lat'],  'o', color='0.75', zorder=1, markersize=10, transform=ccrs.PlateCarree(), 
                alpha=0.5) 
            ax.text(wp['lon'], wp['lat'], f'    {k} {distline[nn]:1.1f} km',
                    transform=ccrs.PlateCarree(), fontsize=8, fontstyle='italic', color='1', alpha=0.5, #jpnote maybe this?
                    verticalalignment='center')
          #  print(wp.lon, wp.lat, 'd')
       # ax = fig.add_subplot(gs[1, 0])
       # ax.plot(ds.time[ind:], np.polyval(p, tt), label=f'{p[0]:1.1f} km/d (last 36h)', color='0.5', ls='--')
       # ax.plot(ds.time, ds.alongx)
        #ax.legend()
        #ax.set_ylabel('Along-line [km]')

  #  ax = fig.add_subplot(gs[1, 1])
    #with xr.open_dataset('/Users/juliaputko/processing/dfo-bb046/dfo-bb046-20210511/rawnc_realtime/dfo-bb046-rawgli.nc') as gli:  #jpnote ? 
    with xr.open_dataset(rawnc) as gli:
        gli=gli.isel(time=slice(100,None))
       # ax.plot(gli.time, gli.Voltage, label=f'min: {np.min(gli.Voltage.values):1.1f} V')
        #ax.set_ylabel('Voltage [V]')
        #ax.set_ylim(bottom=22.5)
        #ax.axhline(22.8, color='r', linestyle='--', linewidth=2)
        #ax.legend()
        print(gli.Lon)

   # now = np.datetime64(datetime.datetime.now()) + np.timedelta64(7, 'h')
    now = np.datetime64(datetime.now()) + np.timedelta64(7, 'h')
    st = f'{now}'
    st = st[:-10]
    st2 = f'{gli.time[-1].values}'
    st2 = st2[:-10]
    fig.suptitle(f'Processed at {st}, data: {st2}', fontsize=8)
    fig.canvas.draw()
    fig.savefig(figfile, dpi=200)

## end of map test 1 


## map test 2


def test2():
    fig = plt.figure(figsize=(14,7), constrained_layout=False)
    #fig = plt.figure(figsize=(12,9), constrained_layout=False)
    #gs = fig.add_gridspec(3, 2, height_ratios=[0.7, 0.3, 1.1], bottom=0.02, top=0.98, left=0.06, right=0.98, hspace=0.02, wspace=0.02)
    gs = fig.add_gridspec(3, 2, height_ratios=[9, 0.5, 1.1], bottom=0.02, top=0.90, left=0.01, right=0.98, hspace=0.02, wspace=0.02)
    #gs = fig.add_gridspec(3,2,height_ratios=[0.7, 0.3, 1.1],bottom=0.02, top=0.98, left=0.06, right=0.98, hspace=0.02, wspace=0.02)
    ax = fig.add_subplot(gs[0, :])
    #ax = fig.subplots(1)
    
    if 0:
        with xr.open_dataset('/Users/cproof/processing/deployments/dfo-eva035/dfo-eva035-20190718/rawnc_realtime/dfo_eva035-rawgli.nc') as ds:  #jp changed path '../../dfo_eva035/dfo_eva035-20190718/rawnc_realtime/dfo_eva035-rawgli.nc'
            ind = np.unique(ds.Lon, return_index=True)[1]
            ind = np.sort(ind)

            lon = pygutils.nmea2deg(ds.Lon[ind])
            lat = pygutils.nmea2deg(ds.Lat[ind])
           # ax.plot(lon, lat, '.', markersize=5, alpha=0.8, label='', mec='none', color='C2')
           # ax.plot(lon, lat, label='', color='C2', alpha=0.5, linewidth=0.7)
            for ind in range(1, 2):
                label = ''
                if ind == 1:
                    timest = f'{ds.time.values[-1]}'[:-13]
                    label = f'dfo-eva035 ' + timest
                size = (10 - ind)/10*8 + 1.3
                size=8
                print(size)
             #   p, = ax.plot(lon[-ind], lat[-ind], '.',
             #           markersize=size, alpha=1, color='C2', label=label, mec='none')
             #   ax.plot(lon[-1], lat[-1], 'o', color='C2', alpha=0.5,
             #           markersize=10, label='', markeredgecolor='none', zorder=-10)
                


    xmln = '/Users/jklymak/gliderdata/slocum_dockserver/wall_e_652/gliderState.xml'
    from bs4 import BeautifulSoup
    with open(xmln, 'r') as fin:
        y=BeautifulSoup(fin, features="xml")
        time = None
        for a in y.find_all('valid_location'):
            try:
                dtime = np.datetime64(a.time.text)
                if dtime > np.datetime64('2019-07-19'):
                    #print(dtime, float(a.lat.text), float(a.lon.text))
                    #print(np.array(dtime))
                    if time is not None:
                        time = np.append(time, dtime)
                        lat = np.append(lat, float(a.lat.text))
                        lon = np.append(lon, float(a.lon.text))
                    else:
                        time = np.array(dtime)
                        lat = np.array(float(a.lat.text))
                        lon = np.array(float(a.lon.text))
            except:
                pass
        lon = pygutils.nmea2deg(lon)
        lat = pygutils.nmea2deg(lat)
        timest = f'{time[-1]}'[:-3]
        label = f'dfo-walle652 ' + timest

        ax.plot(lon, lat, '.', color='C1', markersize=10, alpha=0.2, label='', mec='none') #markersize 4
        ax.plot(lon, lat, 'C1', alpha=0.5, linewidth=3, label='')
        ax.plot(lon[-1], lat[-1], '.', color='r', alpha=1, markersize=12, label=label, mec='none')
        ax.plot(lon[-1], lat[-1], 'o', color='C1', alpha=0.5, markersize=10,
                    label='', markeredgecolor='none', zorder=-10)
        wallelon = lon[-1]
        wallelat = lat[-1]
        # make a KML?
        if 0:
            kml = simplekml.Kml()
            coords = []
            for thelon, thelat  in zip(lon, lat):
                coords += [(thelon, thelat)]
            pnt = kml.newpoint(coords=[(thelon, thelat)])
            pnt.style.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/placemark_square_highlight.png'
            pnt.timestamp.when = timest
            kml.newlinestring(coords=coords, name="dfo-walle652")
            kml.save("dfo-walle652.kml")
        # get the waypoint:
        nexty = False
        for ch in y.dataParameters.children:
            for cch in ch.children:
                print(cch.text)
                if nexty:
                    nextwpt = cch.text
                if cch.text == 'Next waypoint coordinates':
                    nexty = True
                else:
                    nexty = False
        wlat, wlon = nextwpt[1:-1].split(',')
        wlat = pygutils.nmea2deg(float(wlat))
        wlon = pygutils.nmea2deg(float(wlon))
        print('LL', wlon, wlat)
       # ax.plot(wlon, wlat, 'd', color='C1', markersize=15, alpha=0.5, mec='g', zorder=-20)

    wallelons = lon
    wallelats = lat
    walletime = time

    #    ax.set_title(ds.time.values[-1], )
    evawps = np.array([[49+08.75/60, -(130+40.6/60)],
    [48+58.50/60, -(131+05.15/60)],
    [48+52.85/60, -(130+43.15/60)],
    [49+09.65/60, -(131+05.35/60)]])

    wallewps =np.array([[49+12.90/60, -(130+54.45/60)],
    [48+52.45/60, -(130+57.05/60)],
    [49+01.35/60, -(130+38.40/60)],
    [49+03.45/60, -(131+13.50/60)]])

    # line P
    inst = """#P1 -12530.0000   4834.5000
    #P2 -12600.0000   4836.0000
    #P3 -12620.0000   4837.5000
    #P4 -12640.0000   4839.0000
    #P5 -12710.0000   4841.5000
    #P6 -12740.0000   4844.6000
    #P7 -12810.0000   4846.6000
    #P8 -12840.0000   4849.0000
    #P9 -12910.0000   4851.4000
    #P10 -12940.0000   4853.6000
    #P11 -13010.0000   4856.0000
    #P12 -13040.0000   4858.2000
    #P13 -13140.0000   4902.6000
    #P14 -13240.0000   4907.4000
    #P15 -13340.0000   4912.0000
    #P16 -13440.0000   4917.0000
    #P17 -13540.0000   4921.0000
    #P18 -13640.0000   4926.0000
    #P19 -13740.0000   4930.0000
    #P20 -13840.0000   4934.0000
    #P21 -13940.0000   4938.0000
    #P22 -14040.0000   4942.0000
    #P23 -14140.0000   4946.0000
    #P24 -14240.0000   4950.2000
    #P25 -14336.3000   5000.0000
    #P35 -14418.2000   5000.0000
    #P26 -14500.0000   5000.0000
    """
    import io
    f = io.StringIO(inst)
    wallewps = None
    wallenm  = []
    for line in f:
        st = line.split(' ')
        print(st)
        print(st[0], st[1], st[-1])
        if st[1] == '':  #jp '' ?? 
            st[1]== 0
            lon = int(0 / 100)
            dlon = (0 - lon * 100)
            lon = lon + dlon / 60
            lat = int(0/ 100)
            dlat = (0 - lat * 100)
        else: 
            lon = int(float(st[1]) / 100)
            dlon = (float(st[1]) - lon * 100)
            lon = lon + dlon / 60
            lat = int(float(st[-1]) / 100)
            dlat = (float(st[-1]) - lat * 100)
        lat = lat + dlat / 60
        print(lat)
        if wallewps is None:
            wallewps = np.array([[lat, lon]])
        else:
            wallewps = np.append(wallewps, np.array([[lat, lon]]), axis=0)
        wallenm += [st[0][1:]]
    print('WALL', wallenm, wallewps)
    # get distances
    dist = np.zeros(len(wallewps))
    dist, ang  = seawater.extras.dist(wallewps[:, 0], wallewps[:, 1])
    dist = np.cumsum(np.append([0], dist))
    walledist = dist - dist[3]
    walledist = -walledist  # make Papa at -1300 km
    print('dist', dist)


    # now figure out a speed along the line....
    # figure out a speed
    dist = np.interp(wallelons, wallewps[:, 1][::-1], walledist[::-1])
    gg = (walletime >= walletime[-1] - np.timedelta64(5,'D'))

    dt = (walletime[gg][-1] - walletime[gg][0])
    print(dt)
    speed = -(dist[gg][0] - dist[gg][-1]) / dt.astype(float) * 24 * 3600
    print(walletime[gg][0])
    print('speed', speed,' km/d')

    gg = (walletime <= walletime[-1] - np.timedelta64(1,'h'))
    dt = (walletime[-1] - walletime[gg][-1])
    lastspeed = (dist[-1] - dist[gg][-1]) / dt.astype(float) * 24 * 3600


    p4dist = -dist[gg][-1]
    p26dist = (walledist[-1] - dist[gg][-1])

    print(time[-1].tolist())
    lastconnect = datetime.utcnow() - time[-1].tolist() #jpnote removed 1 datetime
    print(lastconnect)

    # get time offset...
    timeoff = datetime.utcnow() - datetime.now() #jpnote removed 1 datetime
    print(timeoff)
    utcnow = str(datetime.utcnow())[:-10] #jpnote removed 1 datetime
    now = str(datetime.now())[-15:-10] #jpnote removed 1 datetime


    infost = 'Plot update:    ' + utcnow + ' ('+ now + ')\n'
    surfacing = time[-1]
    surfst = str(surfacing)[:-3].replace('T', ' ')
    localst = str(surfacing.tolist() - timeoff)[11:16]
    infost += 'Last surfacing: ' + surfst + ' (' + localst + ')\n'
    infost += ('Last connect:                     ' + f'{str(lastconnect).split(":")[0]:>2s}' + ':'
            + str(lastconnect).split(':')[1] + '\n')
    infost += f'Dist to P4:  {p4dist:5.0f} km ~{p4dist/speed:4.0f} d\n'
    infost += f'Dist to P26: {p26dist:5.0f} km ~{p26dist/speed:4.0f} d\n'
    infost += f'Speed (last five days): {speed:1.1f} km/d\n'
    infost += f'Speed (last section):   {lastspeed:1.1f} km/d'

#original textbox
   # ax.text(0.03, 0.1, infost, transform=ax.transAxes, fontsize=18,
   #         bbox={'facecolor':'w', 'alpha':0.5}, family='monospace')
#lower textbox
   # ax.text(0.03, -0.53, infost, transform=ax.transAxes, fontsize=18,
   #         bbox={'facecolor':'w', 'alpha':0.5}, family='monospace')
   
    #ax.plot(evawps[:, 1], evawps[:, 0], 'd', color='C2', alpha=0.75)
   # ax.plot(wallewps[:, 1], wallewps[:, 0], 'd', color='C1', alpha=0.75)  #jpnote comment one of these
    #for i in range(len(wallewps[:, 0])):
    #    ax.text(wallewps[i, 1], wallewps[i, 0], ' ' + wallenm[i] + f'; {walledist[i]:1.0f} km', color='C1', fontsize=8, zorder=-2, alpha=0.75, rotation=30)  #jp commented to remove text and points along line 
    #for i in range(len(evawps[:, 0])):
      #  ax.text(evawps[i, 1], evawps[i, 0], f' E{(i+1):d}', color='C2', fontsize=8, zorder=-2, alpha=0.75)
    ax.set_aspect(1/np.cos(np.deg2rad(49)))

    ax.legend(loc='lower right')


    with open('/Users/cproof/processing/deployments/Explorer19.stations.csv') as fin: ##jpnote find path /Users/cproof/processing/deployments  ../../Explorer19.stations.csv'
        reader = csv.DictReader(fin)
        for row in reader:
            print(row)
           # ax.plot(float(row['Long']), float(row['Lat']), 'd', markersize=4, color='0.6', alpha=0.3, zorder=-20)
          #  ax.text(float(row['Long']), float(row['Lat']), ' ' + row['StaName'],
               #     fontsize=6, color='0.6', alpha=0.7, clip_on=True)
    with open('/Users/cproof/processing/deployments/Explorer19.lines1.csv') as fin: # '../../Explorer19.lines1.csv'
        reader = csv.DictReader(fin)
        for row in reader:
            print(row)

            x = [float(row['START_X']), float(row['END_X'])]
            y = [float(row['START_Y']), float(row['END_Y'])]

           # ax.plot(x, y, '-', color='0.2', alpha=0.3, zorder=-20)
    axmap = ax
    ax.set_title('Last update: ' + str(datetime.utcnow())[:-10])

    ## add smith sandwell

    with xr.open_dataset('/Users/jklymak/gliderdata/smith_sandwell_topo_v8_2.nc') as ds:
        # longitude, latitiude, ROSE
        ind = np.where((ds.longitude < -124 + 360) & (ds.longitude > -145.2 + 360))[0]
        ds = ds.isel(longitude=ind)
        ind = np.where((ds.latitude < 52.2) & (ds.latitude > 46.0))[0]
        ds = ds.isel(latitude=ind)
        ds['longitude'] = ds.longitude - 360

        ax.contourf(ds.longitude, ds.latitude, -ds.ROSE, linewidths=0, cmap='Blues',
                levels=np.arange(0, 6000, 100), vmax=1000, vmin=-1000, alpha=0.9, zorder=-100)
    ax.set_ylim([46.1, 52.2])
    ax.set_xlim([-145.2, -124.4])
   
    # Now get some info from the logs....


    #ax = fig.add_subplot(gs[1, 0]) #jp commented
    ## make a time series of amphours....

#    full = glob.glob('/Users/juliaputko/processing/dfo-walle652/dfo-walle652-20210121/realtime_rawnc/dfo-walle652*-rawdbd.nc') #jpnote change /Users/juliaputko/processing/dfo-walle652/dfo-walle652-20210121  'realtime_rawnc/dfo-walle652*-rawdbd.nc'
#    full.sort()
#    with xr.open_dataset(full[-1]) as ds:
#        d = ds.where(np.isfinite(ds['m_coulomb_amphr_total']), drop=True)
#        print(d)
#        ax.plot(d.time, d.m_coulomb_amphr_total,'.')
#        drec =  d.where(ds.time > ds.time[-1] - np.timedelta64(5, 'D'), drop=True)
#        print(drec)
#        dt = (drec.time[-1] - drec.time[0]).astype(float) / 1e9 /24 / 3600
#        print(dt)
#        dah = drec.m_coulomb_amphr_total[-1] - drec.m_coulomb_amphr_total[0]
#        print(dah/dt)
   # ax.set_ylabel('[amp-hours]')
  
#    txt = str(drec.time.values[-1])[:16].replace('T', ' ') + '\n'
#    txt += 'Current consumption: %1.0f amph\n' % (drec.m_coulomb_amphr_total.values[-1])
#    txt += 'Rate (five days):    %1.2f amph/d\n' % (dah/dt)
#    daysp4 = p4dist/speed
#    txt += 'Expected by P4:      %1.0f amph' % (drec.m_coulomb_amphr_total.values[-1] + daysp4 * dah/dt)
   # ax.text(0.02, 0.1, txt,
   #         transform=ax.transAxes, family='monospace', fontsize=14,
   #         bbox={'facecolor':'w', 'alpha':0.5})
 


  #  ax = fig.add_subplot(gs[2, 0])  #jp commented
#    logdir = '/Users/jklymak/gliderdata/slocum_dockserver/wall_e_652/logs/'
#    logs = glob.glob(logdir + '/*')
#    logs.sort(key=lambda x: os.path.getmtime(x))
#    logs = logs[::-1]



    #ax.text(0.02, 0.9, str, transform=ax.transAxes, family='monospace', fontsize=14)

    # not get digifins etc...
 #   out = ''
 #   done = False
 #   for log in logs:
 #       with open(log, 'r') as fin:
 #           found = False
 #           for line in fin:
 #               if 'name       lim' in line:
 #                   found = True
 #                   out = log + '\n'
 #               if found:
 #                   out += line
 #                   if 'devices:(t/m/s) errs:' in line:
 #                       break
 #       if found:
 #           break

    #ax.text(0.02, 0.02, out, transform=ax.transAxes, family='monospace', fontsize=10)
    #ax.axis('off')


    #ax = fig.add_subplot(gs[1, 1]) #jp commented
 #   out0 = subprocess.run(['ls -hlt ../../../slocum_dockserver/wall_e_652/from-glider | head -n 10'], shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8') #jpnote commented
 #   lines = out0.split('\n')
 #   out = ''
 #   for l in lines[1:]:
 #       out += l[35:] + '\n'
 #   print(out)
 #   out = out[:-2]
 #   ax.text(0.02, 0.0, out , transform=ax.transAxes, family='monospace', fontsize=12)
    #ax.axis('off')
    axmap.set_ylim(np.array([-1.5, 1.4]) + wallelat)
    axmap.set_xlim(np.array([-8.5, 1.7]) + wallelon)

   
    fig.savefig('/Users/juliaputko/processing/dfo-walle652/dfo-walle652-20210121/figs/Explorer19MissionMapBig.png', dpi=200) #jpnote change './figs/Explorer19MissionMapBig.png'

    axmap.set_ylim(np.array([-0.75, 0.25]) + wallelat)
    axmap.set_xlim(np.array([-1.7, 1.7]) + wallelon)

    fig.savefig('/Users/juliaputko/processing/dfo-walle652/dfo-walle652-20210121/figs/Explorer19MissionMap.png', dpi=200) #jpnote change './figs/Explorer19MissionMap.png'

