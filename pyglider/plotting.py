import xarray as xr
import logging
import numpy as np
import matplotlib.pyplot as plt
import os
import yaml
import matplotlib.dates as mdates
import matplotlib.units as munits
import seawater

_log = logging.getLogger(__name__)

try:
    converter = mdates.ConciseDateConverter()
    munits.registry[np.datetime64] = converter
except:
    # older matplotlib...
    pass

plt.rcParams['figure.constrained_layout.use']=True

def _autoclim(vals):
    vals = vals.values.flatten()
    vals = vals[~np.isnan(vals)]
    min = np.min(vals)
    max = np.max(vals)

    m15 = np.quantile(vals, 0.15 / 2)
    m85 = np.quantile(vals, 1 - 0.15 / 2)
    d = m85 - m15
    m15 = m15 - d / 10
    m85 = m85 + d / 10
    print('m15, m85', m15, m85, np.max((min, m15)))
    return np.max((min, m15)), np.min((max, m85))

def timeseries_plots(fname, plottingyaml):

    with open(plottingyaml) as fin:
        config = yaml.safe_load(fin)

    try:
        os.mkdir(config['figdir'])
    except:
        pass

    with xr.open_dataset(fname, decode_times=True) as ds:
        # map!
        fig, axs = plt.subplots(3, 1, gridspec_kw={'height_ratios': [3, 1, 1]})
        ax = axs[0]
        ax.plot(ds.longitude, ds.latitude, '.')
        ax.set_aspect(1 / np.cos(np.deg2rad(ds.latitude.mean())))
        ax.set_ylabel('Lat [degrees north]')
        ax.set_xlabel('Lon [degrees east]')

        ax = axs[1]
        ax.plot(ds.time, ds.longitude, '.')
        ax.set_ylabel('Lon [degrees east]')

        ax = axs[2]
        ax.plot(ds.time, ds.latitude, '.')
        ax.set_ylabel('Lat [degrees north]')
        fig.savefig(config['figdir'] + '/map%s.png'%ds.attrs['deployment_name'], dpi=200)

        # timeseries of things....
        _log.info('Plotting timeseries data')
        keys = config['timeseries'].keys()
        N = len(keys)
        if 1:
            fig, axs = plt.subplots(int(N / 2), 2, figsize=(7.5, 7),
                                    sharex=True, sharey=False)
            axs = axs.flat
            for n, k in enumerate(keys):
                print('key', k)
                if config['timeseries'][k] == 'True':
	            ax = axs[n]
                    good = np.where(~np.isnan(ds[k]))[0]
                    pc = ax.plot(ds.time[good], ds[k][good], '.')
                    min, max = _autoclim(ds[k][good])
                    ax.set_ylim(min, max)
                    ax.set_title(ds[k].attrs['long_name'] + ' [' +
                             ds[k].attrs['units'] + ']', loc='left', fontsize=9)

            fig.savefig(config['figdir'] + '/ts_%s.png'%ds.attrs['deployment_name'], dpi=200)

        # colorline?
        fig, axs = plt.subplots(int(N / 2), 2, figsize=(7.5, 7),
                                sharex=True, sharey=True)
        axs = axs.flat
<<<<<<< HEAD
        _log.info('Plotting colorline data')

        for n, k in enumerate(keys):
            print('key', k)
            ax = axs[n]
            locator = mdates.AutoDateLocator(minticks=3, maxticks=7)
            formatter = mdates.ConciseDateFormatter(locator)
            ax.xaxis.set_major_locator(locator)
            ax.xaxis.set_major_formatter(formatter)
            good = np.where(~np.isnan(ds[k]))[0]
            min, max = _autoclim(ds[k][good])
            pc = ax.scatter(ds.time[good], ds.depth[good], s=3, c=ds[k][good],
                rasterized=True, vmin=min, vmax=max)

            fig.colorbar(pc, ax=ax, extend='both', shrink=0.6)
            ax.set_title(ds[k].attrs['long_name'] + ' [' +
                         ds[k].attrs['units'] + ']', loc='left', fontsize=9)
            t0 = ds.time[good][0]
            t1 = ds.time[good][-1]
            ax.set_xlim([t0, t1])
            ax.set_ylim([ds['depth'].max(), ds['depth'].min()])
            if n == 0:
                ax.set_ylabel('DEPTH [m]')
        fig.savefig(config['figdir'] + '/cl_%s.png'%ds.attrs['deployment_name'], dpi=200)

        # depth profiles...
        fig, axs = plt.subplots(int(N / 2), 2, figsize=(7.5, 7),
                                sharex=True, sharey=True)
        axs = axs.flat
        _log.info('Plotting colorline data')

        for n, k in enumerate(keys):
            _log.info('cl %s', k)
            if config['timeseries'][k] == 'True':
                ax = axs[n]
                locator = mdates.AutoDateLocator(minticks=3, maxticks=7)
                formatter = mdates.ConciseDateFormatter(locator)
                ax.xaxis.set_major_locator(locator)
                ax.xaxis.set_major_formatter(formatter)
                good = np.where(~np.isnan(ds[k] + ds.depth))[0]
                min, max = _autoclim(ds[k][good])
                pc = ax.scatter(ds.time[good], ds.depth[good], s=3, c=ds[k][good],
                    rasterized=True)
                fig.colorbar(pc, ax=ax, extend='both', shrink=0.6)
                ax.set_title(ds[k].attrs['long_name'] + ' [' +
                             ds[k].attrs['units'] + ']', loc='left', fontsize=9)
                t0 = ds.time[good][0]
                t1 = ds.time[good][-1]
                ax.set_xlim([t0, t1])
                ax.set_ylim([ds['depth'].max(), ds['depth'].min()])
                if n == 0:
                    ax.set_ylabel('DEPTH [m]')
        fig.savefig(config['figdir'] +
                    '/cl_%s.png'%ds.attrs['deployment_name'], dpi=200)

        # prop_v_prop:
        _log.info('property vs property')

        fig, axs = plt.subplots(1, 2, constrained_layout=True, sharey=True)
        ax = axs[0]
        s = np.linspace(ds.salinity.min(), ds.salinity.max(), 100)
        t = np.linspace(ds.temperature.min(), ds.temperature.max(), 100)
        S, T = np.meshgrid(s, t)
        pd = seawater.eos80.pden(S, T, 0, 0) - 1000
        levels = np.arange(20, 30, 0.5)
        c = ax.contour(s, t, pd, colors='0.5', levels=levels)
        ax.clabel(c, levels[::2], fontsize=8, fmt='%1.1f')
        ax.plot(ds['salinity'], ds['temperature'], '.', markersize=2)
        ax.set_xlabel('$S\ [psu]$')
        ax.set_ylabel('$T\ [^oC]$')

        try:
            ax = axs[1]

            ax.plot(ds['oxygen_concentration'], ds['temperature'],
                    '.', markersize=2)
            ax.set_xlabel('$O^2\ [mmol/L]$')
        except:
            pass
        add_suptitle(fig, ds)
        fig.savefig(config['figdir'] +
            '/pvp_%s.png'%ds.attrs['deployment_name'], dpi=200)

def add_suptitle(fig, ds):
    sst = 'Glider/deployment: %s; ' % ds.attrs['deployment_name']
    sst += 'lat = %1.3f to %1.3f N; ' % (ds.attrs['geospatial_lat_min'],
                                  ds.attrs['geospatial_lat_max'])
    sst += 'lon = %1.3f to %1.3f E ' % (ds.attrs['geospatial_lon_min'],
                                  ds.attrs['geospatial_lon_max'])
    sst += '\n CPROOF: http://cproof.uvic.ca/; Preliminary data'

    fig.suptitle(sst, fontsize=8, fontfamily='courier')


def grid_plots(fname, plottingyaml):
    with open(plottingyaml) as fin:
        config = yaml.safe_load(fin)

    try:
        os.mkdir(config['figdir'])
    except:
        pass

    with xr.open_dataset(fname, decode_times=True) as ds:
        keys = config['timeseries'].keys()
        N = len(keys)
        # get the max depth that data is at:
        tmean = ds.temperature.mean(axis=1)
        indmax = np.where(~np.isnan(tmean))[0][-1]
        depmax = ds.depth[indmax]
        fig, axs = plt.subplots(int(N / 2), 2, figsize=(7.5, 7),
                                sharex=True, sharey=True)
        axs = axs.flat
        for n, k in enumerate(keys):
            print('key', k)
            if config['timeseries'][k] == 'True':
                ax = axs[n]
                locator = mdates.AutoDateLocator(minticks=3, maxticks=7)
                formatter = mdates.ConciseDateFormatter(locator)
                ax.xaxis.set_major_locator(locator)
                ax.xaxis.set_major_formatter(formatter)
                min, max = _autoclim(ds[k])
                print('min, max', min, max)
                pc = ax.pcolormesh(ds.time, ds.depth, ds[k],
                    rasterized=True, vmin=min, vmax=max)
                print(ds[k])

                fig.colorbar(pc, ax=ax, extend='both', shrink=0.6)
                ax.set_title(ds[k].attrs['long_name'] + ' [' +
                         ds[k].attrs['units'] + ']', loc='left', fontsize=9)
                t0 = ds.time[0]
                t1 = ds.time[-1]
                ax.set_xlim([t0, t1])
                ax.set_ylim([maxdepth, 0])
                if n == 0:
                    ax.set_ylabel('DEPTH [m]')
        fig.savefig(config['figdir'] + '/pcolor_%s.png'%ds.attrs['deployment_name'], dpi=200)
