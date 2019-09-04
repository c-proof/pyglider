import yaml
import glob
import os
import xarray as xr
import numpy as np
from jinja2 import Environment, FileSystemLoader
import geojson
import datetime

import logging

"""
Utilities to make smart directories for the websites.  Dumb persons erdapp
server...
"""

_log = logging.getLogger(__name__)

def index_deployments(dir, templatedir='./.templates/'):
    """
    Get useful info from deployments under "dir" and add to an
    index.html in "dir"

    The structure is meant to be simple:

    dir/deployment1/deployment.yml
    dir/deployment2/deployment.yml
    dir/deployment3/deployment.yml

    Makes a file `dir/index.html` that has table of the deployments.
    """

    file_loader = FileSystemLoader(templatedir)
    env = Environment(loader=file_loader)
    template = env.get_template('deploymentsIndex.html')


    subdirs = glob.glob(dir + '/*')
    atts = []
    for d in subdirs:
        if os.path.isdir(d):
            if 1:
                print(d)
                nc = glob.glob(d+'/L1-timeseries/*.nc')
                with xr.open_dataset(nc[0]) as ds:
                    att = ds.attrs
                    atts.append(att)
            else:
                pass
    output = template.render(atts=atts,
        title=atts[-1]['glider_name'] + atts[-1]['glider_serial'])
    with open(dir + '/index.html', 'w') as fout:
        fout.write(output)

    # now make the individual deployment index pages...
    template = env.get_template('deploymentsInfo.html')
    subdirs = glob.glob(dir + '/*')
    atts = []
    for d in subdirs:
        print(d)
        if os.path.isdir(d):
            print(d)
            if 1:
                nc = glob.glob(d+'/L1-timeseries/*.nc')
                figs = glob.glob(d + '/figs/*.png')
                for n, fig in enumerate(figs):
                    figs[n] = './figs/' + os.path.split(fig)[1]
                with xr.open_dataset(nc[0]) as ds:
                    att = ds.attrs
                    print(type(ds.keys()))
                    keys = []
                    units = []
                    for key in ds.keys():
                        print(ds[key].attrs)
                        try:
                            unit = ds[key].attrs['units']
                        except KeyError:
                            unit = 'no units'
                        keys.append(key + ' [' + unit +']')
                depname = att['deployment_name']
                output = template.render(deploy_name=depname,
                    title=att['glider_name'] + att['glider_serial'],
                    figs=figs, att=att, keys=keys)
                with open(d + '/index.html', 'w') as fout:
                    fout.write(output)
            else:
                pass

def geojson_deployments(dir, outfile='cproof-deployments.geojson'):
    props = ['deployment_start', 'deployment_end', 'platform_type',
             'glider_model', 'glider_name', 'glider_serial',
             'deployment_name', 'project', 'institution', 'comment']
    subdirs = glob.glob(dir + '/*')
    features = []
    np.random.seed(20190101)

    for d in subdirs:
        if os.path.isdir(d):
            subdirs2 = glob.glob(d + '/*')
            for d2 in subdirs2:
                if os.path.isdir(d2):
                    try:
                        nc = glob.glob(d2+'/L2-gridfiles/*.nc')
                        with xr.open_dataset(nc[0]) as ds:
                            print(ds)
                            att = ds.attrs
                            line = np.vstack((ds.longitude, ds.latitude)).T
                            ls = geojson.LineString(line.tolist())
                            feat = geojson.Feature(geometry=ls)
                            for prop in props:
                                if prop in ds.attrs.keys():
                                    feat.properties[prop] = ds.attrs[prop]
                                else:
                                    feat.properties[prop] = ''

                            # get URL....
                            feat.properties['url'] = ('' +
                                'http://cproof.uvic.ca/gliderdata/deployments/' +
                                d2[2:])
                            # get color:
                            cols = np.random.randint(0, 200, 3)
                            print(cols)
                            feat.properties['color'] = '#%02X%02X%02X' % (cols[0], cols[1], int(cols[2] / 2))
                            if ds['time'][-1] > np.datetime64(datetime.datetime.now()):
                                feat.properties['active'] = 'True'
                            else:
                                feat.properties['active'] = 'False'

                            features += [feat]


                    except:
                        _log.info(f'Could not find grid file {d2}')
    feature_collection = geojson.FeatureCollection(features)
    with open(outfile, 'w') as fout:
        s = geojson.dumps(feature_collection)
        fout.write(s)
