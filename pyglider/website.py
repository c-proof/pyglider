import yaml
import glob
import os
import xarray as xr
import numpy as np
from jinja2 import Environment, FileSystemLoader
import geojson
import datetime
import simplekml
import pyglider.utils as pygu

import logging

"""
Utilities to make smart directories for the websites.

A driving script may look like

```
import pyglider.website as pyweb
import logging

logging.basicConfig(level=logging.DEBUG)
# uses the templat in `.templates/deploymentsIndex.html` to render
# index.html in each of the subdirectories:

pyweb.geojson_deployments('./')

if 1:
    pyweb.index_deployments('./dfo-walle652/')
    pyweb.index_deployments('./dfo-bb046/')
```

where `dfo-walle652` will contain subdirectories, each with a deployment:

```
./dfo-walle652/dfo-walle652-20210903
./dfo-walle652/dfo-walle652-20220101
...
```

Output are index files in each subdirectory and a geojson of all
the deployments (if that is desired)

```
./deployments.geojson   # from the first script....
./dfo-walle652/index.html
./dfo-walle652/dfo-walle652-20210903/index.html
...

Requires an html template.  The ones used for cproof are in
../example_html_templates of this repo and you need
to supply

`deploymentsIndex.html`
`deploymentsInfoNew.html`
and maybe `deploymentsInfo.html`

These are example templates, and written in jinga templating language.  In
particular, see the strings near the bottom with double braces: eg:
"{{ att['deployment_name'] }}"  These are what get filled by calling
`pyweb.index_deployments`.

Note there are some hardcoded things in this tool, so please read the source
code!

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
                _log.info(d)
                nc = glob.glob(d+'/L0-timeseries/*.nc')
                if len(nc) < 1:
                    nc = glob.glob(d+'/L1-timeseries/*.nc')
                    if len(nc) < 1:
                        continue
                with xr.open_dataset(nc[0]) as ds:
                    att = ds.attrs
                    atts.append(att)
            else:
                pass
    if len(atts) > 0:
        output = template.render(atts=atts,
            title=atts[-1]['glider_name'] + atts[-1]['glider_serial'])
        with open(dir + '/index.html', 'w') as fout:
            fout.write(output)

    # now make the individual deployment index pages...
    templateOld = env.get_template('deploymentsInfo.html')
    templateNew = env.get_template('deploymentsInfoNew.html')
    subdirs = glob.glob(dir + '/*')
    atts = []
    for d in subdirs:
        if os.path.isdir(d):
            _log.info(d)
            if 1:
                figs = glob.glob(d + '/figs/*.png')
                for n, fig in enumerate(figs):
                    figs[n] = './figs/' + os.path.split(fig)[1]
                nc = glob.glob(d+'/L0-timeseries/*.nc')
                template = templateNew

                if len(nc) < 1:
                    # try old style
                    nc = glob.glob(d+'/L1-timeseries/*.nc')
                    template = templateOld
                    if len(nc) < 1:
                        continue
                with xr.open_dataset(nc[0]) as ds:
                    att = ds.attrs
                    _log.debug(type(ds.keys()))
                    keys = []
                    units = []
                    for key in ds.keys():
                        _log.debug(ds[key].attrs)
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

    kml = simplekml.Kml()

    np.random.seed(20190101)
    _log.debug(f'subdirs, {subdirs}')
    colornum = 0;
    for d in subdirs:
        _log.info(d)
        if os.path.isdir(d):
            subdirs2 = glob.glob(d + '/*')
            for d2 in subdirs2:
                _log.info(d2)
                if os.path.isdir(d2):
                    try:
                        nc = glob.glob(d2+'/L0-gridfiles/*.nc')
                        if len(nc) < 1:
                            # old style
                            nc = glob.glob(d2+'/L2-gridfiles/*.nc')

                        with xr.open_dataset(nc[0]) as ds:
                            _log.info(f'opened {nc[0]}')
                            att = ds.attrs
                            good = (ds.longitude < -125)
                            line = np.vstack((ds.longitude[good], ds.latitude[good])).T
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
                            # cols = pygu.get_html_non_blue(colornum)
                            colornum += 1
                            feat.properties['color'] = '#%02X%02X%02X' % (cols[0], cols[1], cols[2])
                            if ds['time'][-1] > np.datetime64(datetime.datetime.now()) - np.timedelta64(2, 'D'):
                                feat.properties['active'] = True
                            else:
                                feat.properties['active'] = False

                            features += [feat]

                            # make the kml:
                            pnt = kml.newpoint(coords=[line[-1]])
                            pnt.style.iconstyle.icon.href = 'http://cproof.uvic.ca/deployments/assets/images/slocum_glider.png'
                            coords = []
                            for thelon, thelat  in zip(ds.longitude.values, ds.latitude.values):
                                coords += [(thelon, thelat)]
                            pnt.timestamp.when = f'{ds.time.values[-1]}'[:-3]
                            ls = kml.newlinestring(coords=coords,
                                name=att['deployment_name'],
                                )
                            ls.timespan.begin = f'{ds.time.values[0]}'[:-3]
                            ls.timespan.end = f'{ds.time.values[-1]}'[:-3]
                            ls.style.linestyle.color = 'ee' + '%02X%02X%02X' %  (cols[2], cols[1], cols[0])
                            ls.style.linestyle.width = 3;
                            kml.save(d2[2:]+'/'+att['deployment_name']+'.kml')

                    except:
                        _log.info(f'Could not find grid file {d2}')
    feature_collection = geojson.FeatureCollection(features)
    with open(outfile, 'w') as fout:
        s = geojson.dumps(feature_collection)
        fout.write(s)
