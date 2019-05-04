import yaml
import glob
import os
import xarray as xr
from jinja2 import Environment, FileSystemLoader

"""
Utilities to make smart directories for the websites.  Dumb persons erdapp
server...
"""


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
                nc = glob.glob(d+'/L1-timeseries/*.nc')
                with xr.open_dataset(nc[0]) as ds:
                    att = ds.attrs
                    atts.append(att)
            else:
                pass
    output = template.render(atts=atts,
        title=atts[0]['glider_name'] + atts[0]['glider_serial'])
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
                depname = att['title']
                output = template.render(deploy_name=depname,
                    title=att['glider_name'] + att['glider_serial'],
                    figs=figs, att=att, keys=keys)
                with open(d + '/index.html', 'w') as fout:
                    fout.write(output)
            else:
                pass
