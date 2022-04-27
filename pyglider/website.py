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

# jpnote addition for argo float calcs 
import argopandas as argo 
import pandas as pd
# end of jp addition 

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

def ticker(dir):

#jpnote: 
# Function to calculate distance traveled by gliders 
# For every mission directory where L0-timeseries exists,
# extract distance_over_ground and sum 
# Written into html file, and put directly into cproofwesbite/_includes/ 
# directory, and then added using jekyll {% include %}

    sump = 0                   # set variable to 0
    subdirs = glob.glob(dir + '/dfo-*') 
   # atts = []
    for d in subdirs:
        if os.path.isdir(d):   #define subdirs so that only go into directories that start with 'dfo-*'
            if 1:
                _log.info(d) 
                nc = sorted(glob.glob(d+'/*/L0-timeseries/*.nc'), key=os.path.getmtime)
                with xr.open_dataset(nc[-1]) as ds:  #open most recent of sorted netcdf files (last, [-1])
                    
                    sump += ds.variables['distance_over_ground'].data[-1]                   # Add last element of dataset to final sum
                    list_ = np.where(ds.variables['distance_over_ground'].data==0)          # Locate where in distance over ground there is a 0, 
                    flattened = [val for sublist in list_ for val in sublist]               # List comprhension to allow traversal through list 
           
                    for i in flattened:
                        sump += ds.variables['distance_over_ground'].data[i-1]              # add number before every 0 appearance to final sum 

   # print('final sum', "{:.f}".format(sump)) #final sum print check 
   # print('final sum', int(sump)) #final sum print check 

    with open('/Users/cproof/cproofwebsite/_includes/inputfile.html', 'w') as output_file:  #output to html file in website directory
        output_file.write(str(int(sump)))

#######


def url_jp(dir):  

    """
    jpnote: Function to create files with urls from L0-timeseries and 
    L0-gridfiles in individual glider directories 
    for  Wget of data from WGET DATA Page 
    accessible from the DATA Page on the C-PROOF website. 

    """

    list_dir = []
    # Inital open and write to overwrite any existing text upon rerun of script. 
    # Otherwise only changes will be ammended 
    with open(dir + '/LineP.txt','w') as f:         # Clear the file so that running the script doest append previous things to it.  
        pass
    with open(dir + '/CalvertLine.txt','w') as f:   # Clear the file so that running the script doest append previous things to it.
        pass
    with open(dir + '/mission_all.txt','w') as f:   # Clear the file so that running the script doest append previous things to it.  
        pass
    # To add a new line, include:  
    #with open(dir + '/NewLineName.txt','w') as f:   # Clear the file so that running the script doest append previous things to it.
        #pass

   
    subdirs = glob.glob(dir + '/dfo-*')             # Only enter into subdirectories beginnign with "dfo-""
    for d in subdirs:
        if os.path.isdir(d):
            list_dir.append(d)                      # The directories 
            
            if 1:   #all 
                _log.info(d)
                tnc = sorted(glob.glob(d+'/*/L0-timeseries/*.nc'), key=os.path.getmtime)    # Extract all files in L0-timeseries, and sort by time created 
                tnc = tnc[-1]                       # Only take path from most recent file created
                gnc = glob.glob(d+'/*/L0-gridfiles/*.nc')    # Extract file from L0-gridfiles directory 

             
                with xr.open_dataset(gnc[0]) as ds:
                    comment = ds.comment             # Extract comment from file in L0-gridfiles directory 
                    
                list_jp = []
                tncq = []
                tncq.append(tnc)                     # Convert to list so that extend works properly (because extracting [-1] changes form) (bug fix)
                list_jp.extend(tncq)
                list_jp.extend(gnc)

               
                mystring = 'https://cproof.uvic.ca/gliderdata/deployments' 

                list_jp = [s[1:] for s in list_jp]
                list_jp = [mystring + x for x in list_jp]



                if 'Line P' in ds.comment:
                    with open(dir + 'LineP.txt', 'a') as fout:
                        fout.write('\n'.join(str(line) for line in list_jp))
                        fout.write('\n')
                if 'Calvert' in ds.comment:
                    with open(dir + 'CalvertLine.txt', 'a') as fout:
                        fout.write('\n'.join(str(line) for line in list_jp))
                        fout.write('\n')
                
                # To add new line, include

                #if 'NewLineName' in ds.comment:
                #   with open(dir + 'NewLineName.txt', 'a') as fout:
                #       fout.write('\n'.join(str(line) for line in list_jp))
                #       fout.write('\n')

                with open(dir + d + '.txt', 'w') as fout:
                    fout.write('\n'.join(str(line) for line in list_jp))
                    fout.write('\n')
                
                with open(dir + '/mission_all.txt', 'a') as fout:
                    fout.write('\n'.join(str(line) for line in list_jp))
                    fout.write('\n')

                    fout.close()

##############

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
                nc = sorted(glob.glob(d+'/L0-timeseries/*.nc'), key=os.path.getmtime) # jpnote: Sort files by time created
                if len(nc) < 1:
                    nc = glob.glob(d+'/L1-timeseries/*.nc')
                    if len(nc) < 1:
                        continue
                with xr.open_dataset(nc[-1]) as ds:          # jpnote: Open most recently created file. Last file in list after sort
                    att = ds.attrs
                    atts.append(att)
                    #print(ds.attrs)
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
                nc = sorted(glob.glob(d+'/L0-timeseries/*.nc'), key=os.path.getmtime) # jpnote: Sort files by time created
               # nc = glob.glob(d+'/L0-timeseries/*.nc') # Old method 
                template = templateNew

                if len(nc) < 1:
                    # try old style
                    nc = glob.glob(d+'/L1-timeseries/*.nc')
                    template = templateOld
                    if len(nc) < 1:
                        continue
                with xr.open_dataset(nc[-1]) as ds:  # jpnote: Open most recently created file. 
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

    #kml = simplekml.Kml()   #jpnote: commented, open kml file in loop 

    np.random.seed(20190101)
    _log.debug(f'subdirs, {subdirs}')
    colornum = 0;
    for d in subdirs:
        _log.info(d)
        if os.path.isdir(d):
            subdirs2 = glob.glob(d + '/*')
            for d2 in subdirs2:
                _log.info(d2)
                kml = simplekml.Kml()  #jpnote: kml object in for loop
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
                          # jpnote: URL depending on method of testing 
                            # get URL....
                          # For cproof website:
                            feat.properties['url'] = ('' +
                                'http://cproof.uvic.ca/gliderdata/deployments/' +  
                                d2[2:])
                          # For local testing: 
                           # feat.properties['url'] = ('' +
                           #     'http://localhost:4000/gliderdata/deployments/' +
                           #     d2[2:])
                          # For git pages testing: 
                            #feat.properties['url'] = ('' +
                            #    'http://juliaputko.github.io/cproofwebsite/gliderdata/deployments/' +
                            #    d2[2:])
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
######
# jpnote: added to different geojson file: Glider + ARGO for the C-PROOF website home page 
######
    floatsmeds_bgc_floats = argo.bio_prof[argo.bio_prof['file'].str.contains('meds')]
    my_wmos = [4902549,4902550, 4902551, 4902552,4902553,4902554,4902555,4902583,4902584,4902585,4902586,4902587,4902588,4902589]
   # Adding additional ARGO float numbers 

    df = pd.DataFrame()

    coords = []
    for f in argo.float(my_wmos):
        for row in f.prof:
           # if f.prof.latitude > 0:  #make not nan
            lat = f.prof.latitude
            lon = f.prof.longitude
         
            line = np.vstack((f.prof.longitude, f.prof.latitude)).T
              # line = np.vstack((f.prof.longitude, f.prof.latitude)).T
            mask = np.all(np.isnan(line), axis=1)
            line = line[~mask]
            ls = geojson.LineString(line.tolist())
            feat = geojson.Feature(geometry=ls)
               
            # get color:
            cols = np.random.randint(0, 200, 3)
            colornum += 1
            feat.properties['color'] = '#%02X%02X%02X' % (cols[0], cols[1], cols[2])
            feat.properties['name'] = 'argo'
          #  if ds['time'][-1] > np.datetime64(datetime.datetime.now()) - np.timedelta64(2, 'D'):
           #     feat.properties['active'] = True
          #  else:
           #     feat.properties['active'] = False  ?? current status ?? 

        features += [feat]
#### 

    #jpnote: Addition of separate outfile that contains argo float data as well as glider mission data
    feature_collection = geojson.FeatureCollection(features)
    with open('cproof-deployments_all.geojson', 'w') as fout:
        s = geojson.dumps(feature_collection)
        fout.write(s)


## GLOBAL VARIABLES STORY PAGE CREATION
all_c = []
link_c = []
figs_c= []

all_lp = []
link_lp = []
figs_lp= []

#new line addition 
#allall_newline = []
#link_newline = []
#figs_newline = []

allall_lp = []
allall_cal = []

def index_story(dir, templatedir='./.templates/'):

    """
    jpnote: function to create Story pages for each line. 
    For additions to Story pages, if new line is created, 
    Line page is hardcoded, but outline can be copy pasted. 

    """

    global all_c 
    global link_c
    global figs_c
    global all_lp
    global link_lp
    global figs_lp

    global allall_lp
    global allall_cal


    file_loader = FileSystemLoader(templatedir)
    env = Environment(loader=file_loader)
    template = env.get_template('story.html')

    # now make the individual deployment index pages...
    #templateOld = env.get_template('deploymentsInfo.html')
    templateNew = env.get_template('story.html')
    subdirs = glob.glob(dir + '/dfo-*')  #if subdirs are "df0-" missions   #dir is './'
       
   # print('jp_subdirs',subdirs)
    atts_lp = []
    atts_c = []
    for d in subdirs:
        if os.path.isdir(d):  #dir and subdires 
            subdirs2 = glob.glob(d + '/dfo-*')
            for d2 in subdirs2:
                if 1:
                    figs = glob.glob(d2 + '/figs/*.png')
                    for n, fig in enumerate(figs):
                       figs[n] = '/figs/' + os.path.split(fig)[1]

                   
                    nc = sorted(glob.glob(d2+'/L0-timeseries/*.nc'), key=os.path.getmtime)
            
                    template = templateNew
               
                    with xr.open_dataset(nc[-1]) as ds:  

                        if 'Line P' in ds.comment:
                            allall_lp.append(nc[-1])
                            att_lp = ds.attrs
                            _log.debug(type(ds.keys()))
                            keys_lp = []
                            units_lp = []
                            for key in ds.keys():
                                _log.debug(ds[key].attrs)
                                try:
                                    unit_lp = ds[key].attrs['units']
                                except KeyError:
                                    unit_lp = 'no units'
                                keys_lp.append(key + ' [' + unit_lp +']')
                            depname_lp = att_lp['deployment_name']
                            name_lp = att_lp['glider_name'] + att_lp['glider_serial']   # add the + glider serial once bb's naming is fixed
                         
                            all_lp.append(depname_lp)
                            link_lp.append(name_lp)
                        
                ## the bb missions :: 
                        if 'Calvert' in ds.comment:
                            allall_cal.append(nc[-1])
                            att_c = ds.attrs
                            _log.debug(type(ds.keys()))
                            keys_c = []
                            units_c = []
                        
                           
                            for key in ds.keys():
                                _log.debug(ds[key].attrs)
                               
                                try:
                                    unit_c = ds[key].attrs['units']
                                except KeyError:
                                    unit_c = 'no units'
                                keys_c.append(key + ' [' + unit_c +']')
                            depname_c = att_c['deployment_name']
                            name_c = att_c['glider_name']
                         
                            all_c.append(depname_c)
                            link_c.append(name_c)

                            #comment
                            figs_c = glob.glob(d2 + '/figs/*.png')
                        
                            for n, fig in enumerate(figs):
                                figs_c[n] = '/figs/' + os.path.split(fig)[1]  
                       
        # ADDING A NEW LINE: 
        # New Line Addition: if ''
                    #if 'NewLineName' in ds.comment:

                           # att_newline = ds.attrs
                           # _log.debug(type(ds.keys()))
                           # keys_newline = []
                           # units_newline = []
                          #  for key in ds.keys():
                              #  _log.debug(ds[key].attrs)
                           #     try:
                              #      unit_newline = ds[key].attrs['units']
                              #  except KeyError:
                              #      unit_newline = 'no units'
                             #   keys_newline.append(key + ' [' + unit_newline +']')
                         #   depname_newline = att_newline['deployment_name']
                        #    name_newline = att_newline['glider_name']
                         
                        #    all_newline.append(depname_newline)
                        #    link_newline.append(name_newline)
                        #    figs_newline = glob.glob(d2 + '/figs/*.png')
                        #    for n, fig in enumerate(figs):
                        #        figs_newline[n] = '/figs/' + os.path.split(fig)[1]                           
    ### end of new line addition copy paste 


    # extract figures only from most recent figures for each line 
    timelp = []
    timecal = []
    for l in allall_lp:
        with xr.open_dataset(l) as lp: 
            timelp.append((lp.time_coverage_start, l))
    tmax = max(timelp)
    mrecent = timelp.index(tmax)
    path = timelp[mrecent][1]     
    path = path.split('/L0-timeseries')     
    if 1:
        figs_lp = glob.glob(path[0] + '/figs/*.png')
        for n, fig in enumerate(figs_lp):
            figs_lp[n] = '/figs/' + os.path.split(fig)[1]
     
    new = path[0].split('/')
    mission_lp = new[2]
    

#EXTRACT MOST RECENT MISSIONS FROM CALVERT LINE  
    timecal = []
    for l in allall_cal:
        with xr.open_dataset(l) as lp: 
            timecal.append((lp.time_coverage_start, l))
    tmax = max(timecal)
    mrecent = timecal.index(tmax)
    path = timecal[mrecent][1]                          #second element of the tuple <-- the nice and then open that file and           
    path = path.split('/L0-timeseries')
    if 1:
        figs_cal = glob.glob(path[0] + '/figs/*.png')
        for n, fig in enumerate(figs_cal):
            figs_cal[n] = '/figs/' + os.path.split(fig)[1]

    new = path[0].split('/')
    mission_cal = new[2]

#output line p info to template
    output_lp = template.render(deploy_name=depname_lp,
                 title=att_lp['comment'],
                 figs=figs_lp, att=att_lp, keys=keys_lp,all_=all_lp, line_=att_lp['comment'], link_=link_lp[0], mission=mission_lp)  #dir_figs=dirfigs)

    with open(dir+'/linep.html', 'w') as fout:
        fout.write(output_lp)

#output calvert line info to template

    output_c = template.render(deploy_name=depname_c,
        title=att_c['deployment_name'] ,
        figs=figs_c, att=att_c, keys=keys_c, all_=all_c, line_=att_c['comment'], link_=link_c[0], mission=mission_cal)  #+ att_c['glider_serial']
        
    with open(dir+'/calvert.html', 'w') as fout:
        fout.write(output_c)

    ## New Line addition: 

    # If you'd to make additions from the netcdf files to the Story Pages, 
    # this is where the content is passed into the template

    #output_newline = template.render(deploy_name=depname_newline,
    #    title=att_newline['deployment_name'] ,
    #    figs=figs_newline, att=att_newline, keys=keys_newline, all_=all_newline, line_=att_newline['comment'], link_=link_newline[0])  #+ att_c['glider_serial']
        
    #with open(dir+'/NewLineName.html', 'w') as fout:
    #    fout.write(output_newline)

  
