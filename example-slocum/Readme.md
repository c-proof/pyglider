

## Get raw data:
```
% mkdir ../SL713_0001/binary/
% cp ../card_offloads/rosie_713/card_offload_20190313_rosie713/Science/SENTLOGS/*.EBD binary
% cp ../card_offloads/rosie_713/card_offload_20190313_rosie713/Science/LOGS/*.EBD binary
% cp ../card_offloads/rosie_713/card_offload_20190313_rosie713/Main_board/SENTLOGS/*.EBD binary
% cp ../card_offloads/rosie_713/card_offload_20190313_rosie713/Main_board/LOGS/*.EBD binary
```

You can organize these how you like.  I think for realtime, I'm putting in `raw_realtime`


### Get the headers:

```
mkdir cac
cp ../card_offloads/rosie_713/card_offload_20190313_rosie713/Science/STATE/CACHE/*.CAC cac
```

Or get these from SFMC: `Configuration/Cache Files` and put them in `cac`


### Get and edit the sensor filter file:

```
./dfo_rosie713_sensors.txt
```

## Edit deployment info:

Edit this to have correct glider name, and deployment dates...
`deployment.yml`

## Edit `process_deploymentRealtime.py`

Note that you will want to edit this to have correct directories:
Also, realtime and delayed time have different suffixies:

```
scisuffix    = 'tbd'
glidersuffix = 'sbd'
```
