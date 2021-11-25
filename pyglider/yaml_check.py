import yaml
import sys
from datetime import datetime, timedelta
from urllib import request


def check_yaml(yaml_path, check_urls=False):
    """
    Simple yaml checker for expected SeaExplorer deployment yaml.
    THIS IS NOT A SUBSTITUTE FOR CHECKING YOUR YAML!
    But it will flag some values you might have missed.

    Parameters
    ----------
    yaml_path: path to your yaml file
    check_urls: boolean. If True, checks that urls are reachable. Default: False
    """
    failures = 0
    print(f'Checking deployment yaml at:\n{yaml_path}')
    with open(yaml_path) as fin:
        deployment = yaml.safe_load(fin)
    print('read yaml successfully')
    print('Checking top level items')
    for item in ['metadata', 'glider_devices', 'netcdf_variables', 'profile_variables']:
        if item not in deployment.keys():
            print(f'ERROR: {item} not found')
            failures += 1

    print('Checking metadata')
    metadata_keys = (
        'acknowledgement', 'comment', 'contributor_name', 'contributor_role', 'creator_email', 'creator_name',
        'creator_url', 'deployment_id', 'deployment_name', 'deployment_start', 'deployment_end', 'format_version',
        'glider_name', 'glider_serial', 'glider_model', 'glider_instrument_name', 'glider_wmo', 'institution',
        'keywords',
        'keywords_vocabulary', 'license', 'metadata_link', 'Metadata_Conventions', 'naming_authority', 'platform_type',
        'processing_level', 'project', 'project_url', 'publisher_email', 'publisher_name', 'publisher_url',
        'references',
        'sea_name', 'source', 'standard_name_vocabulary', 'summary', 'transmission_system', 'wmo_id')

    for key in metadata_keys:
        if key not in deployment['metadata'].keys():
            print(f'ERROR {key} not found in metadata')
            failures += 1

    print('Checking dates')
    start = deployment['metadata']['deployment_start']
    end = deployment['metadata']['deployment_end']
    try:
        start_time = datetime.strptime(start, "%Y-%m-%d")
    except ValueError:
        print(f'ERROR deployment_start {start} incorrectly formatted. Should be YYYY-MM-DD')
        failures += 1
    try:
        end_time = datetime.strptime(end, "%Y-%m-%d")
    except ValueError:
        print(f'ERROR deployment_end {end} incorrectly formatted. Should be YYYY-MM-DD')
        failures += 1
    try:
        deployment_duration = end_time - start_time
        if deployment_duration < timedelta(0):
            print('ERROR deployment_end date is sooner than deployment_start date')
            failures += 1
        if deployment_duration > timedelta(days=365):
            print('Warning inferred deployment duration > 1 year, please check deployment_start and deployment_end')
    except ValueError:
        pass
    if check_urls:
        print('Checking urls')
        for url_id in ['creator_url', 'project_url', 'publisher_url']:
            url = deployment['metadata'][url_id]
            try:
                http_code = request.urlopen(url).getcode()
                if int(http_code / 100) != 2:
                    print(f'Warning, did not receive 200 html response from {url}')
            except ValueError:
                print(f"ERROR could not reach {url_id}: {url}")
                failures += 1
    else:
        print(
            'Warning, not checking urls. Enable this time-consuming check by calling:\n'
            'python yaml_check.py deployment.yaml True')

    print('Checking glider_devices')
    devices = deployment['glider_devices']
    for name, device in devices.items():
        for field in ('make', 'model', 'serial'):
            if field not in device.keys():
                print(f'ERROR {field} not present for glider_devices: {name}')
                failures += 1

    failures = check_strings(deployment, failures)
    if failures == 0:
        print('\nYour yaml looks pretty good to me :)\nyou should still check it though')
    else:
        print(f'\nYour yaml seems to have {failures} missing/erroneous values, look for ERROR statements printed above')
    return True


def check_strings(d, failures=0):
    for k, v in d.items():
        if isinstance(v, dict):
            failures = check_strings(v, failures)
        else:
            if type(v) is not str:
                print(f'Error, value of {k}: {v} is not a string')
                failures += 1
            if not bool(v):
                print(f'Error, {k} is empty')
                failures += 1
    return failures


if __name__ == '__main__':
    args = sys.argv
    url_check = False
    if len(args) > 1:
        yaml_file = args[1]
    else:
        sys.exit('Must specify yaml file to parse')
    if len(args) > 2:
        url_check = args[2]
    check_yaml(yaml_file, check_urls=url_check)
