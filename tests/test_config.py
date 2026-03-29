from pathlib import Path
from typing import Callable

import pytest
import yaml

from pyglider.config import Deployment

library_dir = Path(__file__).parent.parent.absolute()
example_dir = library_dir / 'tests/example-data/'

VALID_CONFIG_PATHS = [
    example_dir / 'example-slocum/deploymentRealtime.yml',
    example_dir / 'example-seaexplorer/deploymentRealtime.yml',
    # TODO: Add other valid example configs?
]


@pytest.fixture(params=VALID_CONFIG_PATHS)
def valid_config(request):
    return request.param.read_text()


def test_valid_config(valid_config):
    """Checks all valid configs can be loaded."""
    Deployment.load_yaml(valid_config)


def patch_config(yaml_str: str, f: Callable) -> str:
    """Patch a yaml string with a function.
    Function should do an in place operation on a dictionary/list.
    """
    d = yaml.safe_load(yaml_str)
    f(d)
    return yaml.safe_dump(d, default_flow_style=False)


@pytest.mark.parametrize(
    'input_, expected, f',
    [
        (
            'a: 1\n' 'b: 2\n',
            'a: 1\n' 'b: 3\n',
            lambda d: d.update({'b': 3}),
        ),
    ],
)
def test_patch_config(input_, expected, f):
    assert patch_config(input_, f) == expected


# TODO: Stress test the model by taking existing configs and modifying them in breaking ways


def test_incorrect_date_format(): ...
