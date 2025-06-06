import yaml
from pydantic import BaseModel

from pyglider._config_components import (
    GliderDevices,
    Metadata,
    NetCDFVariables,
    ProfileVariables,
)

__all__ = ['Deployment', 'dump_yaml']


class Deployment(BaseModel):
    metadata: Metadata
    glider_devices: GliderDevices
    netcdf_variables: NetCDFVariables
    profile_variables: ProfileVariables

    @classmethod
    def load_yaml(cls, yaml_str: str) -> 'Deployment':
        """Load a yaml string into a Deployment model."""
        return _generic_load_yaml(yaml_str, cls)


def dump_yaml(model: BaseModel) -> str:
    """Dump a pydantic model to a yaml string."""
    return yaml.safe_dump(model.model_dump(), default_flow_style=False)


def _generic_load_yaml(data: str, model: BaseModel) -> BaseModel:
    """Load a yaml string into a pydantic model."""
    return model.model_validate(yaml.safe_load(data))
