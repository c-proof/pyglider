from typing import Optional

from pydantic import BaseModel, HttpUrl


class Metadata(BaseModel):
    acknowledgement: str
    comment: str
    contributor_name: str
    contributor_role: str
    creator_email: str
    creator_name: str
    creator_url: HttpUrl
    deployment_id: str
    deployment_name: str
    deployment_start: str
    deployment_end: str
    format_version: str
    glider_name: str
    glider_serial: str
    glider_model: str
    glider_instrument_name: str
    glider_wmo: str
    institution: str
    keywords: str
    keywords_vocabulary: str
    license: str
    metadata_link: HttpUrl
    Metadata_Conventions: str
    naming_authority: str
    platform_type: str
    processing_level: str
    project: str
    project_url: HttpUrl
    publisher_email: str
    publisher_name: str
    publisher_url: HttpUrl
    references: str
    sea_name: str
    source: str
    standard_name_vocabulary: str
    summary: str
    transmission_system: str
    wmo_id: str


class Device(BaseModel):
    make: str
    model: str
    serial: str
    long_name: Optional[str] = None
    make_model: Optional[str] = None
    factory_calibrated: Optional[str] = None
    calibration_date: Optional[str] = None
    calibration_report: Optional[str] = None
    comment: Optional[str] = None


class GliderDevices(BaseModel):
    pressure: Device
    ctd: Device
    optics: Device
    oxygen: Device


class NetCDFVariable(BaseModel):
    source: str
    long_name: Optional[str] = None
    standard_name: Optional[str] = None
    units: Optional[str] = None
    axis: Optional[str] = None
    coordinates: Optional[str] = None
    conversion: Optional[str] = None
    comment: Optional[str] = None
    observation_type: Optional[str] = None
    platform: Optional[str] = None
    reference: Optional[str] = None
    valid_max: Optional[float] = None
    valid_min: Optional[float] = None
    coordinate_reference_frame: Optional[str] = None
    instrument: Optional[str] = None
    accuracy: Optional[float] = None
    precision: Optional[float] = None
    resolution: Optional[float] = None
    positive: Optional[str] = None
    reference_datum: Optional[str] = None
    coarsen: Optional[int] = None


class NetCDFVariables(BaseModel):
    timebase: Optional[NetCDFVariable] = (
        None  #! Is this required? `example-slocum`` doesn't have it
    )
    time: NetCDFVariable
    latitude: NetCDFVariable
    longitude: NetCDFVariable
    heading: NetCDFVariable
    pitch: NetCDFVariable
    roll: NetCDFVariable
    conductivity: NetCDFVariable
    temperature: NetCDFVariable
    pressure: NetCDFVariable
    chlorophyll: NetCDFVariable
    cdom: NetCDFVariable
    backscatter_700: NetCDFVariable
    oxygen_concentration: NetCDFVariable
    temperature_oxygen: Optional[NetCDFVariable] = (
        None  #! Is this required? `example-slocum`` doesn't have it
    )


class ProfileVariable(BaseModel):
    comment: str
    long_name: str
    valid_max: Optional[float] = None
    valid_min: Optional[float] = None
    observation_type: Optional[str] = None
    platform: Optional[str] = None
    standard_name: Optional[str] = None
    units: Optional[str] = None
    calendar: Optional[str] = None
    type: Optional[str] = None
    calibration_date: Optional[str] = None
    calibration_report: Optional[str] = None
    factory_calibrated: Optional[str] = None
    make_model: Optional[str] = None
    serial_number: Optional[str] = None


class ProfileVariables(BaseModel):
    profile_id: ProfileVariable
    profile_time: ProfileVariable
    profile_time_start: ProfileVariable
    profile_time_end: ProfileVariable
    profile_lat: ProfileVariable
    profile_lon: ProfileVariable
    u: ProfileVariable
    v: ProfileVariable
    lon_uv: ProfileVariable
    lat_uv: ProfileVariable
    time_uv: ProfileVariable
    instrument_ctd: ProfileVariable
