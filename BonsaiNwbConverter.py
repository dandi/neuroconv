"""
Functions modified from NWBRecordingExtractor from spikeextractors

Example
-------
from BonsaiRecordingExtractor import BonsaiRecordingExtractor
from BonsaiNwbConverter import *

bonsai_params = {'bonsai_dir': 'data/',
                'metadata_file': 'aquisition_oni.bonsai',
                'traces_file': 'intan.raw',
                'time':'intan-first-time.csv'}

rx = BonsaiRecordingExtractor(**bonsai_params) 
create_nwb(rx, 'bonsai.nwb')

"""

from BonsaiRecordingExtractor import BonsaiRecordingExtractor
from spikeextractors import NwbRecordingExtractor
from spikeextractors.extractors.nwbextractors import *

import os
import uuid
from datetime import datetime
from collections import defaultdict, abc
from pathlib import Path
import numpy as np
import distutils.version

try:
    import pynwb
    from pynwb import NWBHDF5IO
    from pynwb import NWBFile
    from pynwb.ecephys import ElectricalSeries
    from pynwb.ecephys import ElectrodeGroup

    HAVE_NWB = True
except ModuleNotFoundError:
    HAVE_NWB = False


def add_nwb_devices(recording, nwbfile):
    """ 
    Add relevant device info from Bonsai recording to NWB file & 
    updates NWB metadata in recording  

    Note: only adds ephys devices at the moment

    Parameters
    ----------
    recording: RecordingExtractor
    save_path: NWB object
    """

    if "Ecephys" not in recording.nwb_metadata:
        recording.nwb_metadata["Ecephys"] = dict()

    if "Device" not in recording.nwb_metadata["Ecephys"]:
        recording.nwb_metadata["Ecephys"]["Device"] = []

    # TODO: handle cases with multiple ephys devices?
    for dev in recording.metadata["Device"]:
        # only ephys devices are recorded in NWB
        if dev["type"] == "ephys":
            # nwbfile.create_device(name=dev["name"])
            if "Device" not in recording.nwb_metadata["Ecephys"]:
                recording.nwb_metadata["Ecephys"]["Device"] = []
            # save metadata to recording
            recording.nwb_metadata["Ecephys"]["Device"].append({"name": dev["name"]})

    return NwbRecordingExtractor.add_devices(
        recording=recording, nwbfile=nwbfile, metadata=recording.nwb_metadata
    )


def add_nwb_electrode_groups(recording, nwbfile):
    """ 
    Adds relevant device info from Bonsai recording to NWB file & 
    updates NWB metadata in recording  

    Note: assumes electrode group information not available in recording.metadata

    Parameters
    ----------
    recording: RecordingExtractor
    save_path: NWB object
    """
    metadata = recording.nwb_metadata["Ecephys"]
    channel_ids = recording.get_channel_ids()

    # Electrode groups
    ephys_metadata = recording.nwb_metadata["Ecephys"]
    if "ElectrodeGroup" not in ephys_metadata:
        ephys_metadata["ElectrodeGroup"] = []

        # Check if 'groups' property exists in self._channel_properties
        if "group" in recording.get_shared_channel_property_names():
            RX_groups_names = list(np.unique(recording.get_channel_groups()))
        else:
            RX_groups_names = ["0"]
            # Electrode groups are required for NWB, for consistency we create group for Recording channels
            vals = [0] * len(channel_ids)
            recording.set_channel_groups(channel_ids=channel_ids, groups=vals)
        for grp_name in RX_groups_names:
            ephys_metadata["ElectrodeGroup"].append(
                {
                    "name": grp_name,
                    "description": "electrode_group_description",
                    "location": "electrode_group_location",
                    "device": ephys_metadata["Device"][0]["name"],
                }
            )

    # Tests if electrode groups exist in nwbfile, if not create them from metadata
    for grp in ephys_metadata["ElectrodeGroup"]:
        if str(grp["name"]) not in nwbfile.electrode_groups:
            nwbfile.create_electrode_group(
                name=str(grp["name"]),
                location=grp["location"],
                device=nwbfile.devices[grp["device"]],
                description=grp["description"],
            )

    recording.nwb_metadata["Ecephys"] = ephys_metadata
    return nwbfile


def add_nwb_electrodes(recording, nwbfile):
    """
    Note: NwbRecordingExtractor.add_electrodes doesn't add information to NWB metadata
    """
    return NwbRecordingExtractor.add_electrodes(
        recording=recording, nwbfile=nwbfile, metadata=None
    )


def create_nwb(recording, save_path):
    """
    Use metadata in BonsaiRecordingExtractor to create a NWB files

    Parameters
    ----------
    recording: BonsaiRecordingExtractor
    save_path: str
    nwb_metadata: dict
        extra metadata info for constructing the nwb file (optional).
    """
    assert HAVE_NWB, NwbRecordingExtractor.installation_mesg

    assert (
        distutils.version.LooseVersion(pynwb.__version__) >= "1.3.3"
    ), "'write_recording' not supported for version < 1.3.3. Run pip install --upgrade pynwb"

    if os.path.exists(save_path):
        read_mode = "r+"
    else:
        read_mode = "w"

    # Update any previous metadata with user passed dictionary
    if recording.nwb_metadata is None:
        recording.nwb_metadata = dict()

    with NWBHDF5IO(save_path, mode=read_mode) as io:
        if read_mode == "r+":
            nwbfile = io.read()
        else:
            if "NWBFile" not in recording.nwb_metadata:
                recording.nwb_metadata["NWBFile"] = {
                    "session_description": "no description",
                    "identifier": str(uuid.uuid4()),
                    "session_start_time": recording.session_start_time,
                }
            nwbfile = NWBFile(**recording.nwb_metadata["NWBFile"])

        # Add devices
        nwbfile = add_nwb_devices(recording=recording, nwbfile=nwbfile,)

        # Add electrode groups
        nwbfile = add_nwb_electrode_groups(recording=recording, nwbfile=nwbfile,)

        # Add electrodes
        nwbfile = add_nwb_electrodes(recording=recording, nwbfile=nwbfile,)

        # Add electrical series
        # nwbfile = add_nwb_electrical_series(recording=recording, nwbfile=nwbfile)

        # Write to file
        io.write(nwbfile)
