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


def check_nwb_install():
    assert (
        HAVE_NWB
    ), "To use the Nwb extractors, install pynwb: \n\n pip install pynwb\n\n"


def update_dict(d, u):
    for k, v in u.items():
        if isinstance(v, abc.Mapping):
            d[k] = update_dict(d.get(k, {}), v)
        else:
            d[k] = v
    return d


def set_dynamic_table_property(
    dynamic_table,
    row_ids,
    property_name,
    values,
    index=False,
    default_value=np.nan,
    description="no description",
):
    check_nwb_install()
    if not isinstance(row_ids, list) or not all(isinstance(x, int) for x in row_ids):
        raise TypeError("'ids' must be a list of integers")
    ids = list(dynamic_table.id[:])
    if any([i not in ids for i in row_ids]):
        raise ValueError("'ids' contains values outside the range of existing ids")
    if not isinstance(property_name, str):
        raise TypeError("'property_name' must be a string")
    if len(row_ids) != len(values) and index is False:
        raise ValueError("'ids' and 'values' should be lists of same size")

    if index is False:
        if property_name in dynamic_table:
            for (row_id, value) in zip(row_ids, values):
                dynamic_table[property_name].data[ids.index(row_id)] = value
        else:
            col_data = [default_value] * len(ids)  # init with default val
            for (row_id, value) in zip(row_ids, values):
                col_data[ids.index(row_id)] = value
            dynamic_table.add_column(
                name=property_name, description=description, data=col_data, index=index
            )
    else:
        if property_name in dynamic_table:
            # TODO
            raise NotImplementedError
        else:
            dynamic_table.add_column(
                name=property_name, description=description, data=values, index=index
            )


def add_nwb_device(recording, nwbfile):
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

    # only ephys devices are recorded in NWB
    nwb_devices = recording.nwb_metadata["Ecephys"]["Device"]
    # extract ephys device information from bonsai metadata
    for dev in recording.metadata["Device"]:
        # TODO: handle cases with multiple ephys devices?
        if dev["type"] == "ephys":
            nwbfile.create_device(name=dev["name"])
            nwb_devices.append({"name": dev["name"]})

    # create new device if no ephys devices found in bonsai metadata
    if len(nwb_devices) == 0:
        nwbfile.create_device(name="Device")
        nwb_devices.append({"name": "Device"})

    recording.nwb_metadata["Ecephys"]["Device"] = nwb_devices
    return nwbfile


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
        Auxiliary static method for nwbextractor.
        Adds channels from recording object as electrodes to nwbfile object.
        """
    # Check for existing electrodes
    if nwbfile.electrodes is not None:
        nwb_elec_ids = nwbfile.electrodes.id.data[:]
    else:
        nwb_elec_ids = []

    # Extractors channel groups must be integers, but Nwb electrodes group_name can be strings
    # nwb_groups_names = [str(grp['name']) for grp in metadata['Ecephys']['ElectrodeGroup']]
    nwb_groups_names = list(nwbfile.electrode_groups.keys())

    # For older versions of pynwb, we need to manually add these columns
    if distutils.version.LooseVersion(pynwb.__version__) < "1.3.0":
        if nwbfile.electrodes is None or "rel_x" not in nwbfile.electrodes.colnames:
            nwbfile.add_electrode_column(
                "rel_x", "x position of electrode in electrode group"
            )
        if nwbfile.electrodes is None or "rel_y" not in nwbfile.electrodes.colnames:
            nwbfile.add_electrode_column(
                "rel_y", "y position of electrode in electrode group"
            )

    # add new electrodes with id, (rel_x, rel_y) and groups
    channel_ids = list(recording.get_channel_ids())
    for m in channel_ids:
        if m not in nwb_elec_ids:
            location = recording.get_channel_locations(channel_ids=m)[0]
            grp_name = recording.get_channel_groups(channel_ids=m)[0]
            grp = nwbfile.electrode_groups[nwb_groups_names[grp_name]]
            impedance = -1.0
            nwbfile.add_electrode(
                id=m,
                x=np.nan,
                y=np.nan,
                z=np.nan,
                rel_x=float(location[0]),
                rel_y=float(location[1]),
                imp=impedance,
                location="unknown",
                filtering="none",
                group=grp,
            )
    electrode_table = nwbfile.electrodes

    # add/update electrode properties
    for ch in channel_ids:
        rx_channel_properties = recording.get_channel_property_names(channel_id=ch)
        for pr in rx_channel_properties:
            val = recording.get_channel_property(ch, pr)
            desc = "no description"
            # property 'location' of RX channels corresponds to rel_x and rel_ y of NWB electrodes
            if pr == "location":
                names = ["rel_x", "rel_y"]
                for (nm, v) in zip(names, val):
                    set_dynamic_table_property(
                        dynamic_table=electrode_table,
                        row_ids=[ch],
                        property_name=nm,
                        values=[float(v)],
                        default_value=np.nan,
                        description=nm + " coordinate location on the implant",
                    )
                continue
            # property 'group' of electrodes can not be updated
            if pr == "group":
                continue
            # property 'gain' should not be in the NWB electrodes_table
            if pr == "gain":
                continue
            # property 'brain_area' of RX channels corresponds to 'location' of NWB electrodes
            if pr == "brain_area":
                pr = "location"
                desc = "brain area location"
            set_dynamic_table_property(
                dynamic_table=electrode_table,
                row_ids=[ch],
                property_name=pr,
                values=[val],
                default_value=np.nan,
                description=desc,
            )

    return nwbfile


def create_nwb(
    recording, save_path,
):
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
        nwbfile = add_nwb_device(recording=recording, nwbfile=nwbfile,)

        # Add electrode groups
        nwbfile = add_nwb_electrode_groups(recording=recording, nwbfile=nwbfile,)

        # Add electrodes
        nwbfile = add_nwb_electrodes(recording=recording, nwbfile=nwbfile,)

        # Write to file
        io.write(nwbfile)
