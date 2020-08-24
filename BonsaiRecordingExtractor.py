from spikeextractors import RecordingExtractor, BinDatRecordingExtractor
from bs4 import BeautifulSoup
from pathlib import Path
import dateutil.parser
import lxml
import csv
import re



def find_soup(soup, **kwargs):
    """
    Raise error when no matching text is found in soup

    Parameters
    ----------
    soup : BeautifulSoup object 
    kwargs : parameters to pass into find
    """
    if soup.find(**kwargs) is None: 
        raise ValueError(f"Cannot find text matching {kwargs} in file. Check if specified file path is correct")
    return soup.find(**kwargs)


class BonsaiRecordingExtractor(BinDatRecordingExtractor):
    extractor_name = 'BonsaiRecording'
    has_default_locations = False
    is_writable = False
    mode = 'folder'
    installation_mesg = "For parsing metadata from Bonsai XML outputs"

    """
    Extract metadata from Bonsai XML and csv ouputs

    """
    # TODO

    def __init__(self, bonsai_dir, metadata_file, traces_file, 
                 sampling_frequency=None, numchan=None, dtype=None,
                 time=None, time_file=None, **kwargs):
        """

        Parameters
        ----------
        bonsai_dir: str
            Path of directory that contains all outputs from Bonsai.  
            Files shouldn't be nested
        metatdata_file : str
            Path of bonsai XML file that contains workflow metadata relative to 
            `bonsai_dir`.  Usually ends with ".bonsai".  
        traces_file : str
            Path of binary time series files to read traces from TODO?
        sampling_frequency : int or None, optional
            Directory containing all Bonsai output
        numchan: int, optional
            Number of channels
        dtype : str
            TODO
        """
        self._bonsai_dir = str(Path(bonsai_dir).absolute())
        self._metadata_file = str(Path(self._bonsai_dir) / metadata_file)

        with open(self._metadata_file, 'r') as f:
            soup = BeautifulSoup(f, 'lxml')
        # remove disabled steps
        for tag in soup.find_all('expression', {'xsi:type': "Disable"}):
            tag.decompose()

        self.nodes = find_soup(soup, name='nodes')  # workflow information
        self.edges = find_soup(soup, name='edges')  # how channels connect to each other

        # attributes set at init
        if sampling_frequency is None:
            self.sampling_frequency = self.get_sampling_frequency()
        if numchan is None:
            self.numchan = self.get_num_channels()
        if dtype is None:
            self.dtype = self.get_bin_dat_dtype()
        if self.dtype not in ['uint16', 'float32']:
            raise ValueError(f"dtype {self.dtype} not valid. Choose 'uint16' or 'float32'")
        if time is not None or time_file is not None:
            self.add_session_start_time(time, time_file)
    
        traces_file = str(Path(self._bonsai_dir) / traces_file)
        super().__init__(traces_file, self.sampling_frequency, self.numchan, self.dtype, **kwargs)

        # compile metadata
        self.create_metadata()
        #self.metadata['devices'] = self.add_devices()


    def add_session_start_time(self, time=None, file=None):
        """ 
        Specifc session start time or return first timestamp found in a csv file.
        Needs to specific either `time` or `file`.  

        Parameters
        ----------
        time : datetime object or str that can be converted to datetime object
        file : csv file path  
        """
        if time is None and file is None:
            raise Exception('Specify either `time` or `time_file`')

        dt = None
        if time:
            try: # strings
                dt = dateutil.parser.parse(time).isoformat()
            except parser.ParserError:
                return time.isoformat()
            except Exception:
                raise

        if file and dt is not None:
            fp = str(Path(self._bonsai_dir) / file)
            with open(fp) as f:
                reader = csv.reader(f)
                for row in reader:
                    for col in row:
                        try: # returns first time stamp found
                            dt = dateutil.parser.parse(col).isoformat() 
                            break
                        except Exception:
                            continue
                raise ValueError(f'Timestamp not found in file: {file}')

        self.session_start_time = dt


    def extract_device_metadata(self, ephys_device ='rhd', file=None):
        """ 
        Add device information from Bonsai XML metadata 


        Parameters
        ----------
        ephys_device: str, optional
        file: str, optional
            Path of XML file to grab metadata from relative to self._bonsai_dir
            If None, defaults to self._metadata_file
        """
        # clean file
        if file is None: 
            soup = self.nodes  # grab all
        else:
            with open(str(Path(self._bonsai_dir) / file), 'r') as f:
                soup = BeautifulSoup(f, 'lxml')
            # remove disabled steps
            for tag in soup.find_all('expression', {'xsi:type': "Disable"}):
                tag.decompose()

        # ecephys device only (be consistent with NWB extractor metadata)
        self.metadata['Ecephys']['Device'] = []  
        # all devices, including ephys device. 
        self.metadata['Device'] = []


        for d in soup.find_all(re.compile('.*:deviceindex')):
            device_md = dict()
            
            # device name is in the parent tag  e.g. <Combinator xsi:type="q1:InfoDevice">
            device_md['name'] = list(d.parent.attrs.values())[0].split(':')[-1]  
            d_id = find_soup(d,name=re.compile('.*:selectedindex'))
            device_md['id'] = int(d_id.get_text())

            #TODO raise Exception('device already in index')    
            # loop over attribute information in a device
            for s in d.next_siblings:
                if s is not None and s.name:
                    attr_name = d.next_sibling.name.split(':')[-1]
                    text = d.next_sibling.text
                    # convert strings to bool when necessary
                    if text == 'true':
                        text = True
                    elif text == 'false':
                        text = False
                    device_md[attr_name] = text
                d = d.next_sibling

            if re.search(ephys_device, device_md['name'], re.IGNORECASE) is not None:
                self.metadata['Ecephys']['Device'].append(device_md) 
                device_md['ephys_device'] = True
            else:
                device_md['ephys_device'] = False
            
            self.metadata['Device'].append(device_md)


    def get_all_valid_files(self):
        # filter out empty files, grab metadata for each file
        # distinguish between inputs and outputs?

        # binary files
        bin_md = {}
        for s in soup.find_all('combinator', {'xsi:type': re.compile('.*matrixwriter', re.IGNORECASE)}):
            file_md = {}
            for attr in s.children:
                if attr.name is not None:
                    file_md[attr.name.split(':')[-1]] = attr.string
            if 'path' in file_md:
                if 'suffix' in file_md:
                    file_md['prefix'], file_md['ext'] = file_md['path'].rsplit('.')
            
            for f in all_files:
                if f.startswith(file_md['prefix']) and f.endswith(file_md['ext']):
                    bin_md[f] = file_md
        return(bin_md)


    # TODO - not hard code?
    def create_metadata(self, **kwargs): 
        self.metadata = dict()
        self.metadata['Ecephys'] = {
            'name': 'ephys_data',
            'data': self.get_traces(), # TODO add args
             #electrodes=electrode_table_region,
            #'starting_time': self.add_session_start_time(),
            'sampling_frequency': self.get_sampling_frequency(),
            'comments': 'Generated from SpikeInterface::BonsaiRecordingExtractor',
            'description': 'acquisition_description',
         }
        # add self.metadata['Ecephys']['Device'] and self.metadata['Devices']
        self.extract_device_metadata()
        

    # TODO: make file or device specific?
    def get_bin_dat_dtype(self): 
        """ 
        Determine data type of binary data that contains time series.

        Default is `uint16`.  For any binary data saved after an `AdcScale` operator, which 
        applies the 0.195 scale factor for spike detection, dtype is `float32`

        """
        pat = re.compile(r'.*adcscale', re.IGNORECASE)
        try:
            adc = find_soup(self.nodes, name='combinator', attrs={'xsi:type': pat})
        except ValueError:  # no scaling in workflow
            return('uint16')
        return('float32')



    # TODO: check all instances of channel ?  
    # QUESTION: how to figure out channel count when it's not speciftied
    def get_num_channels(self, device=None): 
        """ 
        Return number of channels in int 

        Notes
        --------
        This assumes channel counts are in tag 'dsp:channelcount'
        """
        text = find_soup(self.nodes, name=re.compile(r'.*channelcount$'))
        numchan = int(text.get_text())
        return numchan


    # TODO: get device names from <q1:DeviceIndex> and id # from <q1:SelectedIndex>
    def get_channel_ids(self):
        """ Return list of channel ids. If channel ids are not specified, return range(num_channels) """
        channel_ids = range(self.get_num_channels())
        return channel_ids


    def get_sampling_frequency(self, device='RHD', index=None):
        """ 
        Return sampling rate of recordings in Hertz 

        Parameters
        ----------
        device : str, optional
            Name of device. Case insensitive, supports regular expression.

        """
        pat = re.compile(r'.*' + device + '.*', re.IGNORECASE)
        device_info = self.nodes.find_all('combinator', {'xsi:type': pat})
        if len(device_info) != 1:
            raise Exception('More than one RHD device, revise code')

        sf = find_soup(device_info[0], name=re.compile(r'.*samplerate$')).get_text()
        sf = re.findall(r'\d+', sf)[0]  # samplerate can be a string
        return float(sf)  


