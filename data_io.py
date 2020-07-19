import numpy as np
import h5py
import yt
from fields_defination import \
        setup_fluid_fields
class ScorpioLoader():
    def __init__(self,flnm,variables,dimensions,periodicity):
        """
        Input:
            flnm: file name of the h5 data to be load
            variables: an list with length of 8, describing which fields are included
            dimensions: the initial cell number on each direction.
            periodicity: boolean array with length of 3. describing which direction has period boundary condition
        """
        self.flnm = flnm
        self.variables=variables
        self.periodicity = periodicity
        self.ndim = len(dimensions)
        self.dimensions=np.ones(3,dtype=int)
        self.dimensions[:self.ndim]=dimensions
        self.file_opened=False
        self.file_closed=True
        # the variable index <-> name table
        self.index_name_table={
            0 :'density', 7 :'ene',
            1 :'mom1',    2 :'mom2', 3  :'mom3',
            4 :'bl1',     5 :'bl2',  6  :'bl3',
            8 :'br1',     9 :'br2',  10 :'br3'}
        revd=dict([reversed(i) for i in self.index_name_table.items()])
        self.index_name_table.update(revd)

        code_density  = "code_mass / code_length**3"
        code_momentum = "code_mass / code_length**2 / code_time"
        code_pressure = "code_mass / code_length / code_time**2"
        self.name_unit_table={
                'density': code_density    , 'ene' : code_pressure   ,
                'mom1'   : code_momentum   , 'mom2': code_momentum   , 'mom3' : code_momentum   ,
                'bl1'    : "code_magnetic" , 'bl2' : "code_magnetic" , 'bl3'  : "code_magnetic" ,
                'br1'    : "code_magnetic" , 'br2' : "code_magnetic" , 'br3'  : "code_magnetic" }
    def open_file(self):
        if not self.file_opened:
            self.file = h5py.File(self.flnm,'r')
            self.file_opened = True
            self.file_closed = False
        else:
            pass
        return 

    def close_file(self):
        if self.file_opened and not self.file_closed:
            self.file.close()
            self.file_closed=True
            self.file_opened=False
        else:
            pass
        return 

    def _preload(self):
        """
        preload maximum level, number of buffer, number of variables etc.
        """
        self.lmax = int(self.file['lmax'][()])
        self.block_dimensions = np.ones(3,dtype=np.int)
        self.block_dimensions[:self.ndim] = self.file['blockSize'][()].astype(int)
        self.nbuf = int(self.file['nbuf'][()])
        self.nvar = int(self.file['nvar'][()])
        num_var = (np.asarray(self.variables)*np.array([1,1,1,1,2,2,2,1])).sum()
        if num_var != self.nvar:
            sys.exit("Error: 'variables' do not match nvar")
        self.time = self.file['t'][()]
        self.lv_blockID=[]
        for i in range(0,self.lmax+1):
            dsetname = 'lv{:02d}_gridID'.format(i)
            self.lv_blockID.append(self.file[dsetname][()])
        self.domain = [] 
        for dim in range(self.ndim):
            self.domain.append(range(self.nbuf,self.nbuf+self.block_dimensions[dim]))
        while len(self.domain) < 3:
            self.domain.append([0])
        _buf_list = np.zeros(3,dtype='int32')
        _buf_list[:self.ndim] = 2*self.nbuf
        self.block_data_shape = tuple(([i+j for i,j in zip(self.block_dimensions,_buf_list)]+[self.nvar])[::-1])
        return

    def _read_one_block_coord(self,level,blockID):
        """
        Get the left and right edge of a block
        """
        one_block = {}
        xl=np.zeros(3)
        xr=np.ones(3)
        for dim in range(0,self.ndim):
            dsetname = 'lv{:02d}_{:010d}_xc{:d}'.format(level,blockID,dim+1)
            xc = self.file[dsetname][self.domain[dim]]
            dxl = xc[1]-xc[0]
            dxr = xc[-1]-xc[-2]
            xl[dim] = (xc[0] -dxl/2)
            xr[dim] = (xc[-1]+dxr/2)
        one_block['left_edge']  = np.array(xl)
        one_block['right_edge'] = np.array(xr)
        return one_block
        
    def _read_one_block_variable(self,level,blockID,variable_indese):
        """
        Get the data of one block
        """
        one_block = {}
        dsetname = 'lv{:02d}_{:010d}_q1'.format(level,blockID)
        def _grid1d(variable_indese):
            grid   = np.ix_(variable_indese,*self.domain[::-1]) 
            grid1d = np.ravel_multi_index(grid,self.block_data_shape)
            return grid1d
        # This indexing method only valid for numpy array, 
        # thus need to convert hdf5 data to numpy array anyway.
        # a potential improvment is adding reference when output data in scorpio
        variables = self.file[dsetname][()][_grid1d(variable_indese)]
        for variable_index in variable_indese:
            variable_name=self.index_name_table[variable_index]
            unit = self.name_unit_table[variable_name]
            one_block[variable_name] = (variables[variable_index].T, unit)
        return one_block

    def load(self,request_fields=None):
        """
        Load all blocks' information, pack into a dictionary and feed to YT.
        """
        block_data = []
        if request_fields is None:
            request_fields = np.zeros(11)
            request_fields[:8] = self.variables
            request_fields[8:] = request_fields[4:7]
        variable_indese = np.flatnonzero(request_fields)
        self.open_file()
        self._preload()
        xl_domain = [0.]*3
        xr_domain = [0.]*3
        for level in range(0,self.lmax+1):
            for blockID in self.lv_blockID[level].astype(np.int):
                one_block = self._read_one_block_coord(level,blockID)
                one_block.update(self._read_one_block_variable(level,blockID,variable_indese))
                one_block['level']      = level+1
                one_block['dimensions'] = self.block_dimensions
                block_data.append(one_block)
                xl_domain = np.minimum(xl_domain,one_block['left_edge'])
                xr_domain = np.maximum(xr_domain,one_block['right_edge'])
        self.close_file()
        bbox = np.vstack([xl_domain,xr_domain]).T.astype(np.float32)
        ds = yt.load_amr_grids(block_data, self.dimensions, bbox,
                sim_time = self.time,
                unit_system='code',
                periodicity=self.periodicity)
        setup_fluid_fields(ds) 
        return ds
