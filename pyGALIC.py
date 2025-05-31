from random import uniform, gauss 

def find_closest(array, value):
    # find index of 'array' entry which is closest to 'value'
    diff = array - value
    idx = np.abs(diff).argmin()
    return idx

import h5py
import numpy as np

def write_IC_hdf5(filename, data):
 def ifset(d,key):
  if key in d.keys():
   result=d[key]
  else:
   result=None
   if key in ['lt','fmt']:  #for plot
    result=''
  return result

 def ifset2(d,key,value):
  if value is None:
   result=ifset(d,key)
  else:
   result=value
  return result

 if isinstance(data, dict):
  data=[data]

 BoxSize = None
 NumPartType = 6
 MassTable = np.zeros(NumPartType, dtype = float)
 NumPart = np.zeros(NumPartType, dtype = int)
 
 for d in data:
  BoxSize = ifset2(d, 'BoxSize', BoxSize)
  i = d['PartType']
  MassTable[i] = d['PartMass']
  NumPart[i] = d['count']
  
 file = h5py.File(filename+'.hdf5','w')

 group = file.create_group("Header")
 if BoxSize is not None:
  group.attrs.create("BoxSize", BoxSize, shape=None, dtype=h5py.h5t.IEEE_F64LE)
 else:
  group.attrs.create("BoxSize", 0, shape=None, dtype=h5py.h5t.IEEE_F64LE)
 group.attrs.create("Flag_Cooling", 0, shape=None, dtype=h5py.h5t.STD_I32LE)
 group.attrs.create("Flag_DoublePrecision", 0, shape=None, dtype=h5py.h5t.STD_I32LE)
 group.attrs.create("Flag_Feedback", 0, shape=None, dtype=h5py.h5t.STD_I32LE)
 group.attrs.create("Flag_IC_Info", 0, shape=None, dtype=h5py.h5t.STD_I32LE)
 group.attrs.create("Flag_Metals", 0, shape=None, dtype=h5py.h5t.STD_I32LE) #in makegal ics is 1
 group.attrs.create("Flag_Sfr", 0, shape=None, dtype=h5py.h5t.STD_I32LE)
 group.attrs.create("Flag_StellarAge", 0, shape=None, dtype=h5py.h5t.STD_I32LE)
 group.attrs.create("HubbleParam", 1, shape=None, dtype=h5py.h5t.IEEE_F64LE)
 group.attrs.create("MassTable", MassTable, shape=None,
                                                     dtype=h5py.h5t.IEEE_F64LE)
 group.attrs.create("NumFilesPerSnapshot", 1, shape=None, dtype=h5py.h5t.STD_I32LE)
 group.attrs.create("NumPart_ThisFile", NumPart, shape=None, dtype=h5py.h5t.STD_I32LE)
 group.attrs.create("NumPart_Total", NumPart, shape=None, dtype=h5py.h5t.STD_U32LE)
 group.attrs.create("NumPart_Total_HighWord", (0,0,0,0,0,0), shape=None, dtype=h5py.h5t.STD_U32LE)
 group.attrs.create("Omega0", 0, shape=None, dtype=h5py.h5t.IEEE_F64LE)
 group.attrs.create("OmegaLambda", 0, shape=None, dtype=h5py.h5t.IEEE_F64LE)
 group.attrs.create("Redshift", 0, shape=None, dtype=h5py.h5t.IEEE_F64LE)
 group.attrs.create("Time", 0, shape=None, dtype=h5py.h5t.IEEE_F64LE)

 ID_offset = 0
 for d in data:
  group = file.create_group("PartType"+str(d['PartType']))
  dataset = group.create_dataset("Coordinates", (d['count'],3), data=d['Coordinates'],
                                                            dtype=h5py.h5t.IEEE_F32LE)
  dataset = group.create_dataset("ParticleIDs", (d['count'],),
                   data=np.array(range(ID_offset,ID_offset+d['count'])), dtype=h5py.h5t.STD_I32BE)
  ID_offset += d['count']
  dataset = group.create_dataset("Velocities", (d['count'],3), data=d['Velocities'],
                                                            dtype=h5py.h5t.IEEE_F32LE)
  if d["PartType"] == 0:
   dataset = group.create_dataset("InternalEnergy", (d['count'],),
                   data=d["InternalEnergy"], dtype=h5py.h5t.IEEE_F32LE)
 file.close()

class IC:
    """
    Use to sample initial conditions for one particle species. 
    Usage: 
    IC(profile, N_part, r_max, species, args=(), G=43000, external_profile=lambda r:0, external_args=())
    Arguments:
        profile: function which returns the density at a given radius. Should be of the form profile(r, *args)
        N_part: number of particles to sample
        r_max: maximum radius to sample particles
        species: particle species (0: gas, 1: DM)
        args: extra arguments for the profile function
        G: gravitational constant in desired units. Default is 43000 (Gadget units)
        external_profile: function which returns the density of an external potential. 
                          E.g. the DM potential when sampling gas particles. 
                          Should be of the form external_profile(r, *external_args)
        external_args: extra arguments for the external_profile function
    
    Methods:
        sample(rmax, N_table=100): sample particles and return data 
            Arguments:
                rmax: maximum radius to sample particles
                N_table: number of points in the integration tables. Default is 100
                N_part: number of particles to sample. 
            Returns:
                sampled data in dictionary format. 
    """

    def __init__(self, profile, species, args=(), G=43000, external_profile=lambda r:0, external_args=()):
        self.profile = profile
        self.args = args
        self.G = G
        self.external_profile = external_profile
        self.external_args = external_args
        self.species = species

    def density(self, r):
        return self.profile(r, *self.args)
    
    def total_density(self, r):
        return self.profile(r, *self.args) + self.external_profile(r, *self.external_args)

    def encl_mass(self, r):
        return quad(lambda r: 4*np.pi * self.density(r) * r**2, 0, r)[0]
    
    def encl_mass_total(self, r):
        return quad(lambda r: 4*np.pi * self.total_density(r) * r**2, 0, r)[0]
    
    def temperature(self, r):
        # for collisionless particles, this is actually the squared velocity dispersion * 3/2
        p = quad(lambda r: self.G * self.density(r) * self.encl_mass_total(r)/r**2, r, np.inf)[0]
        return (3/2) * p / self.density(r)

    def sample(self, rmax, N_part, N_table=1000, mirrored=False):
        r_arr = np.logspace(-9+np.log10(rmax), np.log10(rmax), N_table)
        encl_mass_arr = np.array([self.encl_mass(r) for r in r_arr])
        CPD_arr = encl_mass_arr/np.amax(encl_mass_arr)
        temp_arr = np.array([self.temperature(r) for r in r_arr])

        if mirrored and N_part % 2 != 0:
            print('Warning: N_part must be an even number for mirrored sampling. Adding one particle.')
            N_part += 1
        N_sample = N_part if not mirrored else int(N_part/2)
        coords = np.zeros([N_part, 3])
        vels = np.zeros([N_part, 3])
        temps = np.zeros(N_part) if self.species == 0 else None

        for i in range(N_sample):
            # get particle radius by commulative invers sampling. I.e., sample a random number between 0 and 1
            # and solve M(r)/M(rmax) = x for r, where M(r) is the enclosed mass at r
            x = uniform(0, 1)
            r = r_arr[find_closest(CPD_arr, x)]
            # sample random angular position 
            phi = uniform(0, 1) * 2 * np.pi
            x = uniform(0.0,1.0)-0.5
            theta = np.arccos(-2.0*x)
            # set coordinates
            coords[i][0] = r*np.sin(theta)*np.cos(phi)
            coords[i][1] = r*np.sin(theta)*np.sin(phi)
            coords[i][2] = r*np.cos(theta)

            if self.species == 0:
                temps[i] = temp_arr[find_closest(r_arr, r)]
            else:
                sigma2 = temp_arr[find_closest(r_arr, r)] * 2/3 # velocity dispersion^2 is calculated as 2/3 gas temperature
                sigma = np.sqrt(sigma2)
                vels[i][0] = gauss(0, sigma)
                vels[i][1] = gauss(0, sigma)
                vels[i][2] = gauss(0, sigma)
        
        if mirrored:
            # find phase space coordinates of mirror particles
            coords[N_sample:] = -coords[:N_sample]
            vels[N_sample:] = -vels[:N_sample]
            if self.species == 0:
                temps[N_sample:] = temps[:N_sample]
        
        total_mass = self.encl_mass(rmax)
        data = {}
        data['count'] = N_part
        data['PartMass'] = total_mass/N_part
        data['PartType'] = self.species
        data['Coordinates'] = coords
        data['Velocities'] = vels
        if self.species == 0:
            data['InternalEnergy'] = temps 
        return data