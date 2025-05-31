# Initial Condition Generator with `pyGALIC`

## 📦 Dependencies

This code relies on the following Python packages:

- `numpy`
- `scipy`
- `h5py`

---

## 🔧 Usage

### 1. Define a Density Profile

To sample ICs for a certain particle species, first define the desired density distribution. For example, for an **NFW** profile:

```python
def NFW(r, r_s, rho_s):
    x = r / r_s
    return rho_s / (x * (1 + x)**2)
```

> The first argument must be the radius in your desired units.  
> The profile must be non-zero everywhere, and ensure that ρ(r)·r² does not diverge within the sampling radius.

---

### 2. Initialize the IC Generator

```python
from pyGALIC import IC
ic_DM = IC(NFW, species=1, args=(10, 1), G=43000)
```

- `species=1`: DM particles (per Gadget-3 convention)
- `species=0`: Gas particles
- `args`: Tuple of arguments for the profile function (after `r`)
- `G`: Gravitational constant (default is Gadget-3 units)

---

### 3. Sample the Initial Conditions

```python
data = ic_DM.sample(rmax=100, N_part=int(1e7), N_table=1000, mirrored=False)
```

- `rmax`: Maximum sampling radius
- `N_part`: Number of particles
- `N_table`: Number of entries in the lookup table
- `mirrored`: Whether to sample mirrored particle pairs (makes total momentum exactly zero)

---

### 4. Save Data in HDF5 Format

```python
from pyGALIC import write_IC_hdf5
filename = 'halo'
write_IC_hdf5(filename, [data])
```

---

## 🌌 Multiple Species in Equilibrium

### 1. Define Density Profiles

```python
# Gas: Hernquist profile
rho_s_gas = 3.5e-2
r_s_gas = 0.77
def density_gas(r):
    x = r / r_s_gas
    return rho_s_gas / (x * (1 + x)**3)

# DM: NFW profile
rho_s_DM = 6.9e-4
r_s_DM = 9.1
def density_DM(r):
    x = r / r_s_DM
    return rho_s_DM / (x * (1 + x)**2)
```

---

### 2. Initialize ICs with External Potentials

```python
from pyGALIC import IC

ic_DM = IC(density_DM, species=1, external_profile=density_gas)
ic_gas = IC(density_gas, species=0, external_profile=density_DM)
```

> `external_profile` accounts for an external potential affecting the Jeans equation but does **not** contribute sampled particles.  
> If it requires arguments beyond radius, pass them with `external_args`.

---

### 3. Sample and Save

```python
from pyGALIC import write_IC_hdf5

N_part_DM = int(1e7)
N_part_gas = int(1e5)
r_max_DM = 10 * r_s_DM
r_max_gas = 10 * r_s_gas

data_DM = ic_DM.sample(N_part=N_part_DM, rmax=r_max_DM)
data_gas = ic_gas.sample(N_part=N_part_gas, rmax=r_max_gas)

write_IC_hdf5('halo', [data_DM, data_gas])
```

---
