import numpy as np
from scipy.io import FortranFile
from iris import save
from iris.coords import DimCoord
from iris.cube import Cube, CubeList
from glob import glob
import argparse


# Create Iris cube from given data and name with the given latitudes
def create_cube(data, name, lats):
    ntime, _, nx = data.shape
    xs = np.linspace(0, float(nx), nx, dtype=np.float32)
    time = DimCoord(np.arange(ntime, dtype=np.float32), standard_name='time', var_name='time',\
        units=1)
    x = DimCoord(xs, standard_name='longitude', long_name='longitude', var_name='x')
    y = DimCoord(lats, standard_name='latitude', long_name='latitude', var_name='y')

    return Cube(data, dim_coords_and_dims=[(time,0), (y, 1), (x, 2)], long_name=name, var_name=name)


# Process reduced Gaussian grid data
def proc_red_grid(data, nlats, lons, lonmax):
    data_procd = np.full((nlats,lonmax), np.nan)
    j = 0
    # Iterate over each latitude, putting the raw data from that latitude in
    # the centre of the container array
    for i, nlon in enumerate(lons):
        if nlon < lonmax:
            data_procd[i,int((lonmax-nlon)/2):int((lonmax+nlon)/2)] = data[j:j+nlon]
        else:
            data_procd[i,:] = data[j:j+nlon]
        j += nlon

    return data_procd


# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("field", help="The field to process (S/U/V/T)", type=str)
parser.add_argument("gridtype", help="The type of grid (O/F/H)", type=str)
parser.add_argument("ngauss", help="The guassian number", type=int)
parser.add_argument("prec",help="The precision (sp/dp)",type=str)
parser.add_argument("endian",help="The endian (big/little)",type=str)
args = parser.parse_args()
field = args.field
gridtype = args.gridtype
nlats = 2*(args.ngauss)
prec=args.prec
if args.endian=="big":
    endian=">"
elif args.endian=="little":
    endian="<"

# Calculate maximum number of longitudes in a latitude band (i.e. at the
# Equator)
if gridtype == "O":
    lonmax = 20 + 2*nlats-4
    # Define number of longitudes at each latitude band for the cubic octahedral reduced grid
    northern_hemisphere_lats = [20 + 4*i for i in range(nlats//2)]
    lons = northern_hemisphere_lats + northern_hemisphere_lats[::-1]
elif gridtype == "F":
    lonmax = 2*nlats
    # Define number of longitudes at each latitude band for the cubic full grid
    northern_hemisphere_lats = [lonmax for i in range(nlats//2-1)]
    lons = northern_hemisphere_lats + northern_hemisphere_lats[::-1]
elif gridtype == "H":
    lonmax = nlats
    # Define number of longitudes at each latitude band for the cubic HEALPix reduced grid
    northern_polar_lats = [4+4*i for i in range(nlats//4)]
    northern_tropic_lats = [nlats for i in range(nlats//4)] 
    northern_hemisphere_lats = northern_polar_lats + northern_tropic_lats
    lons = northern_hemisphere_lats + northern_hemisphere_lats[::-1]
else:
    raise ValueError(f"Unsuported gridtype ${gridtype}. It must one of O/F/H.")


# Define regular latitudes (I'm too lazy to do it properly for a Gaussian grid)
lats = np.linspace(90, -90, nlats)

# Process output data for given rundir
time_slices = []
files = sorted(glob(f"./{field}.???.0001.dat"))
if len(files) == 0:
    raise ValueError(f"Current directory does not seem to contain any {field}.???.????.dat files")

for filename in files:
    # Open binary file
    print(f"Opening {filename}")
    if prec == "sp":
        f = FortranFile(filename,'r',endian+'u4')
        data = f.read_record(dtype=endian+'f4')
        f.close()
    elif prec == "dp":
        f = FortranFile(filename,'r',endian+'u4')
        data = f.read_record(dtype=endian+'f8')
        f.close()

    # Process data
    data_procd = proc_red_grid(data, nlats, lons, lonmax)

    # Parse 1D array containing reduced Gaussian grid into lon-lat array, without doing anything
    # clever about the reduced latitudes (they are just centred in a full-width array of NaNs)
    time_slices.append(data_procd)

# Create Iris cube from processed data
cube = create_cube(np.array(time_slices), field, lats)

# Save to file
save(CubeList([cube]), f"{field}.nc", fill_value=np.nan)
