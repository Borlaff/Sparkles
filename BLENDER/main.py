originalStar = bpy.data.objects['Star'].data # input mesh. Should it be a lamp? Or both?

scene = bpy.context.scene

# Duplicated object
ob = bpy.data.objects.new('Duplicate_Linked', originalStar)

# Link originalStar mesh to this duplicate, or something like that.
scene.collection.objects.link(ob)


# Next: import the csv into a pandas dataframe

import pandas as pd

#df = pd.read_csv("low_res_euclid_gaia_map_v1.csv", index_col = 0) # Can't find this from blender

fullPath = "/home/mario/Work/Sparkles/BLENDER/low_res_euclid_gaia_map_v1.csv"
df = pd.read_csv(fullPath, index_col = 0) # Can't find this from blender
print(df.head())

# Drop some columns.
df = df.drop(columns=["ra", "dec", "ipix"])
print(df.head())

# Now convert to cartesian
from astropy import units as u
from astropy.coordinates import SkyCoord

# Testing with first coordinates
raTest = df.iloc[0]["ra_bary"] * u.degree
decTest = df.iloc[0]["dec_bary"] * u.degree
celestialSphereRadius = 100000000 * u.km

c = SkyCoord(ra = raTest, dec = decTest, distance = celestialSphereRadius, frame='icrs') # ICRS? Sure?

x = c.cartesian.x
y = c.cartesian.y
z = c.cartesian.z

print(x, y, z)

ob.location = (2,2,2) # place it somewhere else...
