import pandas as pd

originalStar = bpy.data.objects['Star'].data # input mesh. Should it be a lamp? Or both?

scene = bpy.context.scene

# Duplicated object
# These should all go to a collection
ob = bpy.data.objects.new('DuplicateStar', originalStar)

# Link originalStar mesh to this duplicate, or something like that.
scene.collection.objects.link(ob)


# Next: import the csv into a pandas dataframe


filepath = "/home/mario/Work/Sparkles/BLENDER/low_res_euclid_gaia_map_v1.csv" # Has to be full path!
df = pd.read_csv(filepath, index_col = 0)
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

xCoord = c.cartesian.x
yCoord = c.cartesian.y
zCoord = c.cartesian.z

print(xCoord, yCoord, zCoord)
print(type(xCoord.value))

ob.location = (xCoord.value, yCoord.value, zCoord.value)    # This works
