
celestialRadius = 1 # Gotta decide this value.

originalStar = bpy.data.objects['Star'].data # input mesh. Should it be a lamp? Or both?

scene = bpy.context.scene

# Duplicated object
ob = bpy.data.objects.new('Duplicate_Linked', originalStar)

# Link originalStar mesh to this duplicate, or something like that.
scene.collection.objects.link(ob)


# more code if you want...
ob.location = (2,2,2) # place it somewhere else...
