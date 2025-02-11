import skeletor as sk # Version 1.2.1
import fafbseg # Version 3.0.10
from fafbseg import flywire

### Note for the code to work you will need access to the flywire dataset, which is accessible after training: 
# The below commmand only needs to be done once, where you set your token, this creates a seperate cloudvolume file storing your token for future use. 
#fafbseg.flywire.set_chunkedgraph_secret("INSERT YOUR TOKEN HERE")

#Before getting the mesh of your neurons assure that it is the most up to date version within the dataset
x = flywire.update_ids(720575940606954507)
print(x)

#Here insert the ID of the neuron of interest
neuron_mesh = flywire.get_mesh_neuron(720575940606954507, threads=1)

# Fixing any small disconnections 
fixed_mesh = sk.pre.fix_mesh(neuron_mesh, remove_disconnected=25, inplace=False)

#Keeping the step_size low usually creates more nodes, but maintains a high level of detail for the dendrites of the neuron, 
# thicker processes will end up having issues and need to be edited by hand
skel = sk.skeletonize.by_wavefront(fixed_mesh, waves=1, step_size=1) 

#Cleans the skeleton from excess branching, or recenters nodes outside the mesh
skel_clean= sk.post.clean_up(skel, mesh = fixed_mesh)

## Saving the skeleton in a form to then edit, note that skeleton is not in microns and will need to be adjusted post process, but can be done using NeuTube.
skel_clean.save_swc('720575940606954507_skeletonized.swc')
