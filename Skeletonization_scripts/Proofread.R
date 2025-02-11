#Pulling KC cells morphologies from flywire, this is mainly for getting the meshes themselves and then comparing the generated skeletons back to the original mesh.
{
library(fafbseg)
library(reticulate)
library(elmr)
library(gargle)
library(nat.templatebrains)
}
choose_segmentation(release = 'flywire31')

#Downloading of KC_g meshes to visualize in meshlab
KC_g <- c("720575940606954507")
save_cloudvolume_meshes(KC_g, savedir="KC_Github")

KC_g_1_mesh=read_cloudvolume_meshes("720575940606954507")

# Read skeletons and plotting (get new codex skeletons they have done a better job and now updated it.)
KC_Gamma_skel_1 = read.neurons("720575940606954507_skeletonized_model.swc")

#Getting the calyx and peduncle to estimate SIZ location
Ca_L = as.mesh3d(FAFB14NP.surf, "CA_R")
PED_L = as.mesh3d(FAFB14NP.surf, "PED_R")

# Plot the combined mesh with different colors
plot3d(Ca_L, col = 'green', add = TRUE, alpha = 0.2)

# Plot the PED mesh in red on the same plot
plot3d(PED_L, col = 'red', add = TRUE, alpha = 0.2)
plot3d(FAFB14.surf, col = 'gray', alpha = 0.2)
plot3d(KC_g_1_mesh, col = 'black', add = TRUE)

plot3d(KC_Gamma_skel_1, col = 'blue', WithNodes =F, size =2)