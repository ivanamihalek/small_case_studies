bg_color white
cmd.set("surface_quality", 1)
cmd.alter("all", "b=50")
cmd.alter("all", "q=1")
cmd.set("gaussian_resolution", 7)

#############################
# p75
load 3buk_to_1bndA.pdb,  p75
remove hetatm
remove p75 and (chain A or chain B)
select p75C, p75 and chain C
select p75D, p75 and chain D
hide everything, p75

# 1.5 is grid spacing
cmd.map_new("p75C_map", "gaussian", 1.5, "p75C", 2)
cmd.isosurface("p75C_surf","p75C_map")
cmd.color("tv_orange", "p75C_surf")

cmd.map_new("p75D_map", "gaussian", 1.5, "p75D", 2)
cmd.isosurface("p75D_surf","p75D_map")
cmd.color("tv_orange", "p75D_surf")

#############################
# BDNF
load 1bndA.pdb, bdnfA
load 1bndA_to_1bndB.pdb, bdnfB
remove hetatm
hide everything, bdnfA or bdnfB

cmd.map_new("bdnfA_map", "gaussian", 1.5, "bdnfA", 2)
cmd.isosurface("bdnfA_surf","bdnfA_map")
cmd.color("firebrick", "bdnfA_surf")

cmd.map_new("bdnfB_map", "gaussian", 1.5, "bdnfB", 2)
cmd.isosurface("bdnfB_surf","bdnfB_map")
cmd.color("palecyan", "bdnfB_surf")


#############################
set_view (\
     0.714501262,    0.678893685,    0.169089660,\
     0.109646305,   -0.347351998,    0.931302667,\
     0.690988779,   -0.646877348,   -0.322621375,\
    -0.000021860,   -0.000041448, -405.124908447,\
     2.457776070,  -15.563024521,   54.230216980,\
   369.127227783,  441.127288818,  -20.000000000 )
