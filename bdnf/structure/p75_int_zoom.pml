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


cmd.map_new("bdnfB_map", "gaussian", 1.5, "bdnfB", 2)
cmd.isosurface("bdnfB_surf","bdnfB_map")
cmd.color("palecyan", "bdnfB_surf")

extract small_hairpin_A, bdnfA and resi 136-200
extract big_hairpin_A, bdnfA and resi 201-244

cmd.map_new("smallA_map", "gaussian", 1.2, "small_hairpin_A", 2)
cmd.isosurface("smallA_surf", "smallA_map")
set transparency, 0.5, smallA_surf
cmd.color("firebrick", "smallA_surf")
show cartoon, small_hairpin_A
color firebrick, small_hairpin_A


cmd.map_new("bigA_map", "gaussian", 1.2, "big_hairpin_A", 2)
cmd.isosurface("bigA_surf", "bigA_map")
cmd.color("firebrick", "bigA_surf")
set transparency, 0.5, bigA_surf


show cartoon, big_hairpin_A
color blue,  big_hairpin_A



#############################
set_view (\
     0.698906660,    0.615534067,    0.364207923,\
    -0.105547242,   -0.414887071,    0.903730094,\
     0.707381845,   -0.670064569,   -0.224999204,\
    -0.000015363,   -0.000033483, -318.124572754,\
     9.024620056,  -18.997974396,   50.423500061,\
   282.127410889,  354.127471924,  -20.000000000 )
