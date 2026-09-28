bg_color white
cmd.set("surface_quality", 1)
cmd.alter("all", "b=50")
cmd.alter("all", "q=1")
cmd.set("gaussian_resolution", 7)

#############################
# Trk (TrkA in this file, though TrkB is the actual interaction partner)
load 2ifg_to_1bndA.pdb, trk
remove hetatm
# E and F is the original neurotrophic factor
remove trk and (chain E or chain F)
select trkA, trk and chain A
select trkB, trk and chain B
hide everything, trk

# 1.5 is grid spacing
cmd.map_new("trkA_map", "gaussian", 1.5, "trkA", 2)
cmd.isosurface("trkA_surf","trkA_map")
cmd.color("smudge", "trkA_surf")

cmd.map_new("trkB_map", "gaussian", 1.5, "trkB", 2)
cmd.isosurface("trkB_surf","trkB_map")
cmd.color("smudge", "trkB_surf")


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
     0.841271996,    0.527049541,    0.120299637,\
     0.204698130,   -0.516508162,    0.831454337,\
     0.500355840,   -0.674853981,   -0.542414069,\
    -0.000036307,    0.000384539, -402.546997070,\
    -3.137647629,    2.725387573,   75.008193970,\
   294.024200439,  511.087554932,  -20.000000000 )



#select small_hairpin_A, bdnfA and resi 136-228
#select big_hairpin_A, bdnfA and resi 229-244

#select small_hairpin_B, bdnfB and resi 136-228
#select big_hairpin_B, bdnfB and resi 229-244
