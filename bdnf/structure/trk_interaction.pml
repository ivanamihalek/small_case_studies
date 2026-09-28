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

cmd.map_new("bdnfA_map", "gaussian", 1.5, "bdnfA", 2)
cmd.isosurface("bdnfA_surf","bdnfA_map")
cmd.color("firebrick", "bdnfA_surf")

cmd.map_new("bdnfB_map", "gaussian", 1.5, "bdnfB", 2)
cmd.isosurface("bdnfB_surf","bdnfB_map")
cmd.color("palecyan", "bdnfB_surf")


#############################
set_view (\
     0.791300774,    0.604774892,    0.089924708,\
     0.221825644,   -0.421013266,    0.879511356,\
     0.569767118,   -0.676011384,   -0.467303425,\
    -0.000061989,    0.000121191, -447.818328857,\
    -3.098142147,    2.037345886,   72.033493042,\
   339.294982910,  556.358276367,  -20.000000000 )




#select small_hairpin_A, bdnfA and resi 136-228
#select big_hairpin_A, bdnfA and resi 229-244

#select small_hairpin_B, bdnfB and resi 136-228
#select big_hairpin_B, bdnfB and resi 229-244
