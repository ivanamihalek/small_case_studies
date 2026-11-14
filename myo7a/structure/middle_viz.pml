load SM.middle.pdb, SM.middle;

zoom complete=1;
bg_color white;
hide everything;



################################
set surface_quality, 1;
alter SM.middle, b=50;
alter SM.middle, q=1;
set gaussian_resolution, 10;

map_new map_SM_middle, gaussian, 1, SM.middle, 6;
isosurface surf_SM_middle, map_SM_middle
color deepblue, surf_SM_middle;
set transparency, 0.5, surf_SM_middle;
set spec_reflect, 0.3;


################################
show spheres, resi 1327;
color red,resi 1327;

