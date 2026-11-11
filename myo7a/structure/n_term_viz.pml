load SM.Nterm.pdb, SM.Nterm;

# chain D only: the ADP HETATM (chain _) also has resi 8
create n_terminal, SM.Nterm and chain D and resi 7-60;
remove SM.Nterm and chain D and resi 7-59;

zoom complete=1;
bg_color white;
hide everything;



################################
set surface_quality, 1;
alter SM.Nterm, b=50;
alter SM.Nterm, q=1;
set gaussian_resolution, 10;

map_new map_SM_Nterm, gaussian, 1, SM.Nterm, 6;
isosurface surf_SM_Nterm, map_SM_Nterm
color deepblue, surf_SM_Nterm;
set transparency, 0.5, surf_SM_Nterm;
set spec_reflect, 0.3;


################################
show cartoon, n_terminal;
color deepblue, n_terminal;

################################
show spheres, n_terminal and resi 26;
color red, n_terminal and resi 26;

