load SM.Nterm.pdb, SM.Nterm;

# chain D only: the ADP HETATM (chain _) also has resi 8
create n_terminal, SM.Nterm and chain D and resi 7-60;
remove SM.Nterm and chain D; # we do not need it here

zoom complete=1;
bg_color white;
hide everything;



################################
set surface_quality, 1;
alter n_terminal, b=50;
alter n_terminal, q=1;
set gaussian_resolution, 10;

map_new nterm_map, gaussian, 1, SM.Nterm, 6;
isosurface nterm_surf, nterm_map
color deepblue, nterm_surf;
set transparency, 0.5, nterm_surf;
set spec_reflect, 0.3;


