# To run this file from command line:
# > pymol conservation_on_the_structure.pml
# To run from pymol itself, open pymol and find conservation_on_the_structure.pml
# in the 'File'->'Run Script ...' dropdown menu.

# If you do not have pymol installed, check
# https://pymolwiki.org/index.php/Windows_Install
# or the links for other platforms at the bottom of that page.

load 1bnd.renumbered.pdb, the_whole_thing
zoom complete=1
bg_color white
hide everything
select chainA, the_whole_thing and chain A and polymer
color white,  chainA
show cartoon, chainA
show spheres, chainA
set_color c0 = [1, 0.83, 0.17]
set_color c1 = [1, 0, 0]
set_color c2 = [0.73, 0, 0]
set_color c3 = [0.47, 0, 0]
set_color c4 = [0.21, 0, 0]
set_color c5 = [0.0625, 0.0625, 0.0625]
set_color c6 = [0.125, 0.125, 0.125]
set_color c7 = [0.1875, 0.1875, 0.1875]
set_color c8 = [0.25, 0.25, 0.25]
set_color c9 = [0.3125, 0.3125, 0.3125]
set_color c10 = [0.375, 0.375, 0.375]
set_color c11 = [0.4375, 0.4375, 0.4375]
set_color c12 = [0.5, 0.5, 0.5]
set_color c13 = [0.5625, 0.5625, 0.5625]
set_color c14 = [0.625, 0.625, 0.625]
set_color c15 = [0.6875, 0.6875, 0.6875]
set_color c16 = [0.75, 0.75, 0.75]
set_color c17 = [0.8125, 0.8125, 0.8125]
set_color c18 = [0.875, 0.875, 0.875]
set_color c19 = [0.9375, 0.9375, 0.9375]
set_color c20 = [1, 1, 1]
color c4, resid 217 and chain A
color c4, resid 161 and chain A
color c20, resid 146 and chain A
color c20, resid 162 and chain A
color c8, resid 230 and chain A
color c4, resid 148 and chain A
color c13, resid 187 and chain A
color c4, resid 197 and chain A
color c9, resid 244 and chain A
color c9, resid 139 and chain A
color c16, resid 223 and chain A
color c14, resid 207 and chain A
color c14, resid 157 and chain A
color c4, resid 181 and chain A
color c10, resid 182 and chain A
color c4, resid 234 and chain A
color c4, resid 145 and chain A
color c6, resid 212 and chain A
color c16, resid 211 and chain A
color c6, resid 149 and chain A
color c4, resid 167 and chain A
color c20, resid 151 and chain A
color c10, resid 202 and chain A
color c4, resid 136 and chain A
color c20, resid 201 and chain A
color c4, resid 152 and chain A
color c14, resid 240 and chain A
color c13, resid 173 and chain A
color c16, resid 192 and chain A
color c16, resid 191 and chain A
color c14, resid 138 and chain A
color c4, resid 239 and chain A
color c4, resid 235 and chain A
color c4, resid 183 and chain A
color c10, resid 213 and chain A
color c16, resid 144 and chain A
color c4, resid 227 and chain A
color c4, resid 203 and chain A
color c9, resid 153 and chain A
color c14, resid 193 and chain A
color c20, resid 172 and chain A
color c14, resid 171 and chain A
color c4, resid 238 and chain A
color c4, resid 140 and chain A
color c16, resid 163 and chain A
color c10, resid 236 and chain A
color c14, resid 177 and chain A
color c20, resid 222 and chain A
color c9, resid 221 and chain A
color c9, resid 178 and chain A
color c20, resid 150 and chain A
color c4, resid 229 and chain A
color c4, resid 200 and chain A
color c14, resid 225 and chain A
color c4, resid 241 and chain A
color c13, resid 176 and chain A
color c14, resid 190 and chain A
color c16, resid 242 and chain A
color c4, resid 237 and chain A
color c4, resid 180 and chain A
color c4, resid 164 and chain A
color c9, resid 210 and chain A
color c4, resid 228 and chain A
color c4, resid 179 and chain A
color c7, resid 175 and chain A
color c20, resid 194 and chain A
color c20, resid 154 and chain A
color c13, resid 226 and chain A
color c4, resid 204 and chain A
color c13, resid 160 and chain A
color c6, resid 214 and chain A
color c9, resid 143 and chain A
color c4, resid 232 and chain A
color c4, resid 184 and chain A
color c4, resid 231 and chain A
color c9, resid 199 and chain A
color c16, resid 174 and chain A
color c4, resid 195 and chain A
color c9, resid 155 and chain A
color c13, resid 209 and chain A
color c16, resid 220 and chain A
color c16, resid 159 and chain A
color c9, resid 205 and chain A
color c4, resid 137 and chain A
color c4, resid 219 and chain A
color c4, resid 141 and chain A
color c4, resid 142 and chain A
color c4, resid 166 and chain A
color c4, resid 215 and chain A
color c8, resid 185 and chain A
color c4, resid 233 and chain A
color c16, resid 189 and chain A
color c7, resid 168 and chain A
color c4, resid 198 and chain A
color c4, resid 206 and chain A
color c4, resid 156 and chain A
color c14, resid 224 and chain A
color c4, resid 196 and chain A
color c16, resid 243 and chain A
color c9, resid 170 and chain A
color c4, resid 158 and chain A
color c4, resid 208 and chain A
color c4, resid 218 and chain A
color c4, resid 186 and chain A
color c6, resid 216 and chain A
color c14, resid 165 and chain A
color c7, resid 169 and chain A
color c20, resid 188 and chain A
color c4, resid 147 and chain A
select heteroatoms,  hetatm and not solvent
select other_chains, not chain A
select struct_water, solvent and chain A
select metals, symbol  mg+ca+fe+zn+na+k+mn+cu+ni+cd+i
# cartoon putty 
show  cartoon,  other_chains
hide  spheres,   heteroatoms
show  sticks,   heteroatoms
show  spheres,  struct_water
show  spheres,  metals
color palecyan, struct_water
color lightteal, other_chains or heteroatoms
color magenta, metals
zoom  chain A
select poorly_scoring, resid 146+162+187+223+207+157+182+211+151+202+201
select poorly_scoring, poorly_scoring or resid 240+173+192+191+138+213+144+193+172+171+163
select poorly_scoring, poorly_scoring or resid 236+177+222+150+225+176+190+242+194+154+226
select poorly_scoring, poorly_scoring or resid 160+174+209+220+159+189+224+243+165+188
select poorly_scoring, poorly_scoring and chainA
deselect
