{

'bilayer_control_multiply':{
#####
####
###
##
#
'metarun':[
{'step':'bilayer','do':'bilayer_control'},
{'step':'large','do':'multiply','settings':"""
step: large
requires: multiply
equilibration: npt-bilayer-short,npt-bilayer
minimize: True
proceed: True
genconf gap: 0.3
nx: 2
ny: 2
"""}]},

'bilayer_control_flat_multiply':{
#####
####
###
##
#
'metarun':[
{'step':'bilayer','do':'bilayer_control_flat'},
{'step':'large','do':'multiply','settings':"""
step: large
requires: multiply
equilibration: npt-bilayer-short,npt-bilayer
minimize: True
proceed: True
genconf gap: 0.3
nx: 2
ny: 2
"""}]},

'bilayer_control':{
#####
####
###
##
#
'script':'bilayer.py',
'params':'parameters.py',
'tags':['cgmd','bilayer','free'],
'extensions':[
    '@extras/*.py',
    '@extras/geometry_tools/*.py',
    'codes/*.py'],
'settings':"""

USAGE NOTES:|
    use this procedure to make a "free" bilayer
    this method packs lipids in vacuum with restraints
    requires restrains generated via a run named "generate_lipidome_restraints"

step:               bilayer

#---SHAPE AND COMPOSITION
shape:              flat            # initial mesh shape flat
aspect:             1.0             # XY proportion for flat bilayers
binsize:            1.2             # grid spacing for initial lipid configuration
monolayer offset:   1.5             # initial distance between leaflets 
monolayer top:      90              # number of lipids in the top leaflet
monolayer bottom:   None            # number of lipids in the bottom leaflet (none for symmetric)

lipid structures:   inputs/martini/library-lipidome-structs  # folder for lipid structures and yaml metadata
landscape metadata: inputs/martini/auto_ff/landscape.json    # colloquial types for different molecules

#---COMPOSITIONS (propotional -- no need to sum to unity)
composition top:    {'DOPC':0.64,'DOPS':0.16,'POP2':0.2}
composition bottom: {'POPC':1.0}

#---SOLVATION
cation:             NA+             # residue name for the cation (must be found in the ff)
anion:              CL-             # residue name for the anion (must be found in the ff)
ionic strength:     0.150           # molar ionic strength
sol:                W               # residue name for water
atom resolution:    cgmd            # either cgmd or aamd
water buffer:       8               # water-other gap distance in Angstroms (avoid waters in bilayer!)
solvent:            martini-water   # water box (must be copied via files)
thickness:          14              # thickness of the box at the solvate step

#---COPY DEPENDENCIES
files:              ['inputs/martini/library-general-structs/martini-water.gro']
sources:            ['inputs/martini/martini-sources.ff','inputs/martini/auto_ff/martini_upright.ff']

#---FORCE FIELD
force field:           martini-sources     # specify the name of the force field (minus ".ff" suffix)
                                           # note non-standard force fields must be copied via "sources"
force field upright:   martini_upright     # force field with "upright" vacuum pack restraints

#---EQUILIBRATION
equilibration:      npt-bilayer
mdp specs:|{
    'group':'cgmd',
    'mdps':{
        'input-em-steep-in.mdp':['minimize'],
        'input-md-vacuum-pack1-eq-in.mdp':['vacuum-packing',{'ref_p':'500 1','nsteps':100000}],
        'input-md-vacuum-pack2-eq-in.mdp':['vacuum-packing',{'nsteps':100000}],
        'input-md-npt-bilayer-short-eq-in.mdp':['npt-bilayer',
            {'pressure':'npt-semiisotropic-weak','dt':0.001,'nsteps':200000}],
        'input-md-npt-bilayer-eq-in.mdp':['npt-bilayer',
            {'restrain':'posre-com-only','pressure':'npt-semiisotropic-weak','dt':0.01}],
        'input-md-in.mdp':[{'pressure':'npt-semiisotropic-weak'}],
        },
    }

"""},

'bilayer_control_flat':{
#####
####
###
##
#
'script':'bilayer-flat.py',
'params':'parameters.py',
'tags':['cgmd','bilayer','free','flat'],
'extensions':[
    '@extras/*.py',
    '@extras/geometry_tools/*.py',
    'codes/*.py'],
'settings':"""

USAGE NOTES:|
    use this procedure to make a "free" bilayer
    this method packs lipids in vacuum with restraints
    requires restrains generated via a run named "generate_lipidome_restraints"
    run the full test with multiply via: 
        make quick clear_lipidome_restraints && make prep generate_lipidome_restraints && make run && \
        make clean sure && make prep bilayer_control_flat_multiply && make metarun

step:               bilayer

#---SHAPE AND COMPOSITION
shape:              flat            # initial mesh shape flat
aspect:             1.0             # XY proportion for flat bilayers
binsize:            1.2             # grid spacing for initial lipid configuration
monolayer offset:   1.5             # initial distance between leaflets 
monolayer top:      90              # number of lipids in the top leaflet
monolayer bottom:   None            # number of lipids in the bottom leaflet (none for symmetric)

lipid structures:   inputs/martini/library-lipidome-structs/  # folder for lipid structures and yaml metadata
landscape metadata: inputs/martini/auto_ff/landscape.json     # colloquial types for different molecules

#---COMPOSITIONS (propotional -- no need to sum to unity)
composition top:    {'DOPC':0.64,'DOPS':0.16,'POP2':0.2}
composition bottom: {'POPC':1.0,'DPPC':0.5}

#---SOLVATION
cation:             NA+             # residue name for the cation (must be found in the ff)
anion:              CL-             # residue name for the anion (must be found in the ff)
ionic strength:     0.150           # molar ionic strength
sol:                W               # residue name for water
atom resolution:    cgmd            # either cgmd or aamd
water buffer:       8               # water-other gap distance in Angstroms (avoid waters in bilayer!)
solvent:            martini-water   # water box (must be copied via files)
thickness:          14              # thickness of the box at the solvate step

#---COPY DEPENDENCIES
files:    ['inputs/martini/library-general-structs/martini-water.gro']
sources:| [
    'inputs/martini/martini-sources.ff',
    'inputs/martini/auto_ff/martini_upright.ff',
    'inputs/martini/auto_ff/martini_upright_alt.ff']

#---FORCE FIELD
force field:           martini-sources      # specify the name of the force field (minus ".ff" suffix)
                                            # note non-standard force fields must be copied via "sources"
force field upright:   martini_upright      # force field with "upright" vacuum pack restraints
force field flat:      martini_upright_alt  # force field with alternate lipids for remaining flat

#---EQUILIBRATION
equilibration:      npt-bilayer
mdp specs:|{
    'group':'cgmd',
    'mdps':{
        'input-em-steep-in.mdp':['minimize'],
        'input-md-vacuum-pack1-eq-in.mdp':['vacuum-packing',{'ref_p':'500 1','nsteps':100000}],
        'input-md-vacuum-pack2-eq-in.mdp':['vacuum-packing',{'nsteps':100000}],
        'input-md-npt-bilayer-short-eq-in.mdp':['npt-bilayer',
            {'pressure':'npt-semiisotropic-weak','dt':0.001,'nsteps':200000}],
        'input-md-npt-bilayer-eq-in.mdp':['npt-bilayer',
            {'restrain':'posre-com-only','pressure':'npt-semiisotropic-weak','dt':0.01}],
        'input-md-in.mdp':[{'restrain':'posre-com-only','pressure':'npt-semiisotropic-weak',
            'nsteps':500000}],
        },
    }

"""},

}
