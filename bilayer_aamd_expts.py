{

'bilayer_control_aamd':{
#####
####
###
##
#
'script':'scripts/bilayer.py',
'params':'parameters-aamd.py',
'tags':['aamd','bilayer','free'],
'extensions':[
    '@extras/*.py',
    '@extras/geometry_tools/*.py',
    'codes/*.py'],
'settings':"""

USAGE NOTES:|
	you MUST run `make go bilayer_control_aamd_restrain clean` first to generate restraints
	use this procedure to make a "free" bilayer
	this method packs lipids in vacuum with restraints
	requires restrains generated via a run named "generate_lipidome_restraints"

step:               bilayer

#---SHAPE AND COMPOSITION
shape:              flat            # initial mesh shape flat
aspect:             1.0             # XY proportion for flat bilayers
binsize:            1.2             # grid spacing for initial lipid configuration
monolayer offset:   1.5             # initial distance between leaflets 
monolayer top:      36              # number of lipids in the top leaflet
monolayer bottom:   None            # number of lipids in the bottom leaflet (none for symmetric)

lipid structures:   @structure-repo/bilayers-aamd/lipid-structures  # folder for lipid structures
landscape metadata: @charmm/landscape.json                          # colloquial types for different molecules

#---COMPOSITIONS (propotional -- no need to sum to unity)
composition top:    {'DOPC':1,'DOPS':1,'PI2P':1,'CHL1':1}
composition bottom: {'POPC':3,'CHL1':1}

#---SOLVATION
cation:             NA              # residue name for the cation (must be found in the ff)
anion:              CL              # residue name for the anion (must be found in the ff)
ionic strength:     0.150           # molar ionic strength
sol:                SOL             # residue name for water
atom resolution:    cgmd            # either cgmd or aamd
water buffer:       8               # water-other gap distance in Angstroms (avoid waters in bilayer!)
solvent:            spc216          # water box (must be copied via files)
thickness:          10              # thickness of the box at the solvate step
                                    # ...be careful with this. you can get widely varying levels of water

#---COPY DEPENDENCIES
files:              []
sources:            ['@charmm/charmm36.ff','@charmm/auto_ff/charmm36_upright.ff']

#---FORCE FIELD
force field:           charmm36            # specify the name of the force field (minus ".ff" suffix)
                                           # note non-standard force fields must be copied via "sources"
force field upright:   charmm36_upright    # force field with "upright" vacuum pack restraints

#---EQUILIBRATION
equilibration: npt-bilayer-short,npt-bilayer
mdp specs:|{
	'group':'aamd',
	'mdps':{
		'input-em-steep-in.mdp':['minimize',{'potential':'verlet','emtol':10}],
		'input-md-vacuum-pack1-eq-in.mdp':[
		'vacuum-packing',{'nsteps':10000,'dt':0.001}],
		'input-md-vacuum-pack2-eq-in.mdp':['vacuum-packing',{'ref_p':'500.0 1.0'}],
		'input-md-vacuum-pack3-eq-in.mdp':['vacuum-packing'],
		'input-md-npt-bilayer-eq-in.mdp':['npt-bilayer-simple',],
		'input-md-npt-bilayer-short-eq-in.mdp':['npt-bilayer',{'dt':0.0002,'nsteps':10000}],
		'input-md-in.mdp':['npt-bilayer'],
		},
	}

"""},

'bilayer_control_aamd_restrain':{
#####
####
###
##
#
'tags':['cgmd','lipidome'],
'params':'@bilayers/parameters.py',
'extensions':[
	'@extras/geometry_tools/*.py'],
'quick':"""

from amx import *
init()
#---note that restraint_maker is in amx/topology_tools.py loaded with amx
restraint_maker()

""",
'settings':"""

USAGE NOTES:|
	this method is designed to pre-make any lipid restraints you might want
	its product can be used by 
		(1) the bilayer maker's vacuum packing
		(2) the flat bilayer maker's leaflet-specific restraints
	the "wants" dictionary specifies the outputs and describes their restraints
	the deposit site holds the automatically generated force fields
	this method completely avoids using the "define posre" flags in GROMACS
	all restraints are explicit, but this means you should avoid "define posre" which will restrain water
	note that this experiment was copied nearly verbatim from the lipidome_expts.py
	the quick script restraint function has been generalized 
		and moved from the martini bundle to topology_tools.py

force field: charmm                      # we always include the force field for the Landscape class
base force field: @charmm/charmm36.ff    # source force field to modify
deposit site: @charmm/auto_ff            # where to write new force fields (keys in wants)

#---specify transformed force field copies
wants:|{
	'charmm36_upright.ff':[
		{'restraints':{'charmm_glycerol':{'z':1000},'charmm_tails':{'z':1000}},
		'naming':'same','which':'lipids'},
		{'restraints':{'sterol_out':{'z':1000},'sterol_in':{'z':1000}},
		'naming':'same','which':'sterols'},
		],
	}

"""},


}