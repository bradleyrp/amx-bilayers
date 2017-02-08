{

'bilayer_protein_aamd_banana':{
#####
####
###
##
#
'metarun':[
{'step':'protein','do':'bilayer_protein_topology_only','settings':"""

#---PROTEIN START STRUCTURE
start structure: inputs/structure-repo/proteins/gel-prepped.pdb

step: protein                       # name of the folder is s01-protein
force field: charmm27               # which gromacs-standard force-field to use (see pdb2gmx list)
water: tip3p                        # which water model (another questi

"""},
{'step':'protein','do':'bilayer_protein_adhesion_aamd','settings':"""

bilayer structure: inputs/structure-repo/bilayers-aamd/bilayer-pip2-20pct-top-ow-dopc-dops-4to1.gro

placement method: banana                 # several distinct methods to choose from (see other expts)
group up: resid 161                      # banana method: this group points up *after* laying the banana flat
group down: ['resid 150','resid 169']    # banana method: this group points down *after* laying the banana flat
group origin: ['resid 150','resid 169']  # reference/pivot point for the protein (set to origin)

"""},
]},

'bilayer_protein_adhesion_aamd':{
#####
####
###
##
#
'script':'scripts/bilayer-protein.py',
'params':'parameters-aamd.py',
'tags':['aamd','bilayer','protein-bilayer'],
'extensions':[
    '@extras/*.py',
    '@extras/geometry_tools/*.py',
    'codes/*.py'],
'settings':"""

USAGE NOTES:|
	requires a pre-equilibrated, flat bilayer
    gets the protein details automatically from state.protein_prepared
    bilayer structure and force field are coded below, however you can override them in a metarun
    since you can override in a metarun, we have only one adhesion routine
	
step: adhere

#---STARTING STRUCTURES
bilayer structure: None

#---COPY DEPENDENCIES
files:              []
sources:            ['@charmm/charmm36.ff','@charmm/lipid-tops']

#---FORCE FIELD
force field: charmm # note non-standard force fields must be copied via "sources"
extra itps: 'lipid-tops/lipid*itp'
landscape metadata: @charmm/landscape.json # colloquial types for different molecules

#---PROTEIN POSITION INSTRUCTIONS
placement method: banana                 # several distinct methods to choose from (see other expts)
group up: resid 161                      # banana method: this group points up *after* laying the banana flat
group down: ['resid 150 or resid 169']   # banana method: this group points down *after* laying the banana flat
group origin: ['resid 150 or resid 169'] # reference/pivot point for the protein (set to origin)

#---PROTEIN LATTICE DEFINITIONS
protein_lattice:|{
	'nrows':1,'ncols':1,
	'lattice_type':'square',
	'space_scale':20,
	'total_proteins':1,
	'protein_shift_up':0.0,
	}

sol: SOL                       # !!! remove this to the landscape
ionic strength: 0.150          # box is reionized after adding protein to ensure neutral
protein water gap: 3           # distance in Angstroms between protein in water after addition
reionize method: ions          # select "ions" to replace only ions and "solvent" to replace water too
cation: NA                     # choose a new cation for adding back counterions
anion: CL                      # choose a new anion for adding back counterions

mdp specs:|{
    'group':'aamd',
    'mdps':{
         'input-em-steep-in.mdp':['minimize',{'potential':'verlet','emtol':10}],
         'input-md-npt-bilayer-eq-in.mdp':['npt-bilayer-protein-simple'],
         'input-md-in.mdp':['npt-bilayer-protein'],
        },
    }

"""},

}