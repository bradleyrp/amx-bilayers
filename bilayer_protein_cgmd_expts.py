{

'bilayer_protein_dev':{
#####
####
###
##
#
'metarun':[
{'step':'protein','do':'martinize'},
{'step':'adhere','do':'bilayer_protein_adhesion_dev'}]},

'bilayer_protein_adhesion_dev':{
#####
####
###
##
#
'script':'scripts/bilayer-protein.py',
'params':'parameters.py',
'tags':['cgmd','bilayer','protein-bilayer'],
'extensions':[
    '@extras/*.py',
    '@extras/geometry_tools/*.py',
    'codes/*.py'],
'settings':"""

USAGE NOTES:|
	requires a pre-equilibrated, flat bilayer
	...

step: adhere

########### HACKED. GET THE PROTEIN FROM THE LAST STEP AUTOMATICALLY vvvvvvv
protein structure: s01-protein/protein.gro      # starting protein structure 
protein topology: s01-protein/protein.top       # starting protein topology
bilayer structure: inputs/previous/bilayer.gro  # starting bilayer structure

landscape metadata: @martini/auto_ff/landscape.json    # colloquial types for different molecules

#---COPY DEPENDENCIES
files:              []
sources:            ['@martini/auto_ff/martini_upright_alt.ff']

#---FORCE FIELD
force field:  martini_upright_alt # specify the name of the force field (minus ".ff" suffix)
                                  # note non-standard force fields must be copied via "sources"

#---PROTEIN POSITION INSTRUCTIONS
placement method: banana                # several distinct methods to choose from (see other expts)
group up: resid 258                     # banana method: this group points up *after* laying the banana flat
group down: ['resid 240','resid 652']   # banana method: this group points down *after* laying the banana flat
group origin: ['resid 240','resid 652'] # reference/pivot point for the protein (set to origin)

#---PROTEIN LATTICE DEFINITIONS
protein_lattice:|{
	'nrows':1,'ncols':1,
	'lattice_type':'square',
	'space_scale':20,
	'total_proteins':1,
	'protein_shift_up':3.0,
	}

sol: W                         # !!! remove this to the landscape
ionic strength: 0.150          # box is reionized after adding protein to ensure neutral
protein water gap: 3           # distance in Angstroms between protein in water after addition
reionize method: ions          # select "ions" to replace only ions and "solvent" to replace water too
cation: NA+                    # choose a new cation for adding back counterions
anion: CL-                     # choose a new anion for adding back counterions

mdp specs:|{
    'group':'cgmd',
    'mdps':{
        'input-em-steep-in.mdp':['minimize'],
        'input-md-in.mdp':[{'restrain':'posre-com-only','pressure':'npt-semiisotropic-weak',
            'nsteps':500000,'groups':'protein','temperature':'protein'}],
        },
    }

"""},

}