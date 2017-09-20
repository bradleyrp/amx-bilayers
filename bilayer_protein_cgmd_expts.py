{

'bilayer_protein_adhesion':{
#####
####
###
##
#
'script':'scripts/bilayer-protein.py',
'params':'parameters.py',
'tags':['cgmd','bilayer','protein-bilayer','tag_support','tested_2017.09.14'],
'extensions':[
    '@extras/*.py',
    '@extras/geometry_tools/*.py',
    'codes/*.py'],
'settings':"""

USAGE NOTES:|
    this is a CORE automacs feature used by several metaruns
    see also "bilayer_protein_adhesion_custom" for the standalone version
	requires a pre-equilibrated, flat bilayer
    gets the protein details automatically from state.protein_prepared
    bilayer structure and force field are coded below, however you can override them in a metarun
    since you can override in a metarun, we have only one adhesion routine

step: adhere # name the step (folder)
bilayer structure: inputs/previous/bilayer.gro # starting bilayer structure
landscape metadata: @martini/auto_ff/landscape.json # colloquial types for different molecules

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
	'protein_shift_up':3.0,}

sol: W                         # !!! remove this to the landscape
ionic strength: 0.150          # box is reionized after adding protein to ensure neutral
protein water gap: 3           # distance in Angstroms between protein in water after addition
reionize method: ions          # select "ions" to replace only ions and "solvent" to replace water too
cation: NA+                    # choose a new cation for adding back counterions
anion: CL-                     # choose a new anion for adding back counterions

#---EQUILIBRATION
equilibration: ['npt-bilayer']
mdp specs:|{
    'group':'cgmd',
    'mdps':{
        'input-em-steep-in.mdp':['minimize'],
        'input-md-npt-bilayer-eq-in.mdp':[{'restrain':'posre-com-only','pressure':'npt-semiisotropic-weak',
            'nsteps':500000,'groups':'protein','temperature':'protein','dt':0.01}],
        'input-md-in.mdp':[{'restrain':'posre-com-only','pressure':'npt-semiisotropic-weak',
            'nsteps':500000,'groups':'protein','temperature':'protein'}],},}

"""},

}