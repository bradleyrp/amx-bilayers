{

'enth_demo':{
#####
####
###
##
#
'tags':['cgmd','tested_2018.04.24','note_EXTRA_REQS'],
'prelude':'make go lipidome clean && make clean sure',
'metarun':[

# step 1: coarse-grain the protein using MARTINI
{'step':'protein','do':'martinize','settings':"""

USAGE NOTES:|
    this is the simplest demo for CGMD protein-bilayer simulations
    it is based on the epsin N-terminal homology domain (ENTH)
    it requires a bilayer_structure and a start_structure
    martinize flags
        the martinize_flags string gets passed directly to martinize.py 
        to add an eleastic network add e.g. "-elastic -ef 500 -el 0.5 -eu 0.9 -ea 0 -ep 0"
        for a custom secondary structure add "-ss CCHHH..." to the martinize_flags string
    choose the martinize force field using the martinize_ff flag below
        by default we use the martini22 force field with extended dihedrals
        however several versions are available (run "python inputs/martini/bin/martinize.py -h" for docs)
    tested with a random secondary structure with elastic restains or with extended dihedrals and DSSP

#! structure-special is a non-public repository
start structure: @structure-special/1H0A-prepped.pdb
martinize flags: -ed
martinize ff: martini22

"""},
# step 2: run the bilayer adhesion procedure
{'step':'adhere','do':'bilayer_protein_adhesion','settings':"""
# run the 288 bilayer first or use a different starting structure
bilayer structure: inputs/bilayer-cgmd-288.gro
placement method: globular_up_down
group up: all
group down: ['resid 1-22','resid 68-72']
protein water gap: 5
protein_lattice:|{
    'nrows':1,'ncols':1,
    'lattice_type':'square',
    'space_scale':20,
    'total_proteins':1,
    'protein_shift_up':1.2,}
equilibration: ['npt-bilayer-short','npt-bilayer']
mdp specs:|{
    'group':'cgmd',
    'mdps':{
        'input-em-steep-in.mdp':['minimize'],
        'input-md-npt-bilayer-short-eq-in.mdp':[{
            'restrain':'posre-com-only','pressure':'npt-semiisotropic-weak',
            'nsteps':50000,'dt':0.001,'groups':'protein','temperature':'protein'}],
        'input-md-npt-bilayer-eq-in.mdp':[{'restrain':'posre-com-only',
            'pressure':'npt-semiisotropic-weak',
            'nsteps':50000,'dt':0.01,'groups':'protein','temperature':'protein'}],
        'input-md-in.mdp':[{'restrain':'posre-com-only','pressure':'npt-semiisotropic-weak',
            'nsteps':50000,'dt':0.04,'groups':'protein','temperature':'protein'}],},}
"""},
]},

'enth_demo_flat_prep':{
#####
####
###
##
#
'tags':['cgmd','tested_2018.04.27'], # 14.4 min, 
# tested with enth_demo_flat as a preface in megatest
'prelude':"make go lipidome clean && make clean sure",
'metarun':[
{'step':'bilayer','do':'bilayer_control_flat','settings':"""
step: bilayer
monolayer top: 90
monolayer bottom: 90
composition top: {'DOPC':0.64,'DOPS':0.16,'POP2':0.2}
composition bottom: {'POPC':1.0}
"""},
{'quick':'table','settings':"""
ready: s01-bilayer/md.part0003.gro
store: inputs/bilayer-cgmd-small-flat.gro
"""},
]},

'enth_demo_flat':{
#####
####
###
##
#
'tags':['cgmd','note_structure_repo','tested_2018.04.27'],
'prelude':"make go lipidome clean && make clean sure",
'metarun':[
{'step':'protein','do':'martinize','settings':"""
#! structure-special is a non-public repository
start structure: @structure-special/1H0A-prepped.pdb
"""},
{'step':'adhere','do':'bilayer_protein_adhesion','settings':"""
USAGE NOTES:|
	Requires inputs/bilayer-cgmd-small.gro via enth_demo_flat_prep
	This procedure is based on testset_bilayer_protein_flat. 	
force field: martini_upright_alt
sources: ['@martini/auto_ff/martini_upright_alt.ff']
bilayer structure: inputs/bilayer-cgmd-small-flat.gro
placement method: globular_up_down
group up: all
group down: ['resid 1-22','resid 68-72']
protein water gap: 5
protein_lattice:|{
    'nrows':1,'ncols':1,
    'lattice_type':'square',
    'space_scale':20,
    'total_proteins':1,
    'protein_shift_up':1.2,}

# EQUILIBRATION
equilibration: ['npt-bilayer-short','npt-bilayer']
mdp specs:|{
    'group':'cgmd',
    'mdps':{
        'input-em-steep-in.mdp':['minimize'],
        'input-md-npt-bilayer-short-eq-in.mdp':[{'restrain':'posre-com-only',
            'pressure':'npt-semiisotropic-weak',
            'nsteps':500000,'groups':'protein','temperature':'protein','dt':0.001}],
        'input-md-npt-bilayer-eq-in.mdp':[{'restrain':'posre-com-only',
            'pressure':'npt-semiisotropic-weak',
            'nsteps':500000,'groups':'protein','temperature':'protein','dt':0.01}],
        'input-md-in.mdp':[{'restrain':'posre-com-only','pressure':'npt-semiisotropic-weak',
            'nsteps':500000,'groups':'protein','temperature':'protein'}],},}

"""},
]},

'bilayer_protein_adhesion_custom':{
#####
####
###
##
#
'script':'scripts/bilayer-protein.py',
'params':'parameters.py',
'tags':['cgmd','tested_2017.09.15',],
'extensions':[
    '@extras/*.py',
    '@extras/geometry_tools/*.py',
    'codes/*.py'],
'settings':"""

USAGE NOTES:|
    requires custom top, itp, and gro files supplied to the "prepared protein" section below
    requires a pre-equilibrated, flat bilayer added in the bilayer_structure flag below
    run "make go lipidome" once before using this step
    this is adapted from a run called bilayer_protein_adhesion 
        which typically follows a protein coarse-graining or topology generation step in a metarun
        see the notes below for how we replace the previous step in this settings block
    tested via: ran enth_demo and took key files from the martinize step (s01-protein) and moved them 
        to inputs/custom_protein to simulate a custom construction of a new protein system 
    note that almost all options for creating a MARTINI protein on a bilayer are easily implemented with
        the enth_demo above, which passes any flags you want through to the martinize step
        which is a crude wrapper around the martinize.py script supplied by Marrink et al

step: adhere # step folder name
bilayer structure: @structure-repo/bilayers-cgmd/bilayer-cgmd-288.gro # starting bilayer structure
landscape metadata: @martini/auto_ff/landscape.json # colloquial types for different molecules
sources: ['@martini/auto_ff/martini_upright_alt.ff'] # requires files from "make go lipidome"

#---PREPARED PROTEIN
itp: ['Protein.itp']
files: ['inputs/custom_protein/Protein.itp']
protein prepared:|{
    'gro':'inputs/custom_protein/protein.gro',
    'top':'inputs/custom_protein/protein.top',
    }

#---FORCE FIELD
force field:  martini_upright_alt # specify the name of the force field (minus ".ff" suffix)
                                  # note non-standard force fields must be copied via "sources"

#---PROTEIN POSITION INSTRUCTIONS
placement method: globular_up_down          # several distinct methods to choose from (see other expts)
group up: all                               # group of selections that points up
group down: ['resid 1-22','resid 68-72']    # group which faces the bilayer
group origin: ['resid 1-22','resid 68-72']  # reference or pivot for other use-cases

#---PROTEIN LATTICE DEFINITIONS
protein_lattice:|{
    'nrows':1,'ncols':1,
    'lattice_type':'square',
    'space_scale':20,
    'total_proteins':1,
    'protein_shift_up':1.0,}

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
            'nsteps':500000,'groups':'protein','temperature':'protein'}],},}

"""},

'banana_demo':{
#####
####
###
##
#
'tags':['cgmd','tested_2017.09.14','dev','note_EXTRA_REQS',
# 'note_exo70',
],
'metarun':[
{'step':'protein','do':'martinize','settings':"""
USAGE NOTES:|
    still under active development
    the quintessential banana simulation
    requires additional files not included in the structure-repo
        inputs/structure-repo/bilayers-cgmd/v813.gro
        inputs/structure-repo/proteins/IRSp53_modeller_updated.pdb

start structure: @structure-repo/proteins/IRSp53_modeller_updated.pdb
martinize path: @martini/bin/martinize.py
dssp path: @martini/bin/dssp-2.0.4-linux-amd64
step: protein
martinize flags: -ed
martinize ff: martini22
"""},
{'step':'adhere','do':'bilayer_protein_adhesion','settings':"""
force field: martini_upright_alt
sources: ['@martini/auto_ff/martini_upright_alt.ff']
placement method: banana
group up: ['resid 470']
group down: ['resid 110','resid 280']
group origin: ['resid 110','resid 280']
bilayer structure: @structure-repo/bilayers-cgmd/v813.gro
protein water gap: 3
protein_lattice:|{'banana_direction':'x','nrows':1,'ncols':1,'lattice_type':'square',
    'space_scale':20,'total_proteins':1,'protein_shift_up':3.5}
equilibration: short
mdp specs:|{
    'group':'cgmd',
    'mdps':{
        'input-em-steep-in.mdp':['minimize'],
        'input-md-short-eq-in.mdp':[{'restrain':'posre-com-only','pressure':'npt-semiisotropic-weak',
            'nsteps':10000,'groups':'protein','temperature':'protein','dt':0.005}],
        'input-md-in.mdp':[{'restrain':'posre-com-only',
            'nsteps':50000,'groups':'protein','temperature':'protein'}],
        },
    }
two protein hack: True
"""}]},

'bilayer_control_cholesterol_demo':{
#####
####
###
##
#
'tags':['cgmd','dev'],
'metarun':[
{'step':'bilayer','do':'bilayer_control_flat','settings':"""

USAGE NOTES:|
    under development
        currently experiencing stability problems at npt-bilayer
        increased from 196 to 625

monolayer top:      196
monolayer bottom:   196

composition top:    {'DOPC':0.5,'DOPS':0.1,'POP2':0.2,'CHOL':0.2}
composition bottom: {'POPC':0.8,'CHOL':0.2}

#---EQUILIBRATION
equilibration: ['npt-bilayer-short1','npt-bilayer-short2','npt-bilayer']
#---EQUILIBRATION
mdp specs:|{
    'group':'cgmd',
    'mdps':{
        'input-em-steep-in.mdp':['minimize'],
        'input-md-vacuum-pack1-eq-in.mdp':['vacuum-packing',{'ref_p':'500 1','nsteps':100000}],
        'input-md-vacuum-pack2-eq-in.mdp':['vacuum-packing',{'nsteps':200000}],
        'input-md-npt-bilayer-short1-eq-in.mdp':['npt-bilayer',{
            'restrain':'posre-com-only','virtual-sites':'standard','pressure':'npt-semiisotropic',
            'dt':0.001,'nsteps':100000}],
        'input-md-npt-bilayer-short2-eq-in.mdp':['npt-bilayer',{
            'restrain':'posre-com-only','virtual-sites':'standard','pressure':'npt-semiisotropic',
            'dt':0.005,'nsteps':100000}],
        'input-md-npt-bilayer-eq-in.mdp':['npt-bilayer',{
            'restrain':'posre-com-only','virtual-sites':'standard','pressure':'npt-semiisotropic',
            'dt':0.01}],
        'input-md-in.mdp':[{
            'virtual-sites':'standard','restrain':'posre-com-only','pressure':'npt-semiisotropic-weak',
            'nsteps':500000}],},}

"""},
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

}
