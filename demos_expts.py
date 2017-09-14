{

'helix0_demo':{
#####
####
###
##
#
'tags':['cgmd','tag_structure_repo'],
'metarun':[
{'step':'bilayer','do':'bilayer_control_cgmd','settings':"""
step: bilayer
monolayer top: 90
composition top: {'DOPC':0.64,'DOPS':0.16,'POP2':0.2}
composition bottom: {'POPC':1.0}
"""},
{'step':'protein','do':'martinize','settings':"""
start structure: @structure-repo/proteins/helix0.pdb
"""},
{'step':'adhere','do':'bilayer_protein_adhesion','settings':"""
force field: martini-sources
sources: ['@martini/martini-sources.ff']
placement method: banana
group up: resid 19
group down: resid 7
group origin: resid 7
bilayer structure: s01-bilayer/md.part0001.gro
protein_lattice:|{
	'nrows':1,'ncols':1,
	'lattice_type':'square',
	'space_scale':20,
	'total_proteins':1,
	'protein_shift_up':1.0,
	}
"""}]},

'helix0_flat_demo':{
#####
####
###
##
#
'tags':['cgmd','tag_structure_repo'],
'metarun':[
{'step':'bilayer','do':'bilayer_control_flat','settings':"""
step: bilayer
monolayer top: 90
composition top: {'DOPC':0.64,'DOPS':0.16,'POP2':0.2}
composition bottom: {'POPC':1.0}
"""},
{'step':'protein','do':'martinize','settings':"""
start structure: @structure-repo/proteins/helix0.pdb
"""},
{'step':'adhere','do':'bilayer_protein_adhesion','settings':"""
placement method: banana
group up: resid 19
group down: resid 7
group origin: resid 7
bilayer structure: s01-bilayer/md.part0001.gro
protein_lattice:|{
	'nrows':1,'ncols':1,
	'lattice_type':'square',
	'space_scale':20,
	'total_proteins':1,
	'protein_shift_up':1.0,
	}
"""}]},

'banana_demo':{
#####
####
###
##
#
'tags':['cgmd','tested_2017.09.14','tag_exo70','tag_structure_repo','tag_structure_repo_MISSING'],
'metarun':[
{'step':'protein','do':'martinize','settings':"""
USAGE NOTES:|
    the quintessential banana simulation
    !!! needs:
        starting bilayer v813
        incorporation into a test set that makes the starting bilayer
        switch to multimeric protein class
        ...
    !!! requirements updated on 2017.09.14 during testing
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
"""}]}

}