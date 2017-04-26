{

'test-h0':{
#####
####
###
##
#
'metarun':[
{'step':'bilayer','do':'bilayer_control','settings':"""
step: bilayer
monolayer top: 90
composition top: {'DOPC':0.64,'DOPS':0.16,'POP2':0.2}
composition bottom: {'POPC':1.0}
"""},
{'step':'protein','do':'martinize','settings':"""
start structure: inputs/helix0.pdb
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

'test_helix0_flat':{
#####
####
###
##
#
'metarun':[
{'step':'bilayer','do':'bilayer_control_flat','settings':"""
step: bilayer
monolayer top: 90
composition top: {'DOPC':0.64,'DOPS':0.16,'POP2':0.2}
composition bottom: {'POPC':1.0}
"""},
{'step':'protein','do':'martinize','settings':"""
start structure: inputs/helix0.pdb
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

}