{

'banana_v1':{
#####
####
###
##
#
'metarun':[
{'step':'protein','do':'martinize','settings':"""
USAGE NOTES:|
    the quintessential banana simulation
    !!! needs:
        starting bilayer v813
        incorporation into a test set that makes the starting bilayer
        switch to multimeric protein class
        ...
start structure: @structure-repo/proteins/IRSp53_modeller_updated.pdb
martinize path: @martini/bin/martinize.py
dssp path: @martini/bin/dssp-2.0.4-linux-amd64
step: protein
martinize flags: -ed
martinize ff: martini22
"""},
{'step':'adhere','do':'bilayer_protein_adhesion','settings':"""
USAGE NOTES:|
force field: martini_upright_alt
sources: ['@martini/auto_ff/martini_upright_alt.ff']
placement method: banana
group up: ['resid 470']
group down: ['resid 110','resid 280']
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
            'nsteps':50000,'groups':'protein','temperature':'protein','dt':0.005}],
        'input-md-in.mdp':[{'restrain':'posre-com-only',
            'nsteps':100000,'groups':'protein','temperature':'protein'}],
        },
    }
two protein hack: True
"""}]}

}