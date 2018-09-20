#!/usr/bin/env python

from amx import *

init()
make_step(settings.step)
write_mdp()
structure = GMXStructure(state.here+state.incoming)
state.composition = structure.detect_composition()
detect_lipids(state.incoming)
restuff(
	structure='vacuum-packed',
	gro='solvate-dry',
	tpr=get_last('tpr'),
	ndx=get_last('ndx'))
solvate(
	structure='solvate-dry',
	gro='solvate-big')
water_solvate_exact('solvate-big.gro','solvate.gro')
structure = GMXStructure(state.here+'solvate.gro')
state.composition = structure.detect_composition()

write_top('solvate.top')
minimize('solvate',restraints=True)
remove_jump(
	structure='solvate-minimized',
	tpr='em-solvate-steep',
	gro='solvate-nojump')

remove_jump(
	structure='solvate-minimized',
	tpr='em-solvate-steep',
	gro='solvate-nojump')

counterions('solvate-nojump','solvate',restraints=True,
	counts=state.genion_overrides)

counterion_renamer('counterions')
write_top('counterions.top')
minimize('counterions',restraints=True)
remove_jump(
	structure='counterions-minimized',
	tpr='em-counterions-steep',
	gro='counterions-nojump')
bilayer_middle(
	structure='counterions-nojump',
	gro='system')
bilayer_sorter(
	structure='system',
	ndx='system-groups')
# relax with fixed area and upright restraints
write_top('system.top')
# note no persistent restraint positions
equilibrate(
	seq=state.equilibrate_restrain,
	restraints=True,
	groups='system-groups',
	stages_only=True)
# release restraints and continue
state.force_field = settings.force_field
move_file('system.top','system-restrain.top')
write_top('system.top')
equilibrate(
	structure=state.equilibrate_restrain_final,
	seq=state.equilibrate_free,
	groups='system-groups',
	stages_only=False)
