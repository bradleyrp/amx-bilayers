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
