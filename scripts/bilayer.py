#!/usr/bin/env python

from amx import *

init()
make_step(settings.step)
write_mdp()
write_continue_script()
build_bilayer(name='vacuum-bilayer')
copy_file('vacuum-bilayer.gro','vacuum.gro')
state.force_field = state.force_field_upright
write_top('vacuum.top')
minimize('vacuum')
vacuum_pack_loop(
	structure='vacuum-minimized',
	tpr='em-vacuum-steep',
	gro='vacuum-packed')
restuff(
	structure='vacuum-packed',
	gro='solvate-dry',
	tpr=get_last('tpr'),
	ndx=get_last('ndx'))
solvate(
	structure='solvate-dry',
	gro='solvate')
write_top('solvate.top')
minimize('solvate')
remove_jump(
	structure='solvate-minimized',
	tpr='em-solvate-steep',
	gro='solvate-nojump')
counterions('solvate-nojump','solvate')
counterion_renamer('counterions')
write_top('counterions.top')
minimize('counterions')
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
state.force_field = settings.force_field
write_top('system.top')
equilibrate(groups='system-groups')
