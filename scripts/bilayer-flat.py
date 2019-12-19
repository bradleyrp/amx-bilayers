#!/usr/bin/env python

from amx import *

make_step(settings.step)
write_mdp()
write_continue_script()
build_bilayer(name='vacuum-bilayer.gro')
copy_file('vacuum-bilayer.gro','vacuum.gro')
state.force_field = settings.force_field_upright
write_top('vacuum.top')
minimize('vacuum',restraints='vacuum.gro')
vacuum_pack_loop(
	structure='vacuum-minimized',
	tpr='em-vacuum-steep',
	gro='vacuum-packed')
restuff(
	structure='vacuum-packed.gro',
	gro='solvate-dry-stuff',
	tpr=get_last('tpr'),
	ndx=get_last('ndx'))
remove_jump(
	structure='solvate-dry-stuff',
	tpr=get_last('tpr'),
	gro='solvate-dry')
solvate_bilayer(
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
write_top('system.top')
equilibrate(groups='system-groups')
move_file('md-npt-bilayer.gro','md-npt-bilayer-std.gro')
distinguish_leaflets(
	structure='md-npt-bilayer-std',
	gro='system-leaflets',
	indices='dat-monolayer-indices.py')
move_file('system.top','system-std.top')
state.force_field = settings.force_field_flat
write_top('system.top')
bilayer_flatten_for_restraints(
	structure='system-leaflets',
	gro='system-leaflets-flat')
register_file('system-leaflets-flat.gro')
gmx_register_call(
	command='grompp',flag='r',
	value='system-leaflets-flat.gro')
for rnum,mdp in enumerate(state.flat_restart_mdps):
	restart_clean(part=rnum+2,mdp=mdp,
		structure='system-leaflets',
		groups='system-groups')
