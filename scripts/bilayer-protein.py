#!/usr/bin/env python

from amx import *

init()
make_step(settings.step)
write_continue_script()
write_mdp()
place_protein()
adhere_protein_bilayer('prepped-charged')
remove_ions(structure='prepped-charged',gro='prepped-water')
counterions(structure='prepped-water',gro='prepped')
counterion_renamer('prepped')
write_topology('system-minimize.top')
minimize('prepped',top='system-minimize')
detect_lipids(structure='prepped')
detect_ions(structure='prepped')
recenter_protein_bilayer(structure='prepped-minimized',gro='system')
bilayer_sorter(structure='prepped-minimized',ndx='system-groups',protein=True)
switch_itp()
write_topology('system.top')
equilibrate(groups='system-groups')
