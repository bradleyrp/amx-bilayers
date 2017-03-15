#!/usr/bin/env python

from amx import *

init()
make_step(settings.step)
write_continue_script()
release_restraints()
if not state.mdp_specs:
	state.mdp_specs = state.before[-1]['settings']['mdp_specs']
write_mdp()
write_topology('system.top')
grouper(protein=False)
equilibrate(groups='system-groups')
