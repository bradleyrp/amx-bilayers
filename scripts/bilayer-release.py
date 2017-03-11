#!/usr/bin/env python

from amx import *

init()
make_step(settings.step)
write_continue_script()
release_restraints()
#state.mdp_specs = state.before[-1]['settings']['mdp_specs']
write_mdp()
write_topology('system.top')
restart_clean(part=1,
	structure='system',
	groups='system-groups')
