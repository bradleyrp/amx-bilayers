#!/usr/bin/env python

from amx import *

init()
if not state.force_field: state.force_field = state.before[-1]['settings']['force_field']
if not state.mdp_specs: state.mdp_specs = state.before[-1]['settings']['mdp_specs']
if not state.sources: state.sources = state.before[-1]['settings']['sources']
make_step(settings.step)
write_continue_script()
release_restraints()
write_mdp()
write_topology('system.top')
grouper(protein=protein_laden(structure='system'))
equilibrate(groups='system-groups')
