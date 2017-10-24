{
'cgmd':{
	#---defaults for redundant options (use None if only one)
	'defaults':{
		'potential':'verlet',
		'constrain':'none',
		'restrain':'none',
		'continue':'continue',
		'integrate':'fast',
		'output':'standard',
		'temperature':'none',
		'pressure':'standard',
		'groups':'none',
		'screening':'standard',
		'virtual-sites':'off',
		},
	#---standard options
	'potential':{
		'original':{
			'cutoff-scheme':'group',
			'nstlist':10,
			'ns_type':'grid',
			'pbc':'xyz',
			'rlist':'1.2',
			'coulombtype':'Shift',
			'rcoulomb_switch':0.0,
			'rcoulomb':1.2,
			'vdw_type':'Shift',
			'rvdw_switch':0.9,
			'rvdw':1.2,
			'DispCorr':'No',
			},
		'verlet':{
			'cutoff-scheme':'verlet',
			'nstlist':10,
			'ns_type':'grid',
			'pbc':'xyz',
			'coulombtype':'cut-off',
			'coulomb-modifier':'Potential-shift-verlet',
			'rcoulomb':1.2,
			'vdw_type':'cut-off',
			'vdw-modifier':'Potential-shift-verlet',
			'rvdw':1.2,
			'DispCorr':'No',
			'verlet-buffer-tolerance':0.005,
			},
		},
	'screening':{
		'standard':{'epsilon_r':15},
		'off':{'epsilon_r':0},
		},
	'virtual-sites':{'standard':{'lincs-order':6},'off':{}},
	'restrain':{
		'posre':{'define':'-DPOSRES'},
		#---"posre-com-only": keep position restraints off because we set them manually in the ITP
		#---...and we also turn on refcoord-scaling to avoid gromacs warnings with pressure coupling
		'posre-com-only':{'refcoord-scaling':'com'},
		'none':{},
		},			
	'constrain':{
		'none':{
			'constraints':'none',
			'constraint_algorithm':'Lincs',
			},
		},
	'continue':{
		'continue':{'continuation':'yes','gen_vel':'no','gen_seed':123123,},
		'start':{'continuation':'no','gen_vel':'yes','gen_seed':123123,},
		},
	'integrate':{
		'medium':{
			'integrator':'md',
			'tinit':0.0,
			'dt':0.02,
			'nsteps':50000,
			'nstcomm':100,
			},
		'fast':{
			'integrator':'md',
			'tinit':0.0,
			'dt':0.04,
			'nsteps':50000,
			'nstcomm':100,
			},
		'npt-bilayer':{
			'integrator':'md',
			'tinit':0.0,
			'dt':0.01,
			'nsteps':100000,
			'nstcomm':100,
			},
		},
	'groups':{
		'none':{
			'comm-grps':'LIPIDS SOLVENT',
			'energygrps':'LIPIDS SOLVENT',
			},
		'protein':{
			'comm-grps':'LIPIDS SOLVENT PROTEIN',
			'energygrps':'LIPIDS SOLVENT PROTEIN',
			},
		'protein-water':{
			'comm-grps':'SOLVENT PROTEIN',
			'energygrps':'SOLVENT PROTEIN',
			},
		'blank':{},
		},
	'output':{
		'standard':{
		 	'nstxout':-1,
			'nstvout':-1,
			'nstfout':0,
			'nstlog':1000,
			'nstenergy':200,
			'nstxtcout':2000,
			'xtc_precision':100,
			},
		'npt-bilayer':{
			'nstxout':-1,
			'nstvout':-1,
			'nstfout':0,
			'nstlog':100,
			'nstenergy':100,
			'nstxtcout':1000,
			'xtc_precision':100,
			},
		},
	'temperature':{
		'none':{
			'tcoupl':'v-rescale',
			'tc-grps':'LIPIDS SOLVENT',
			'tau_t':'1.0 1.0',
			'ref_t':'320 320',
			},
		'protein':{
			'tcoupl':'v-rescale',
			'tc-grps':'LIPIDS SOLVENT PROTEIN',
			'tau_t':'1.0 1.0 1.0',
			'ref_t':'320 320 320',
			},
		'protein-water':{
			'tcoupl':'v-rescale',
			'tc-grps':'SOLVENT PROTEIN',
			'tau_t':'1.0 1.0',
			'ref_t':'320 320',
			},
		'protein-vacuum':{
			'tcoupl':'v-rescale',
			'tc-grps':'system',
			'tau_t':'1.0',
			'ref_t':'320',
			},
		},
	'pressure':{
		'standard':{
			'Pcoupl':'parrinello-rahman',
			'Pcoupltype':'semiisotropic',
			'tau_p':12.0,
			'compressibility':'3e-4 3e-4',
			'ref_p':'1.0 1.0',
			},
		'standard-isotropic':{
			'Pcoupl':'parrinello-rahman',
			'Pcoupltype':'isotropic',
			'tau_p':12.0,
			'compressibility':'3e-4',
			'ref_p':'1.0',
			},
		'npt-isotropic-weak':{
			'Pcoupl':'Berendsen',
			'Pcoupltype':'isotropic',
			'tau_p':2.0,
			'compressibility':'3e-4',
			'ref_p':'1.0',
			},
		'npt-semiisotropic':{
			'Pcoupl':'Berendsen',
			'Pcoupltype':'semiisotropic',
			'tau_p':5.0,
			'compressibility':'3e-5 3e-5',
			'ref_p':'1.0 1.0',
			},
		'npt-semiisotropic-weak':{
			'Pcoupl':'Berendsen',
			'Pcoupltype':'semiisotropic',
			'tau_p':5.0,
			'compressibility':'3e-4 3e-4',
			'ref_p':'1.0 1.0',
			},
		'npt-semiisotropic-weakest':{
			'Pcoupl':'Berendsen',
			'Pcoupltype':'semiisotropic',
			'tau_p':1.0,
			'compressibility':'5e-5 5e-5',
			'ref_p':'1.0 1.0',
			},
		'nvt':{'Pcoupl':'no'}
		},
	#---protein-in-water settings
	'nvt-protein':{
		'pressure':'nvt',
		'continue':'continue',
		'integrate':'npt-bilayer',		
		'output':'npt-bilayer',
		'groups':'protein-water',
		'temperature':'protein-water',
		},
	'npt-protein':{
		'pressure':'npt-isotropic-weak',
		'continue':'continue',
		'integrate':'npt-bilayer',		
		'output':'npt-bilayer',
		'groups':'protein-water',
		'temperature':'protein-water',
		},
	'npt-protein-production':{
		'pressure':'standard-isotropic',
		'continue':'continue',
		'integrate':'npt-bilayer',		
		'output':'npt-bilayer',
		'groups':'protein-water',
		'temperature':'protein-water',
		},
	#---override for NPT bilayer and solvent
	'npt-bilayer':{
		'continue':'continue',
		'integrate':'npt-bilayer',		
		'output':'npt-bilayer',
		},
	'npt-bilayer-protein':{
		'continue':'continue',
		'integrate':'npt-bilayer',		
		'output':'npt-bilayer',
		'groups':'protein',
		'temperature':'protein',
		},
	#---override for minimization
	'minimize':{
		'groups':'blank',
		'temperature':{'tcoupl':'no'},
		'pressure':{'Pcoupl':'no'},
		'output':{
		 	'nstxout':500,
			'nstvout':500,
			'nstfout':0,
			'nstlog':500,
			'nstenergy':500,
			'nstxtcout':500,
			'xtc_precision':100,
			},
		'integrate':{
			'integrator':'steep',
			'tinit':0.0,
			'dt':0.01,
			'nsteps':100000,
			'nstcomm':1,
			'emtol':10,
			'emstep':0.01,
			},
		},		
	#---override for vacuum packing steps			
	'vacuum-packing':{
		'groups':'blank',
		'restrain':'posre',
		'screening':'off',
		'integrate':{
			'integrator':'md',
			'tinit':0.0,
			'dt':0.001,
			'nsteps':100000,
			'nstcomm':10,
			},
		'temperature':{
			'tcoupl':'Berendsen',
			'tc-grps':'SYSTEM',
			'tau_t':1.0,
			'ref_t':310,
			'gen_temp':310,
			},
		'pressure':{
			'Pcoupl':'Berendsen',
			'Pcoupltype':'semiisotropic',
			'tau_p':0.5,
			'compressibility':'3e-4 0',
			'ref_p':'1.0 1.0',
			},
		},
	#---override for vacuum packing steps			
	'invisible':{
		'groups':'blank',
		'restrain':'none',
		'potential':{
			'cutoff-scheme':'group',
			'nstlist':10,
			'ns_type':'grid',
			'pbc':'no',
			'rlist':1.0,
			'coulombtype':'Cut-off',
			'rcoulomb_switch':0.0,
			'rcoulomb':1.0,
			'vdw_type':'Cut-off',
			'rvdw':1.0,
			'DispCorr':'No',
			},
		'screening':'off',
		'integrate':{
			'integrator':'md',
			'tinit':0.0,
			'dt':0.001,
			'nsteps':100000,
			'nstcomm':10,
			},
		'temperature':{
			'tcoupl':'Berendsen',
			'tc-grps':'SYSTEM',
			'tau_t':1.0,
			'ref_t':310,
			'gen_temp':310,
			},
		'pressure':{'Pcoupl':'no'},
		},
	},
}
