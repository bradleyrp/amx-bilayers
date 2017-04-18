{

#---MDP settings group for atomistic CHARMM simulations
'aamd':{
#---defaults for redundant options (use None if only one)
'defaults':{
	'potential':'verlet',
	'continue':'continue',
	'temperature':'bussi',
	'pressure':'isotropic',
	'integrate':'production',
	'constrain':'hydrogen',
	'output':'long',
	'restrain':'none',
	'comm':'on',
	},
#---standard options groupings
'potential':{
	'original':{
		},
	'verlet':{
		'cutoff-scheme':'verlet',
		'nstlist':20,
		'ns_type':'grid',
		'coulombtype':'PME',
		'pme_order':4,
		'fourierspacing':0.1125,
		'rcoulomb':0.9,
		'rlist':0.9,
		'rvdw':0.9,
		'pbc':'xyz',
		'dispcorr':'EnerPres',
		},
	'bilayer':{
		'cutoff-scheme':'verlet',
		'nstlist':5,
		'ns_type':'grid',
		'coulombtype':'PME',
		'pme_order':4,
		'fourierspacing':0.1125,
		'rcoulomb':1.2,
		'rlist':1.2,
		'rvdw':1.2,
		'pbc':'xyz',
		'dispcorr':'no',
		},
	'bilayer-pack':{
		'vdw_type':'cut-off',
		'coulombtype':'cut-off',
		'nstlist':20,
		'vdw-modifier':'Potential-shift-verlet',
		'DispCorr':'No',
		'coulomb-modifier':'Potential-shift-verlet',
		'cutoff-scheme':'verlet',
		'pbc':'xyz',
		'rcoulomb':1.1,
		'ns_type':'grid',
		'rvdw':'1.1',
		},
	},
'temperature':{
	'bussi':{
		'tcoupl':'V-rescale',
		'tc_grps':'Protein Non-Protein',
		'tau_t':'0.1 0.1',
		'ref_t':'300 300',
		},
	'bussi-bilayer':{
		'tcoupl':'V-rescale',
		'tc_grps':'LIPIDS SOLVENT',
		'tau_t':'0.1 0.1',
		'ref_t':'300 300',
		},
	'bussi-bilayer-protein':{
		'tcoupl':'V-rescale',
		'tc_grps':'PROTEIN LIPIDS SOLVENT',
		'tau_t':'0.1 0.1 0.1',
		'ref_t':'300 300 300',
		},
	'off':{
		'tcoupl':'no',
		},
	},
'pressure':{
	'isotropic':{
		'pcoupl':'Parrinello-Rahman',
		'pcoupltype':'isotropic',
		'tau_p':2.0,
		'ref_p':1.0,
		'compressibility':'4.5e-5',
		},
	'isotropic-berendsen':{
		'pcoupl':'Berendsen',
		'pcoupltype':'isotropic',
		'tau_p':1.0,
		'ref_p':1.0,
		'compressibility':'4.5e-5',
		},
	'bilayer-berendsen':{
		'pcoupl':'Berendsen',
		'pcoupltype':'semiisotropic',
		'tau_p':1.0,
		'ref_p':"1.0 1.0",
		'compressibility':'4.5e-5 4.5e-5',
		},	
	'bilayer':{
		'pcoupl':'Parrinello-Rahman',
		'pcoupltype':'semiisotropic',
		'tau_p':5.0,
		'ref_p':"1.0 1.0",
		'compressibility':'4.5e-5 4.5e-5',
		},	
	'off':{
		'pcoupl':'no',
		},
	},
'constrain':{
	'hydrogen':{
		'constraints':'h-bonds',
		'constraint_algorithm':'lincs',
		},
	'all':{
		'constraints':'all-bonds',
		'constraint_algorithm':'lincs',
		},
	'off':{
		'constraints':'none',
		},
	},
'continue':{
	'continue':{'continuation':'yes'},
	'start':{'continuation':'no'},
	},
'output':{
	'long':{
		'nstxout':20000,
		'nstvout':20000,
		'nstlog':1000,
		'nstenergy':1000,
		'nstxtcout':1000,		
		},
	'short':{
		'nstxout':1000,
		'nstvout':1000,
		'nstlog':100,
		'nstenergy':100,
		'nstxtcout':500,		
		},
	'long-bilayer':{
		'nstxout':50000,
		'nstvout':50000,
		'nstlog':500,
		'nstenergy':1000,
		'nstxtcout':1000,		
		},
	},
'integrate':{
	'production':{
		'integrator':'md',
		'tinit':0,
		'dt':0.002,
		'nsteps':500000,
		'nstcalcenergy':100,
		},
	'equilibrate':{
		'integrator':'md',
		'tinit':0,
		'dt':0.001,
		'nsteps':100000,
		'nstcalcenergy':100,
		},
	'short':{
		'integrator':'md',
		'tinit':0,
		'dt':0.0001,
		'nsteps':100000,
		'nstcalcenergy':100,
		},
	'bilayer-short':{
		'integrator':'md',
		'tinit':0,
		'dt':0.0005,
		'nsteps':50000,
		'nstcalcenergy':100,
		},
	},
'comm':{
	'on':{
		'nstcomm':100,
		'comm-mode':'Linear',
		},
	'off':{
		'nstcomm':0,
		},
	},
'restrain':{
	'posre':{'define':'-DPOSRES'},
	'posre-com-only':{'refcoord-scaling':'com'},
	'none':{},
	},	
'screening':{
	'standard':{'epsilon_r':15},
	'off':{'epsilon_r':0},
	},
'groups':{
	'none':{
		'comm-grps':'LIPIDS SOLVENT',
		'energygrps':'LIPIDS SOLVENT',
		},
	'bilayer':{
		'comm-grps':'LIPIDS SOLVENT',
		'energygrps':'LIPIDS SOLVENT',
		},
	'bilayer-protein':{
		'comm-grps':'LIPIDS_PROTEIN SOLVENT',
		'energygrps':'LIPIDS SOLVENT PROTEIN',
		},
	'blank':{},
	},
#---recipes
'minimize':{
	'temperature':'off',
	'pressure':'off',
	'constrain':'off',
	'comm':'off',
	'integrate':{
		'integrator':'steep',
		'nsteps':50000,
		'emtol':10.0,
		'emstep':0.01,
		},
	},
'nvt-protein':{
	'pressure':'off',
	'temperature':'bussi',
	'output':'short',
	'continue':'start',
	'constrain':'all',
	'integrate':'equilibrate',
	},
'nvt-protein-short':{
	'continue':'start',
	'output':'short',
	'pressure':'off',
	'temperature':'bussi',
	'constrain':'hydrogen',
	'integrate':'short',
	},		
'npt-protein':{
	'continue':'start',
	'output':'short',
	'temperature':'bussi',
	'pressure':'bilayer-berendsen',
	'constrain':'all',
	'integrate':'equilibrate',
	},
'npt-bilayer':{
	'groups':'bilayer',
	'continue':'start',
	'potential':'bilayer',
	'output':'long-bilayer',		
	'temperature':'bussi-bilayer',
	'pressure':'bilayer',
	'constrain':'hydrogen',
	'integrate':'production',
	},
'npt-bilayer-simple':{
	'groups':'bilayer',
	'continue':'start',
	'potential':'bilayer',
	'output':'long-bilayer',		
	'temperature':'bussi-bilayer',
	'pressure':'bilayer-berendsen',
	'constrain':'off',
	'integrate':'bilayer-short',
	},
'npt-bilayer-protein':{
	'groups':'bilayer-protein',
	'continue':'start',
	'output':'short',
	'temperature':'bussi-bilayer-protein',
	'pressure':'bilayer-berendsen',
	'constrain':'all',
	'integrate':'equilibrate',
	},
'npt-bilayer-protein-simple':{
	'groups':'bilayer-protein',
	'continue':'start',
	'potential':'bilayer',
	'output':'long-bilayer',
	'temperature':'bussi-bilayer-protein',
	'pressure':'bilayer-berendsen',
	'constrain':'off',
	'integrate':'bilayer-short',
	},
#---override for vacuum packing steps (nearly verbatim from CGMD)
'vacuum-packing':{
	'groups':'blank',
	'restrain':'posre',
	'screening':'off',
	'potential':'bilayer-pack',
	'constrain':'hydrogen',
	'integrate':{
		'integrator':'md',
		'tinit':0.0,
		'dt':0.001,
		'nsteps':100000,
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
},
#---MDP settings group for CGMD MARTINI simulations
'cgmd':{
#---defaults for redundant options (use None if only one)
'defaults':{
	'potential':'verlet',
	'constrain':'none',
	'restrain':'none',
	'continue':'continue',
	'integrate':'medium',
	'output':'standard',
	'temperature':'none',
	'pressure':'standard',
	'groups':'none',
	'screening':'standard',
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
		'nstlist':20,
		'ns_type':'grid',
		'pbc':'xyz',
		'coulombtype':'cut-off',
		'coulomb-modifier':'Potential-shift-verlet',
		'rcoulomb':1.1,
		'vdw_type':'cut-off',
		'vdw-modifier':'Potential-shift-verlet',
		'rvdw':1.1,
		'DispCorr':'No',
		'verlet-buffer-tolerance':0.005,
		},
	},
'screening':{
	'standard':{'epsilon_r':15},
	'off':{'epsilon_r':0},
	},
'restrain':{
	'posre':{'define':'-DPOSRES'},
	'posre-pressure':{'define':'-DPOSRES','refcoord-scaling':'com'},
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
	'blank':{},
	},
'output':{
	'standard':{
	 	'nstxout':-1,
		'nstvout':-1,
		'nstfout':0,
		'nstlog':2000,
		'nstenergy':2000,
		'nstxtcout':10000,
		'xtc_precision':1000,
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
	'nvt':{'Pcoupl':'no'}
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
},		

}
