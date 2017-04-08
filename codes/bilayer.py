#!/usr/bin/env python

"""
Bilayer construction routines.
"""

import os,re,shutil
import numpy as np
import random
from copy import deepcopy

_not_reported = ['makemesh','rotation_matrix','dotplace','distinguish_leaflets','deepcopy']
_shared_extensions = ['bilayer_sorter']

def random_lipids(total,composition,binsize):
	"""
	Generate a random 2D grid of lipids in an arrangement according to the aspect ratio.
	Note that this is currently designed for flat bilayers so the binsize spacing is on the XY plane.
	"""
	aspect = state.q('aspect')
	names,complist = zip(*composition.items())
	if sum(complist) != 1.0: complist = [float(i)/sum(complist) for i in complist]
	counts = [int(round(total*i)) for i in complist]
	arrange = np.concatenate([np.ones(j)*jj for jj,j in enumerate(counts)])
	random.shuffle(arrange)
	#---morph to fit the aspect ratio
	side_y = int(np.ceil(np.sqrt(total/aspect)))
	side_x = aspect*side_y
	pts = binsize*np.concatenate(np.array(np.meshgrid(np.arange(side_x),np.arange(side_y))).T)
	#---our grid will always be slightly too large so we trim it randomly
	selects = np.array(sorted(sorted(range(len(pts)),key=lambda *args:random.random())[:total]))
	vecs = np.array([side_x*binsize,side_y*binsize])
	return pts[selects],vecs

def makeshape():
	"""
	Generate the midplane points for various bilayer shapes.
	Historical note: comes from amx/procedures/bilayer.py
	"""
	#---this function uses settings from the state
	shape = state.q('shape')
	lz = state.q('solvent_thickness')
	binsize = state.q('binsize')
	#---are the monolayers symmetric?
	monolayer_other = state.q('monolayer_bottom',None)
	composition_other = state.q('composition_bottom',None)
	if not monolayer_other and composition_other or not composition_other and monolayer_other:
		raise Exception('you must specify both "monolayer bottom" and "composition bottom"')
	monolayers = [[state.q('monolayer_top'),state.q('composition_top')]]

	if monolayer_other: monolayers += [[monolayer_other,composition_other]]
	else: monolayers += [monolayers[0]]
	#---compute spots on the 2D grid where we will place the lipids
	#---! note that we are still in 2D so the grid-spacing is only accurate for the flat bilayer
	spots,vecs = zip(*[random_lipids(total,composition,binsize) for total,composition in monolayers])
	#---for asymmetric bilayers we choose the larger set of box vectors
	vecs = np.transpose([max(i) for i in np.transpose(vecs)])
	pts = [np.concatenate(([s[:,0]],[s[:,1]],[np.zeros(len(s))])).T for s in spots]
	#---! non-flat points will be stretched in Z hence not evenly spaced in 3-space
	#---! needs 2 monolayers
	if shape == 'saddle':

		def bump(x,y,x0,y0,height,width):
			"""
			General function for producing a 2D Gaussian "dimple" or "bump".
			"""
			zs = height*np.exp(-(x-x0)**2/2/width**2)*\
				np.exp(-(y-y0)**2/2/width**2)
			return zs

		offsets = vecs[:2]/2.-2.
		bumpspots = [(0,1,1),(0,-1,1),(1,0,-1),(-1,0,-1)]
		for spot in bumpspots:
			pts[:,2] += bump(xys[:,0],xys[:,1],
				x0+offsets[0]*spot[0],y0+offsets[1]*spot[1],
				height*spot[2],width)
		ptsmid = np.concatenate(np.reshape(pts,(3*m,3*n,3))[m:2*m,n:2*n])
		ptsmid[:,0]-=vecs[0]
		ptsmid[:,1]-=vecs[1]
	
		if 0: meshplot(ptsmid,show='surf')
	
	#---! needs 2 monolayers
	elif shape == 'buckle':

		def buckle(x,y,height):
			zs = height*np.sin(x*2*pi/lx)
			return zs
		
		pts[:,2] += buckle(xys[:,0],xys[:,1],height)

	elif shape == 'flat': pts = [p+0 for p in pts]
	else: raise Exception('\n[ERROR] unclear bilayer topography: %s'%shape)
	#---previously used PBCs and selected the middle tile here before makemesh and then shifted to origin
	monolayer_meshes = [makemesh(p,vecs,debug=False,curvilinear=False) for p in pts]
	return pts,monolayer_meshes,np.array([v for v in vecs]+[lz])

def build_bilayer(name,random_rotation=True):
	"""
	Create a new bilayer according to a particular topography.
	Historical note: comes from amx/procedures/bilayer.py
	"""
	#---collect the bilayer topography and the lipid points
	ptsmid,monolayer_mesh,vecs = makeshape()

	#---infer composition of the other monolayer
	if type(state.q('composition_top'))==str: monolayer0 = {state.q('composition_top'):1.0}
	else: monolayer0 = state.q('composition_top')
	if not state.q('composition_bottom') or not state.q('composition_bottom'): 
		monolayer1 = dict(monolayer0)
	else: monolayer1 = state.q('composition_bottom')
	nlipids0 = state.q('monolayer_top')
	if not state.q('monolayer_bottom') or not state.q('monolayer_bottom'): nlipids1 = nlipids0
	else: nlipids1 = state.q('monolayer_bottom')
	lipid_resnames = list(state.q('composition_top').keys())
	if state.q('composition_bottom'): lipid_resnames += state.q('composition_bottom').keys()
	#---save for bilayer_sorter
	state.lipids = list(set(lipid_resnames))

	lipids,lipid_order = {},[]
	lnames = list(set(list(monolayer0.keys())+list(monolayer1.keys())))
	for key in lnames:
		#---previously used read_molecule
		incoming = read_gro(os.path.join(state.q('lipid_structures'),key+'.gro'),cwd='./')
		lpts,atomnames = np.array(incoming['points']),np.array(incoming['atom_names'])
		lipids[key] = {'lpts':lpts,'atomnames':atomnames}
		lipid_order.append(key)

	#---prepare random grids
	identities = [np.zeros([nlipids0,nlipids1][mn]) for mn in range(2)]
	for mn,composition in enumerate([monolayer0,monolayer1]):
		names,complist = zip(*composition.items())
		#---allow non-unity complist sums so you can use ratios or percentages
		if sum(complist) != 1.0: complist = [float(i)/sum(complist) for i in complist]
		nlipids = [nlipids0,nlipids1][mn]
		counts = np.array([int(round(i)) for i in np.array(complist)*nlipids])
		identities[mn] = np.concatenate([np.ones(counts[ii])*lipid_order.index(lname) 
			for ii,lname in enumerate(names)])[np.array(sorted(range(nlipids),
			key=lambda *args:random.random()))]

	mono_offset = state.q('monolayer_offset')
	monolayer_indices,resnum = [[],[]],0
	#---we wish to preserve the ordering of the lipids so we write them in order of identity
	placements = [[] for l in lipid_order]
	#---loop over lipid types
	for lipid_num,lipid in enumerate(lipid_order):
		lpts = lipids[lipid]['lpts']
		#---loop over monolayers
		for mn in range(2):
			zvec = np.array([0,0,1]) if mn==0 else np.array([0,0,-1])
			#---loop over positions of this lipid type
			indices = np.where(identities[mn]==lipid_num)[0]
			for ii,index in enumerate(indices):
				status('placing %s'%lipid,i=ii,looplen=len(indices))
				#---begin move routine
				point = ptsmid[mn][index]
				offset = [1,-1][mn]*mono_offset*monolayer_mesh[mn]['vertnorms'][index]
				if random_rotation:
					random_angle = np.random.uniform()*2*np.pi
					lpts_copy = np.dot(rotation_matrix(zvec,random_angle),lpts.T).T
				else: lpts_copy = np.array(lpts)
				xys = point+offset+[1,-1][mn]*np.dot(
					rotation_matrix(np.cross(zvec,monolayer_mesh[mn]['vertnorms'][index]),
					np.dot(zvec,monolayer_mesh[mn]['vertnorms'][index])),lpts_copy.T).T
				#---end move routine
				placements[lipid_num].append(xys)
				monolayer_indices[mn].append(resnum)
				resnum += 1
	natoms = np.sum([np.sum(np.concatenate(identities)==ii)*len(lipids[i]['lpts']) 
		for ii,i in enumerate(lipid_order)])

	#---enforce a minimum z-height in case the solvent thickness is very small
	#---! note that the solvent thickness might be incorrectly applied here
	#---! it might be better to use a fixed thickness for building the planar bilayer anyway
	vecs[2] = 10.0 if not vecs[2] or vecs[2]<10.0 else vecs[2]

	#---output the monolayer indices for posterity
	with open(state.here+'dat-monolayer-indices.py','w') as fp: fp.write(str(monolayer_indices))

	#---write the placed lipids to a file
	resnr = 1
	with open(state.here+name+'.gro','w') as fp:
		fp.write('%s\n'%state.q('system_name')+'%d\n'%natoms)
		#---loop over lipid types
		for lipid_num,resname in enumerate(lipid_order):
			atomnames = lipids[resname]['atomnames']
			for xys in placements[lipid_num]:
				fp.write('\n'.join([''.join([
					str(resnr).rjust(5),
					resname.ljust(5),
					atomnames[i].rjust(5),
					(str((resnr-1)*len(xys)+i+1).rjust(5))[:5],
					''.join([dotplace(x) for x in xys[i]])])
					for i in range(len(xys))])+'\n')
				resnr += 1
		fp.write(' '.join([dotplace(x) for x in vecs])+'\n')

	#---save the slab dimensions for solvation
	boxdims_old,boxdims = get_box_vectors(name)
	state.bilayer_dimensions_slab = boxdims
	#---this function is *always* a start point even in the middle of a test set so we reset composition
	if 'composition' in state: del state['composition']
	#---save composition for topology
	for lipid_num,lipid in enumerate(lipid_order):
		component(lipid,count=len(placements[lipid_num]))

def lipid_upright():
	"""
	Assuming that incoming ITP files do not have position restraints, we add them here.
	These position restraints ensure an upright lipid in a vacuum-packed bilayer.
	"""
	#---select the ITP file containing all lipids
	if not state.lipids_itp: raise Exception('to add upright lipid contraints, we need a lipids_itp to rewrite')
	itps = GMXTopology(os.path.join(state.here,state.force_field+'.ff',state.lipids_itp))
	import ipdb;ipdb.set_trace()
	sys.exit(1)

def distinguish_leaflets(structure='incoming',gro='outgoing',indices='dat-monolayer-indices.py'):
	"""
	Take a gro file and some knowledge about leaflets and distinguish them.
	"""
	#---monolayer indices are saved earlier
	monolayer_inds = eval(open(state.here+indices).read())
	struct = GMXStructure(state.here+structure+'.gro')
	not_lipids = [state.sol,'ION',state.anion,state.cation]
	resnames = [i for i in np.unique(np.array(struct.residue_names,dtype='|S5')) 
		if i.decode() not in not_lipids]
	#---only label the bottom leaflet because that one will start with restraints
	renamer = dict([(i,{'top':i+b'','bot':i+b'R'}) for i in resnames])
	#---rename residues by leaflet
	leaflet_index_convetion = {'top':0,'bot':1}
	residue_indices = np.array(struct.residue_indices)
	residue_names = np.array(struct.residue_names,dtype='|S5')
	for ii,i in enumerate(np.unique(struct.residue_indices)):
		for which in ['top','bot']:
			if ii in monolayer_inds[leaflet_index_convetion[which]]:
				target_resids = np.where(residue_indices==i)
				target_resname = residue_names[target_resids[0][0]]
				residue_names[target_resids] = renamer[target_resname][which].decode()
	struct.residue_names = residue_names
	struct.write(state.here+gro+'.gro')

	if False:
		#---generate new topologies for the restrained lipids
		#---! remove path below
		itp = GMXTopology('s01-bilayer/martini.ff/martini-v2.0-lipids.itp')
		restrain_lipids = [('DOPC','DOPCR'),('DOPS','DOPSR')]
		rl = restrain_lipids[0]

		itp.molecules = dict([(k,v) for k,v in itp.molecules.items() if k in ['DOPC','DOPS']])
		
		#---ensure no position restraints on the regular versions
		#---! save them because no automatic writing yet
		posres = {}
		for mol in itp.molecules.keys():  
			posres[mol] = itp.molecules[mol].pop('position_restraints')
		#---! intervene because we only want to restration atoms 1,9,14
		restrain_res = [1,9,14]
		posres_custom = [{'funct': '1', 'fcy': '0', 'ai': '1', 'fcx': '0', 'fcz': '1000' 
			if ii+1 in restrain_res else '0'} for ii in range(14)]
		#---make restrain versions
		for mol in list(itp.molecules.keys()):
			itp.molecules[mol+'R'] = deepcopy(itp.molecules[mol])
			itp.molecules[mol+'R']['position_restraints'] = posres_custom
			#---change the name of a residue
			itp.molecules[mol+'R']['moleculetype']['molname'] = mol+'R'
			for dd,d in enumerate(itp.molecules[mol+'R']['atoms']):
				itp.molecules[mol+'R']['atoms'][dd]['resname'] = mol+'R'

		#---! hacked
		new_itp = os.path.join(state.here+state.force_field+'.ff','martini-v2.0-lipids-flat-top.itp')
		try: os.remove(new_itp)
		except: pass
		#---write the flat-top ITP file
		itp.write(new_itp)

		#---! hack the force field includes
		state.ff_includes = ['martini-v2.2','martini-v2.0-lipids-flat-top',
			'martini-v2.2-aminoacids','martini-v2.0-ions','PIP2']

	#---GMXStructure already knows how to get the right composition
	state.composition = struct.detect_composition()
	state.lipids = [i for i in list(zip(*state.composition))[0] if i not in [state.sol,state.anion,state.cation]]
	#---force field swapping happens in the parent script

def remove_jump(structure,tpr,gro,pbc='nojump'):
	"""
	Correct that thing where the bilayer crosses the PBCs and gets split.
	"""
	gmx('make_ndx',ndx=structure,structure=structure,inpipe="keep 0\nq\n",log='make-ndx-%s'%pbc)	
	gmx('trjconv',ndx=structure,structure=structure,gro=gro,tpr=tpr,
		log='trjconv-%s-%s'%(structure,pbc),pbc=pbc)
	os.remove(state.here+'log-'+'make-ndx-%s'%pbc)

def vacuum_pack(structure='vacuum',name='vacuum-pack',gro='vacuum-packed',pbc='nojump'):
	"""
	Pack the lipids in the plane, gently.
	"""
	gmx('grompp',base='md-%s'%name,top='vacuum',
		structure=structure,log='grompp-%s'%name,mdp='input-md-%s-eq-in'%name,
		maxwarn=100)
	gmx('mdrun',base='md-%s'%name,log='mdrun-%s'%name,nonessential=True)
	if pbc:
		remove_jump(structure='md-%s'%name,tpr='md-'+name,gro='md-%s-%s'%(name,pbc))
		copy_file('md-%s-%s.gro'%(name,pbc),'%s.gro'%gro)
	else: copy_file('md-%s'%gro,'%s.gro'%gro)
	boxdims_old,boxdims = get_box_vectors(gro)
	state.bilayer_dimensions_slab[:2] = boxdims_old[:2]

def vacuum_pack_loop(structure,gro,tpr):
	"""
	Run multiple vacuum packing loops.
	"""
	#---get vacuum-packing MDP specs
	valid_mdps = filter(lambda x:re.match('^input-md-vacuum',x),state.mdp_specs['mdps'])
	#---make step names
	mdps = dict([(re.match('^input-md-(.*?)-eq-in.mdp',x).group(1),x) for x in valid_mdps])
	names = sorted(mdps.keys())
	struct_in = 'vacuum-nojump'
	remove_jump(structure=structure,tpr=tpr,gro=struct_in)
	for nn,name in enumerate(names):
		mdp = mdps[name]
		vacuum_pack(structure=struct_in,name=name,gro=name)
		struct_in = name
	copy_file(name+'.gro',gro+'.gro')

def solvate_bilayer(structure='vacuum'):
	"""
	Solvate a CGMD bilayer (possibly with proteins) avoiding overlaps.
	"""
	#---check the size of the slab
	incoming_structure = str(structure)
	boxdims_old,boxdims = get_box_vectors(structure)
	#---check the size of the water box
	waterbox = state.water_box
	basedim,_ = get_box_vectors(waterbox)
	if not all([i==basedim[0] for i in basedim]):
		raise Exception('[ERROR] expecting water box to be cubic but boxdims are %s'%str(basedim))
	else: basedim = basedim[0]
	#---make an oversized water box
	newdims = boxdims_old[:2]+[state.solvent_thickness]
	gmx('genconf',structure=waterbox,gro='solvate-empty-uncentered-untrimmed',
		nbox=' '.join([str(int(i/basedim+1)) for i in newdims]),log='genconf')
	#---trim the blank water box
	trim_waters(structure='solvate-empty-uncentered-untrimmed',
		gro='solvate-empty-uncentered',boxcut=True,boxvecs=newdims,
		gap=0.0,method=state.atom_resolution)
	#---update waters
	structure='solvate-empty-uncentered'
	component(state.sol,count=count_molecules(structure,state.sol))
	#---translate the water box
	gmx('editconf',structure=structure,gro='solvate-water-shifted',
		translate='0 0 %f'%(state.bilayer_dimensions_slab[2]/2.),log='editconf-solvate-shift')
	#---combine and trim with new box vectors
	structure = 'solvate-water-shifted'
	boxdims_old,boxdims = get_box_vectors(structure)
	boxvecs = state.bilayer_dimensions_slab[:2]+[state.bilayer_dimensions_slab[2]+boxdims[2]]
	gro_combinator('%s.gro'%incoming_structure,structure,box=boxvecs,
		cwd=state.here,gro='solvate-dense')
	structure = 'solvate-dense'
	#---trim everything so that waters are positioned in the box without steric clashes
	trim_waters(structure=structure,gro='solvate',boxcut=False,
		gap=state.protein_water_gap,method=state.atom_resolution,boxvecs=boxvecs)
	structure = 'solvate'
	nwaters = count_molecules(structure,state.sol)/({'aamd':3.0,'cgmd':1.0}[state.atom_resolution])
	if round(nwaters)!=nwaters: raise Exception('[ERROR] fractional water molecules')
	else: nwaters = int(nwaters)
	component(state.sol,count=nwaters)
	state.bilayer_dimensions_solvate = boxvecs
	state.water_without_ions = nwaters

def counterion_renamer(structure):
	"""
	Fix the ion names for MARTINI.
	"""
	with open(state.here+structure+'.gro') as fp: lines = fp.readlines()
	for lineno,line in enumerate(lines):
		if re.match('.{5}(CL|NA)',line):
			lines[lineno] = re.sub(re.escape(line[5:15]),
				'ION  '+line[5:10].strip().rjust(5),lines[lineno])
	with open(state.here+structure+'.gro','w') as fp:
		for line in lines: fp.write(line)

def bilayer_middle(structure,gro):
	"""
	Move the bilayer to the middle of the z-coordinate of the box.
	Note that the protein adhesion procedure works best on a slab that is centered on z=0.
	This means that the bilayer will be broken across z=0.
	For visualization it is better to center it.
	"""
	gmx('make_ndx',ndx='system-dry',structure='counterions-minimized',
		inpipe="keep 0\nr %s || r ION || r %s || r %s\n!1\ndel 1\nq\n"%(
		state.sol,state.anion,state.cation),
		log='make-ndx-center')
	#---bilayer slab is near z=0 so it is likely split so we shift by half of the box vector
	gmx('trjconv',structure='counterions-minimized',gro='counterions-shifted',ndx='system-dry',
		trans='0 0 %f'%(state.bilayer_dimensions_solvate[2]/2.),pbc='mol',
		tpr='em-counterions-steep',log='trjconv-shift',inpipe="0\n0\n")
	#---! added extra 0 above for shifting -- not sure why this error wasn't caught earlier?
	#---center everything
	gmx('trjconv',structure='counterions-shifted',gro='system',ndx='system-dry',
		tpr='em-counterions-steep',log='trjconv-middle',inpipe="1\n0\n",center=True,pbc='mol')

def bilayer_sorter(structure,ndx='system-groups',protein=False):
	"""
	Divide the system into groups.
	"""
	#---figure out if we are aamd or cgmd
	sol_list = [state.sol,'ION',state.cation,state.anion]
	#---! check for lipids and ions !
	#---! fix the protein_ready reference here
	if 'protein_ready' in state or protein:
		#---in the atomistic method, the make_ndx output will figure out which items are proteins
		if atomistic_or_coarse()=='aamd':
			gmx('make_ndx',structure=structure,ndx='%s-inspect'%structure,
				log='make-ndx-%s-inspect'%structure,inpipe="q\n")
			with open(state.here+'log-make-ndx-%s-inspect'%structure) as fp: lines = fp.readlines()
			#---find the protein group because it may not be obvious in CGMD
			make_ndx_sifter = '^\s*([0-9]+)\s*Protein'
			protein_group = int(re.findall(make_ndx_sifter,
				next(i for i in lines if re.match(make_ndx_sifter,i)))[0])
			group_selector = "\n".join([
				"keep %s"%protein_group,
				"name 0 PROTEIN",
				" || ".join(['r '+r for r in state.lipids]),
				"name 1 LIPIDS",
				" || ".join(['r '+r for r in sol_list]),
				"name 2 SOLVENT",
				"0 | 1 | 2","name 3 SYSTEM","0 | 1","name 4 LIPIDS_PROTEIN","q"])+"\n"
		else:
			#---! hard-coded for MARTINI
			land = Landscape('martini')
			protein_selection = land.protein_selection()
			group_selector = "\n".join([
				"keep 0",protein_selection,"keep 1","name 0 PROTEIN",
				" || ".join(['r '+r for r in state.lipids]),
				"name 1 LIPIDS",
				" || ".join(['r '+r for r in sol_list]),
				"name 2 SOLVENT",
				"0 | 1 | 2","name 3 SYSTEM","0 | 1","name 4 LIPIDS_PROTEIN","q"])+"\n"
	else:
		group_selector = "\n".join([
			"keep 0",
			"name 0 SYSTEM",
			" || ".join(['r '+r for r in state.lipids]),
			"name 1 LIPIDS",
			" || ".join(['r '+r for r in sol_list]),
			"name 2 SOLVENT","q"])+"\n"
	gmx('make_ndx',structure='system',ndx=ndx,log='make-ndx-groups',
		inpipe=group_selector)

def bilayer_flatten_for_restraints(structure,gro):
	"""
	Read a bilayer and replace its coordinates with perfectly flat ones.
	"""
	#---read the gro
	struct = GMXStructure(state.here+structure+'.gro')	
	#---load the itps 
	ff = GMXForceField(state.here+state.force_field+'.ff')
	#---find all molecules with z-restrained atoms and log them in a rules list
	restraint_rules = {}
	for name in ff.molecules_list():
		mol = ff.molecule(name)
		if 'position_restraints' in mol:
			atoms = [mol['atoms'][ii]['atom'] for ii,i in 
				enumerate(mol['position_restraints']) if float(i['fcz'])>0]
			restraint_rules[name] = atoms
	#---loop over all residue names in the restraint list
	for resname in [i for i in np.unique(struct.residue_names) if i in restraint_rules]:
		for atom_name in restraint_rules[resname]:
			#---! note that this will only align atoms by resname,type
			#---get the average position of that resname, atom_type 
			subjects = np.where(np.all((
				struct.residue_names==resname,
				struct.atom_names==atom_name),axis=0))[0]
			mean_z = struct.points[subjects].mean(axis=0)[2]
			struct.points[subjects,2] = mean_z
	struct.write(state.here+gro+'.gro')

def release_restraints():
	"""
	Custom method for removing restraints.
	"""
	#---replace restrained lipid names with the correct topology
	lipid_regex = '^([A-Z]{4})R$'
	renames = []
	for key in state.lipids:
		if re.match(lipid_regex,key):
			index = list(zip(*state.composition))[0].index(key)
			new_name = re.match(lipid_regex,key).group(1)
			state.composition[index][0] = new_name
			renames.append([key,new_name])
	#---copy ITP files and system-groups.ndx
	#---! it might be worth standardizing this step
	#---! need a systematic way to track groups
	for itp in (state.itp or [])+['system-groups.ndx']:
		shutil.copyfile(state.before[-1]['here']+itp,state.here+itp)
	#---get the last frame
	get_last_frame(gro='system-restrained')
	struct = GMXStructure(state.here+'system-restrained.gro')
	for old,new in renames: 
		#---rename residues
		struct.residue_names[np.where(struct.residue_names==old)] = new
		#---rename lipids
		state.lipids[state.lipids.index(old)] = new
	struct.write(state.here+'system.gro')
	#---update the force field according to the settings
	state.force_field = settings.force_field
