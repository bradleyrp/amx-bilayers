#!/usr/bin/env python

import os,shutil,glob,re
import numpy as np
from runner import DotDict
from copy import deepcopy
import scipy

_not_reported = ['lay_coords_flat','make_list']
_not_all = ['deepcopy']
_shared_extensions = ['lay_coords_flat']

make_list = lambda x : x if type(x)==list else [x]

def place_protein():
	"""
	Choose one of the protein placement options partitioned into functions below.
	"""
	#---in test sets the force field might change so we always pick it up from settings
	if not settings.force_field: raise Exception('you must set the force field in the settings block')
	state.force_field = settings.force_field
	if not state.placement_method: raise Exception('need a `placement method` in the settings')
	if state.placement_method == 'globular_up_down': place_protein_globular_up_down()
	elif state.placement_method == 'banana': place_protein_banana()
	else: raise Exception('undeveloped placement method: %s'%state.placement_method)

def place_protein_globular_up_down():
	"""
	Fix the name please!
	"""
	#---for the globular_up_down method, the reference axis is the one the up/down axis is aligned to
	reference_axis = state.q('reference_axis',[0,0,1])

	#---get the protein structure from a dedicated variable
	if not state.protein_prepared: raise Exception('this step requires state.protein_prepared')
	protein_fn = state.protein_prepared['gro']
	protein = GMXStructure(protein_fn)
	
	#---find the down and up groups
	pts_down = protein.select_center(state.group_down)
	pts_up = protein.select_center(state.group_up)
	#---translate the down group to the origin
	protein.points -= pts_down
	#---identify the axis between the up and down groups
	axis = vecnorm(pts_up-pts_down)
	#---identify the orthogonal axis to the protein axis and the reference axis
	refaxis = np.array(reference_axis)
	orthaxis = vecnorm(np.cross(refaxis,axis))
	#---compute angle between reference axis and protein axis and the resulting rotation 
	angle = np.arccos(np.dot(refaxis,axis))
	rotation = rotation_matrix(orthaxis,angle)
	#---apply the rotation
	protein.points = np.dot(protein.points,rotation)
	#---make sure the box is big enough if we want to check it in e.g. VMD
	protein.box = protein.points.ptp(axis=0)-protein.points.min(axis=0)+1.0
	protein.write(state.here+'protein-placed.gro')

def lay_coords_flat(points,direction='y'):
	"""
	Rotate a protein so its first principal axis is aligned with one of the cartesian directions.
	"""
	pts = np.array(points)
	pts -= np.mean(pts,axis=0)
	#---take the first principal component as the axis
	axis = vecnorm(pts[0]-pts[-1])
	eigs = np.linalg.eig(np.dot(pts.T,pts))
	principal_axis_index = np.argsort(eigs[0])[-1]
	axis = vecnorm(eigs[1][:,principal_axis_index])
	refaxis = np.zeros(3)
	refaxis['xyz'.index(direction)] = 1
	xyzs = np.dot(pts,rotation_matrix(vecnorm(np.cross(refaxis,axis)),
		np.arccos(np.dot(refaxis,axis)))) 
	return xyzs

def place_protein_banana():
	"""
	Protein flat-laying banana procedure.
	In the banana procedure we start with a protein with an obvious first principal axis 
	(i.e. an elogated rod of some kind). We align the first principal axis with a cartesian direction.
	Then we select an "up" group and a "down" group. The vector between groups (hereafter: axis) is 
	projected onto the plane normal to the direction vector (now the first principal axis of the protein).
	This projected axis is used to describe a rotation. The final rotation forces two constraints on the 
	protein: its longest axis is aligned in a particular direction, and the vector between up and down groups
	is parallel to the negative z-vector. This is equivalent to enforcing an orthogonal basis for the protein,
	where the second vector is a projection of the up/down axis.
	"""
	#---! hard-coded directions for now, but these can be options, or generalized in a lattice-maker
	direction,direction_down = state.protein_lattice.get('banana_direction','y'),'z'
	ref_axis,down_axis = np.zeros(3),np.zeros(3)
	ref_axis['xyz'.index(direction)] = 1
	down_axis['xyz'.index(direction_down)] = -1
	#---get the protein structure from a dedicated variable
	if not state.protein_prepared: raise Exception('this step requires state.protein_prepared')
	protein_fn = state.protein_prepared['gro']
	#---lay the protein flat along the direction
	protein = GMXStructure(protein_fn)
	protein.points = lay_coords_flat(protein.points,direction=direction)
	downer = protein.select_center(state.group_down)
	centroid = protein.cog(protein.select('all'))
	#---note that previous version of the banana routine mistakenly omitted the mean-centering
	coords -= centroid
	coords = np.array(protein.points)
	#---project the vector between centroid and downward-facing group onto the plane normal to the direction
	axis = vecnorm(downer-centroid)
	principal = vecnorm(principal_axis(coords))
	projected = vecnorm(plane_project(axis,principal))
	#---rotate protein along the direction so the axis points down
	rotation_angle = vecangle(down_axis,projected)
	rotation = rotation_matrix(principal,-1*rotation_angle)
	#---apply the rotation
	coords_rotated = np.dot(coords,rotation)
	#---shift back to the original centroid now that rotation is complete
	coords_rotated += centroid
	#---write the modified structure
	protein.points = coords_rotated
	#---require an origin group to act as the reference point for the protein
	group_origin = state.q('group_origin')
	if not group_origin: raise Exception('banana requires group_origin or group_down in the settings')
	pts_origin = protein.select_center(group_origin)
	protein.points -= pts_origin
	#---make sure the box is big enough if we want to check it in e.g. VMD
	protein.box = protein.points.ptp(axis=0)-protein.points.min(axis=0)+1.0
	protein.write(state.here+'protein-placed.gro')
	#---! still need to apply "center_selection" (lipid pocket) and "over"
	#---! also need to recenter the protein laterally (see original code)

def adhere_protein_bilayer(gro,debug=False,**kwargs):
	"""
	MOVING FROM bilayer.py !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	adhere_protein_cgmd_bilayer(bilayer,combo,protein_complex=None)
	Attach proteins to a CGMD (?) bilayer.
	HEAVILY MODIFIED FROM ORIGINAL AUTOMACS.
	Removed debugging, pocket lipids, etc.
	"""
	#---INPUTS
	lattice = DotDict(**state.protein_lattice)
	ncols,nrows = lattice.ncols,lattice.nrows
	lattice_type = lattice.lattice_type
	space_scale = lattice.space_scale
	total_proteins = lattice.total_proteins
	z_shift = lattice.protein_shift_up

	#---create the mesh of focal points
	grid = [(i%ncols,int(i/nrows)) for i in np.arange(nrows*ncols)]
	if lattice_type == 'square': 
		vert = horz = space_scale
		offset = 0
		while len(grid)>total_proteins: grid.remove(-1)
	elif lattice_type == 'triangle': 
		horz,vert = space_scale,space_scale*sqrt(2.0)/2
		offset = space_scale/2.
		#---if total proteins is lower than the grid we remove from the end of odd rows
		scan = 1
		while len(grid)>total_proteins:
			grid.remove((ncols-1,scan))
			scan += 2
	else: raise Exception('unclear lattice type: %s'%lattice_type)
	grid_space = np.array([(horz*i+j%2.0*offset,vert*j) for i,j in np.array(grid).astype(float)])
	focii = np.concatenate((grid_space.T,[np.zeros(total_proteins)])).T

	#---load the bilayer and protein structures
	bilayer = GMXStructure(state.bilayer_structure)
	protein = GMXStructure(state.here+'protein-placed.gro')
	protein_copy = deepcopy(protein)

	#---composition starts with the free bilayer, and later we add proteins to it
	state.composition = bilayer.detect_composition()

	#---center the lattice in the middle of the XY plane of the box (z_shift is made redundant below)
	center_shift = np.array([j/2. for j in bilayer.box[:2]]+
		[z_shift])-np.concatenate((np.mean(grid_space,axis=0),[0]))

	#---loop over protein placement positions
	for translate in grid_space[::-1]:
		#---we only shift in XY
		np.concatenate((translate,[0]))
		#---use a simple distance finder to determine the z-offset for the desired location in the XY
		drop_spot = np.concatenate((translate,[0])) + center_shift
		#---to avoid the problem of distinguishing monolayers we assume a flat bilayer and drop from box top
		drop_spot[2] = bilayer.box[2]
		lipid_indices = bilayer.select('lipid')
		#---perform search in XYZ but if the bilayer is flat it is as if we searched in XY
		closest = scipy.spatial.KDTree(bilayer.points[lipid_indices]).query(drop_spot)[1]
		#---find the nearest lipid residue number
		closest_resnum = bilayer.residue_indices[closest]
		#---pivot is the height of the top atom for the nearest lipid at the drop_spot in XY
		#---...which is a good minimal definition for the nearest part of the surface
		#---...but probably not a good candidate for which lipid to remove when we do replacements later
		pivot = np.array(drop_spot)
		pivot[2] = bilayer.points[np.where(bilayer.residue_indices==closest_resnum)[0]][:,2].max()
		#---when using the simple distance finder we apply the shift here
		relative_origin = np.concatenate((translate,[0])) + pivot + np.array([0,0,z_shift])
		protein_copy.points = protein.points + relative_origin
		bilayer.add(protein_copy,before=True)

	#---write the full system before running trim_waters
	scale_method = atomistic_or_coarse()
	bilayer.write(state.here+'combo-untrimmed.gro')
	trim_waters(structure='combo-untrimmed',gro=gro,gap=state.protein_water_gap,
		method=scale_method,boxvecs=bilayer.box)

	#---! when to reionize ??
	
	#---after trimming waters we have to fix the composition
	combo = GMXStructure(state.here+gro+'.gro')
	#---! note that detect composition only gets the waters/lipids right so we keep track of proteins
	n_waters = [j for i,j in combo.detect_composition() if i==state.sol]
	if len(n_waters)!=1: raise Exception('incorrect number of water counts from detect_composition')
	component(state.sol,count=n_waters[0])
	if scale_method=='aamd':
		#---we expect that protein.itp is already present in the itp list 
		#---lipid ITP might be in a separate place, so we get this directly from the landscape
		#---the ITP must be copied in via sources or files in the settings
		#---! later the landscape can generate these ITP files automatically
		if 'protein_prepared' in state:
			for item in ['itp','posre']:
				if item in state.protein_prepared:
					shutil.copyfile(state.protein_prepared[item],
						state.here+os.path.basename(state.protein_prepared[item]))
		for fn in [os.path.relpath(i,state.here) for i in glob.glob(state.here+state.extra_itps)]:
			state.itp.append(fn)
	elif scale_method=='cgmd':
		#---handle martini ITP extraction
		#---! note that we already use "protein_prepared" to get the gro file so perhaps we could use it for
		#---! ...topology as well and eliminated the aamd/cgmd distinction here
		if state.martinize_itps:
			if not state.itp: state.itp = []
			for fn in state.martinize_itps: 
				shutil.copyfile(fn,state.here+os.path.basename(fn))
				state.itp.append(os.path.basename(fn))
	else: raise Exception('unclear scale: %s'%scale_method)



	#---! only works for a single incoming protein type and corresponding ITP
	collected_protein_itps = [GMXTopology(state.here+fn) for fn in state.itp]
	molecules = dict([j for k in [i.molecules.items() for i in collected_protein_itps] for j in k])
	if not state.q('two_protein_hack',False):
		"""
		in martini we typically have the lipids in the ff and a single incoming protein.itp
		however in aamd we may have lipid itp files as well. the construction procedure always places proteins
		first, so in the event that we have multiple ITP files in state.itp we fish out the protein one and place
		it first in line in the composition and then just hope for the best
		"""
		if len(molecules)>1:
			molecule_names_protein = [i for i in molecules.keys() if re.search('(P|p)rotein',i)]
			if len(molecule_names_protein)!=1: 
				raise Exception('we need to fish out only a single protein from the molecule list but we got: %s'%
					molecule_names_protein)
			protein_molecule_name = molecule_names_protein[0]
		else: protein_molecule_name = list(molecules.keys())[0]
		state.composition = [[protein_molecule_name,total_proteins]] + state.composition
		land = Landscape()
		state.lipids = [i for i in list(zip(*state.composition))[0] if i in Landscape().lipids()]
	else:
		#---from PT a slightly more elegant hack
		molecule_names_protein = [i for i in molecules.keys() if re.search('(P|p)rotein',i)]
		#if len(molecules)>1:
		#	molecule_names_protein = [i for i in molecules.keys() if re.search('(P|p)rotein',i)]
		#	if len(molecule_names_protein)!=1: 
		#		raise Exception('we need to fish out only a single protein from the molecule list but we got: %s'%
		#			molecule_names_protein)
		#	protein_molecule_name = molecule_names_protein[0]
		#else: protein_molecule_name = list(molecules.keys())[0]
		protein_molecule_name = molecule_names_protein[0]
		#---rpb sets total_proteins below to get the PT hack to work
		total_proteins = 2
		state.composition = [[protein_molecule_name,total_proteins]] + state.composition
		land = Landscape()
		state.lipids = [i for i in list(zip(*state.composition))[0] if i in Landscape().lipids()]

def recenter_protein_bilayer(structure,gro):
	"""
	Run this after any grompp step to put everything back in the box.
	! repetitive with center_bilayer in bilayer.py used on the standard method
	! depends on state.lipids being already identified ... 
	"""
	subjects = state.lipids + Landscape.protein_residues
	selection = ' or '.join(['r %s'%i for i in subjects])
	center_by_group(structure=structure,gro=gro,selection=selection)

def detect_lipids(structure):
	"""
	"""
	#---! curently detect_Lipids and detect ions is only necessary for bilayer sorter -- move it there??
	#---! hard-coded martini landscape
	land = Landscape()
	bilayer = GMXStructure(state.here+structure+'.gro')
	lipids = land.lipids() + land.sterols()
	#---! hack for restraint-named lipids
	lipids += ['%sR'%i for i in lipids if len(i)==4]
	state.lipids = [i for i in lipids if i in np.unique(bilayer.residue_names)]

def detect_ions(structure):
	"""
	"""
	#---! hard-coded martini landscape
	land = Landscape()
	bilayer = GMXStructure(state.here+structure+'.gro')
	#---! identifying ions should be left to the Landscape itself since it can decide on NA+ ION vs NA NA
	cations = [i for i in np.unique(bilayer.atom_names) if i in land.cations()]
	anions = [i for i in np.unique(bilayer.atom_names) if i in land.anions()]
	if len(cations)!=1: raise Exception('cannot identify a single cation: %s'%cations)
	if len(anions)!=1: raise Exception('cannot identify a single anion: %s'%anions)
	state.cation,state.anion = cations[0],anions[0]

###---BELOW IS TRASH

def detect_lipids00000(structure):
	"""
	Referesh state.lipids for other functions e.g. bilayer_sorter.
	"""
	bilayer = GMXStructure(state.here+structure+'.gro')
	land = bilayer.get_landscape()
	#---! repetitive with recenter_protein_bilayer above
	any_lipid_resnames = [k for k,v in land['objects'].items() if 'is' in v and v['is']=='lipid']
	#---! hack to include restrained lipids
	any_lipid_resnames += ['%sR'%l for l in any_lipid_resnames]
	state.lipids = [i for i in any_lipid_resnames if i in np.unique(bilayer.residue_names)]

def detect_ions0000(structure):
	"""
	Referesh state.cation and state.anion for other functions e.g. bilayer_sorter.
	Insists on only one cation and one anion.
	"""
	bilayer = GMXStructure(state.here+structure+'.gro')
	land = bilayer.get_landscape()
	any_ions = [k for k,v in land['objects'].items() if 'is' in v and v['is']=='ion']
	#---! only works for martini right now
	cation = [i for i in np.unique(bilayer.atom_names) if i in any_ions and land['objects'][i]['charge']>0]
	anion = [i for i in np.unique(bilayer.atom_names) if i in any_ions and land['objects'][i]['charge']<0]
	if len(cation)!=1: raise Exception('cannot identify a single cation: %s'%cation)
	if len(anion)!=1: raise Exception('cannot identify a single anion: %s'%anion)
	state.cation,state.anion = cation[0],anion[0]

def remove_ions(structure,gro):
	"""
	HACKISH WAY TO RE-NEUTRALIZE.
	BETTER OPTION MIGHT JUST BE TO RESOLVATE ENTIRELY ???
	#---! handle "reionize method" option ^^^
	"""
	#---! previous method was deleting NC3 from POPC because it's also an ion...
	struct = GMXStructure(state.here+structure+'.gro')
	ion_inds = struct.select('ion')
	struct.remove(ion_inds)
	struct.write(state.here+gro+'.gro')
	struct.detect_composition()
	state.water_without_ions = component(state.sol)
	#---! hard-coded for cgmd
	land = Landscape()
	#---remove ions from the composition by the molecule name returned as a key in Landscape.objects
	for i in land.my(struct,'anions')+land.my(struct,'cations'): 
		if i in list(zip(*state.composition))[0]: component(i,count=0)
