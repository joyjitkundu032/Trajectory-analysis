

class PDB_reader(object):

	def __init__(self):
		return
		
	def read_pdb(pdb_input):
		i =0
		atom_name=[]
		residue_name=[]
		residue_number=[]
		chainID=[]
		atom_pos=[]
		segid=[]
		dimensions=[]
		tmp=pdb_input[0]
		#print(tmp[1:6])
		dimensions.append("".join(tmp[7:16]).strip())
		dimensions.append("".join(tmp[17:25]).strip())
		dimensions.append("".join(tmp[26:34]).strip())
		print(dimensions)
		while i < len(pdb_input)-1:
			temp = pdb_input[i]
			check_ATOM = "".join(temp[0:6]).strip()
			#print(check_ATOM)
			if check_ATOM == 'ATOM' or check_ATOM == 'HETATM':
				atom_name.append("".join(temp[13:16]).strip())
				residue_name.append("".join(temp[17:20]).strip())
				residue_number.append("".join(temp[23:26]).strip())
				chainID.append("".join(temp[21:22]))
				atom_position_x = "".join(temp[29:38]).strip()
				atom_position_y = "".join(temp[38:46]).strip()
				atom_position_z = "".join(temp[46:56]).strip()
				atom_pos.append((float(atom_position_x), float(atom_position_y), float(atom_position_z)))
				segid.append("".join(temp[72:78]).strip())
				i+=1
			else:
				i+=1
		return atom_name, residue_name, residue_number, chainID, atom_pos, segid, dimensions
