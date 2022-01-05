from Bio.PDB import *
from glob import glob
from datetime import datetime
import sys
import array as arr
import pandas as pd
import numpy as np
import os

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

####
#
#       This script will mine the ligand binding pocket of local copies of pdb files, it will count residues are present
#       within n angstoms of any ligand atom and output the results to a csv file
#
####

ligand_id = "ATP"

#~~~~~~~Specify the location of the pdbfiles to be mined


pdb_files = glob(os.path.expanduser("~/Desktop/pdb_files/"+ligand_id+"/*"))

#~~~~~~~ligandc is a counter for how many ligands are present in a structure
ligandc = 0

#~~~~~~~The radius to be searched within
aradius = arr.array('d', [3.5])

AAc = 0
AA = ['ALA',	'ARG',	'ASN',	'ASP',	'CYS',	'GLN',	'GLU',	'GLY',	'HIS',	'ILE',	'LEU',	'LYS',	'MET',	'PHE',	'PRO',	'SER',	'THR',	'TRP',	'TYR',	'VAL']

#~~~~~~~The list must contain the the ligand atoms - generally hydrogens are ignored so they are not on this list.
ligand_atoms = ["PG",	"O1G",	"O2G",	"O3G",	"PB",	"O1B",	"O2B",	"O3B",	"PA",	"O1A",	"O2A",	"O3A",	"O5'",	"C5'",	"C4'",	"O4'",	"C3'",	"O3'",	"C2'",	"O2'",	"C1'",	"N9",	"C8",	"N7",	"C5",	"C6",	"N6",	"N1",	"C2",	"N3",	"C4"]

AA_in_radius = []
data = pd.DataFrame()
radius_df = pd.DataFrame({"Running_Total": [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]})

#~~~~~~~function to screen out cases where multiple ligand atoms see the same residue
def to_uniqueness(x):
  return list(dict.fromkeys(x))

for radius in aradius:
        print('Opening Loop Radius:',radius)
        #ic(radius)
        now = datetime.now()
        dateTimeObj = datetime.now()
        timestampStr = dateTimeObj.strftime("%H%M%S%d%b")

        p = PDBParser(PERMISSIVE=1)

        for filename in pdb_files:
                
                #~~~~~~~use the filename to get the 4 letter pdb code
                structure_id = filename.rsplit('/Desktop/pdb_files/'+ligand_id+'\\', 1)[1][:-4]
                structure = p.get_structure(structure_id, filename)
                models = structure.get_list()
                ligandc = 0
                for aModel in models:

                        #~~~~~~~get all the chains in this model
                        chains = aModel.get_list()

                        #~~~~~~~get all the atoms in this model
                        model_atoms = Selection.unfold_entities(aModel,'A')

                        #~~~~~~~create a NeighborSearch
                        ns = NeighborSearch(model_atoms)

                        #~~~~~~~search the chains for the ligand_id
                        for aChain in chains:
                                residues = aChain.get_list()
                                for aResidue in residues:
                                        if aResidue.get_resname() == ligand_id:
                                                ligandc = ligandc + 1
                                                #~~~~~~~unpack all atoms in ligand when a ligand is found
                                                atom_list = Selection.unfold_entities(aResidue, 'A')

                                                #~~~~~~~pick a center, ie the coordinates of the atom

                                                #~~~~~~~Empty list to accend residues found in neighbor search
                                                residue_list_in_radius = []
                                                
                                                #~~~~~~~cycle thru each atom in ligand and get neighboring atoms
                                                for x in atom_list:
                                                        if any(str(f"{x.name}") in s for s in ligand_atoms):
                                                          center = x.get_coord()
                                                          neighbors = ns.search(center, radius)
                                                          residue_list = Selection.unfold_entities(neighbors, 'R')

                                                          #~~~~~~~Cycle through the list of residues associated with each atom found
                                                          for xResidue in residue_list:
                                                            
                                                                  #~~~~~~~Add the found residue to the list
                                                                  residue_list_in_radius.append(str(f"{xResidue.resname}_{xResidue.id[1]}"))
                                                                  
                                                #~~~~~~~Return the unique residues list from those found
                                                residue_list_in_radius_unique = to_uniqueness(residue_list_in_radius)

                                                #~~~~~~~Count up the AAs in the Unique residue list
                                                for xAA in AA:
                                                    AA_in_radius.append(sum(xAA in s for s in residue_list_in_radius_unique))
                                                        
                                                #~~~~~~~Once all atoms of this ligand have had their neighbors uniquely counted at this radius add it to the temp df      
                                                radius_df['Ligand Count']=AA_in_radius
                                                
                                                #~~~~~~~Clear the AA list ready for the next ligand
                                                AA_in_radius.clear()
                                                residue_list_in_radius_unique.clear()
                                                residue_list_in_radius.clear()
                                                
                                                #~~~~~~~Sum the atoms found in this ligand to the running total for this radius          
                                                radius_df['Running_Total']=radius_df.sum(axis=1)
                                                
                                                #~~~~~~~Clear the column associated with ligand just finished
                                                radius_df.drop('Ligand Count', inplace=True, axis=1)
                                                
        #~~~~~~~At the end of all the structures for this radius add the running total to the final dataframe
        data = pd.concat([data, radius_df], axis=1)
        data = data.rename(columns={"Running_Total": str(radius)})
        radius_df = pd.DataFrame({"Running_Total": [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]})
        
print('Final Dataframe:\n',data)
data.to_csv(r'ATP_'+timestampStr+'.csv')
print('Finished')





