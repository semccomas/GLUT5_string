from modeller import *
from modeller.automodel import *


class MyModel(automodel):
	def select_atoms(self):
		return selection(self.residue_range('27', '31'))
	def special_restraints(self, aln):
		rsr = self.restraints
		rsr.add(secondary_structure.alpha(self.residue_range('450', '461'))) #F to N is a helix


env = environ()
env.io.atom_files_directory = ['.']
env.io.hetatm = True

a = automodel(env, alnfile='rat_cow_in.ali', knowns=('cow_in', 'rat_out_N_C_terms_modeling'), sequence='4YBQ:rat_out_fasta', assess_methods=(assess.DOPE))
#a = automodel(env, alnfile='rat_cow_in.ali', knowns='cow_in', sequence='4YBQ:rat_out_fasta', assess_methods=(assess.DOPE))

a.starting_model = 1
a.ending_model = 200
a.make()
