from modeller import *
from modeller.automodel import *

env = environ()
env.io.atom_files_directory = ['.']
env.io.hetatm = True

a = automodel(env, alnfile='rGLUT5_out_occ.ali', knowns=('4zw9_glut3_out_occ'), sequence='4YBQ:rat_out_fasta', assess_methods=(assess.DOPE))

a.starting_model = 1
a.ending_model = 1
a.make()
