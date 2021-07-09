from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO

muscle_exe = r"/home/hermogene/Documents/bioinformatica/muscle3.8.31_i86linux64"
in_file = r"/home/hermogene/Documents/bioinformatica/opuntia.fasta"
out_file = r"/home/hermogene/Documents/bioinformatica/aligned.fasta"
cline = MuscleCommandline(muscle_exe, input=in_file, out=out_file)
stdout, stderr = cline(in_file)
align = AlignIO.read(out_file, "fasta")
print(align)