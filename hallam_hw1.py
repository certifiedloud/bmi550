################################################################################
# ===================
# Hallam HW 1 BMI 550
# ===================
#
# PseudoCode
# ----------
# read fasta file
# read coding strand
# convert coding strand to template strand
# convert template strand to RNA
# convert RNA to protien starting at 0 index, 1 index, and 2 index
# For each of the three protien sequences print the index where the cleavage is
# and the Amino Acids that are being cleaved.
################################################################################
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
import re

# Strings of basic amino acids to cut
basic_amino_acids = "KK|KR|RK|RR|HH|HR|RH|HK|KH"

# regex method(print index and basic_amino_acids given protein)
def cleave_sites(protein):
    p = re.compile(basic_amino_acids)
    for m in p.finditer(str(protein)):
        print m.start(), m.group()

for seq_record in SeqIO.parse("pa1.fasta.txt", "fasta"):
    coding_strand = seq_record.seq
    template_dna = coding_strand.reverse_complement()
    messenger_rna = coding_strand.transcribe()
    protein0=messenger_rna[0:len(messenger_rna)].translate() 
    protein1=messenger_rna[1:len(messenger_rna)].translate() 
    protein2=messenger_rna[2:len(messenger_rna)].translate() 
    print(seq_record.id)
    print"The coding strand is: \n", repr(coding_strand)
    print"The template strand is: \n", repr(template_dna)
    print"The messenger RNA strand is: \n", repr(messenger_rna)
    print
    print"3 Proteins"
    print
    print"The Protein strand indexed at 0 is: \n", repr(protein0)
    cleave_sites(protein0)   
    print
    print"The Protein strand indexed at 1 is: \n", repr(protein1)
    cleave_sites(protein1)   
    print
    print"The Protein strand indexed at 2 is: \n", repr(protein2)
    cleave_sites(protein2) 
    print
