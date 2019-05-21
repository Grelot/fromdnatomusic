#================================================================================
#INFORMATION
#================================================================================
#Pierre-Edouard GUERIN
#april 2019
#convert FASTA files into midi track
#================================================================================
#USAGE
#================================================================================
"""
python3 main.py fasta_files/mullus_surmuletus_cytb_cds.fasta
"""
#================================================================================
#DESCRIPTION
#================================================================================
"""
A partir d'un fichier fasta,
affiche leur sequence precedee du nom de l'espece correspondant a la sequence
"""
#================================================================================
#MODULES
#================================================================================
#outils que je prends depuis Biopython
import Bio
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#argparse permet de gerer des arguments
import argparse


#================================================================================
#ARGUMENTS
#================================================================================

parser = argparse.ArgumentParser()
#creer un argument positionnel (tu mets le nom du fichier apres le script en gros)
parser.add_argument("fichier_fasta", help="input a FASTA file")

#================================================================================
#MAIN
#================================================================================
#recuperer le nom du fichier en argument
args = parser.parse_args()
ffasta=args.fichier_fasta

#SeqIO.parse convertit un fichier fasta en une liste de sequences avec informations
for seq_record in SeqIO.parse(ffasta, "fasta",alphabet=IUPAC.unambiguous_dna):
    seq_record_information=seq_record.id
    seq_record_sequence=seq_record.seq
    #imprimer la sequence au format FASTA
    #cette ligne de barbare corrige quelques defauts (enlever les guillemets,convertir au format string...)
    maseq=str(repr(str(seq_record_sequence.upper()))).replace("'","")
#================================================================================
#END OF PROGRAM
#================================================================================
for nt in maseq:
    print(nt)