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
python3 main.py fasta_files/mullus_surmuletus_cytb_cds.fasta midi_files/mullus_surmuletus_cytb_cds.midi
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
from Bio.Data import CodonTable
from Bio.Alphabet import generic_dna

## MIDI 
#import mido
#import rtmidi

from midiutil import MIDIFile


import csv





#================================================================================
#ARGUMENTS
#================================================================================

parser = argparse.ArgumentParser()
#creer un argument positionnel (tu mets le nom du fichier apres le script en gros)
parser.add_argument("fichier_fasta", help="input a FASTA file")
parser.add_argument("fichier_midi", help="output a FASTA file")


#================================================================================
#MAIN
#================================================================================
#recuperer le nom du fichier en argument
args = parser.parse_args()
ffasta=args.fichier_fasta
fmidi=args.fichier_midi


#================================================================================
## record sequence
for seq_record in SeqIO.parse(ffasta, "fasta",alphabet=IUPAC.unambiguous_dna):
    seq_record_information=seq_record.id
    seq_record_sequence=seq_record.seq
    #imprimer la sequence au format FASTA
    #cette ligne de barbare corrige quelques defauts (enlever les guillemets,convertir au format string...)
    maseq=str(repr(str(seq_record_sequence.upper()))).replace("'","")
coding_dna = Seq(maseq, generic_dna)
# translate dna
#trans_dna=coding_dna.translate(table=2)



#================================================================================
## load table to convert note name into degree MIDI value
noteToDegree={}
with open("midiDegreeToNote.csv", 'r') as csvFile:
    reader = csv.DictReader(csvFile)
    for row in reader:
        splitNote=row['note'].split("/")    	
        for sn in splitNote:
            noteToDegree[sn]=int(row['degree'])
csvFile.close()

#================================================================================
## load table to convert amino acid to coding note
dnaToNoteDic={}
with open("AminoAcidToMusicalNote.csv", 'r') as csvFile:
    reader = csv.DictReader(csvFile)
    for row in reader:
        dnaToNoteDic[row['aminoAcid']]=row['codingNote']
csvFile.close()

#================================================================================
## convert DNA sequence into musical note sequence
noteSequence=[]

for dna in coding_dna:
    noteDNA=dnaToNoteDic[dna]
    noteSequence.append(noteDNA)

#================================================================================
## convert note sequence into midi degree sequence

degreeSequence=[]
for note in noteSequence:
    degreeNote=noteToDegree[note]
    degreeSequence.append(degreeNote)


#================================================================================
## define track MIDI file

degrees  = degreeSequence
track    = 0
channel  = 0
time     = 0    # In beats
duration = 1    # In beats
tempo    = 60   # In BPM
volume   = 100  # 0-127, as per the MIDI standard


#================================================================================
## write midi file
MyMIDI = MIDIFile(1)  # One track, defaults to format 1 (tempo track is created
                      # automatically)
MyMIDI.addTempo(track, time, tempo)

for i, pitch in enumerate(degrees):
    MyMIDI.addNote(track, channel, pitch, time + i, duration, volume)

with open(fmidi, "wb") as output_file:
    MyMIDI.writeFile(output_file)





#================================================================================
#END OF PROGRAM
#================================================================================