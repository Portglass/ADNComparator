# This page contains functions used to read a fasta file and translate the informations in it
# __________________________________________________________________________________________________
# Librairies import
import pandas as pd
from collections import defaultdict
import random
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
# __________________________________________________________________________________________________


def read_file(path):  # Open the fasta and translate the data into a string
    f=open(path,'r')
    header = f.readline()
    ADN =""
    for line in f:
        ADN += line[0:-1]
    ADN = list(ADN)
    return header, ADN


def nbNucleo(dna):  # Use the .count to return the number of each bases in a string
    return dna.count("A"), dna.count("T"), dna.count("C"), dna.count("G")


def pourcentNucleo(dna):  # Return the frequency of each bases in a string as a pd.dataframe
    a, t, c, g = nbNucleo(dna)
    list_pourcent = [int(a/len(dna)*100), int(t/len(dna)*100), int(c/len(dna)*100), int(g/len(dna)*100)]
    source = pd.DataFrame({
        'Pourcent': list_pourcent,
        'Nucleotide': ['A', 'T', 'C','G']
     })
    return list_pourcent, source


def fromStringToDataframe(seq):  # Transform a string into a dataframe for display in streamlit
    return pd.DataFrame(seq, columns=['Nucleotide'])


def complementaire(seq):  # Give the complementary sequence as a string
    seq_comp = []
    for i in seq:
        if i == 'A':
            seq_comp.append("T")
        elif i == 'T':
            seq_comp.append("A")
        elif i == 'C':
            seq_comp.append("G")
        elif i == 'G':
            seq_comp.append("C")
        else:
            return -1
    return seq_comp


def transposer(seq):  # Give the reverse sequence as a string
    return reversed(seq)


def nbMismatch(seq1, seq2):  # Give the number of non-matching bases between to sequences of the same length
    compt = 0
    if len(seq1) != len(seq2):
        print("Pas la même longueur")
        return -1
    else:
        for i in range(0,len(seq1)):
            if (seq1[i]=='A' and seq2[i]!='T') or (seq1[i]=='T' and seq2[i]!='A') or (seq1[i]=='C' and seq2[i]!='G') or (seq1[i]=='G' and seq2[i]!='C'):
                compt += 1
        return compt


def isCorrectaa(seq):  # Check if the sequence only contains valid amino acids
    for letter in seq:
        if letter != "A" and letter != "D" and letter != "C" and letter != "E" and letter != "F" and letter != "G" and letter != "H" and letter != "I" and letter != "K" and letter != "L" and letter != "*" and letter != "N" and letter != "P" and letter != "Q" and letter != "R" and letter != "S" and letter != "T" and letter != "V" and letter != "W" and letter != "Y" and letter != "*":
            return False
    return True


def isCorrectADN(seq):  # Check if the sequence only contains A,T,C or G
    if sum(nbNucleo(seq)) != len(seq):
        return False
    return True

def perc_CG(seq):  # Return the percentage of CG in a sequence
    a, t, c, g = nbNucleo(seq)
    return int((c + g) / len(seq)*100)


def summaryADNText(seq):  # Give a summary of information about the sequence
    text = "Cette séquence d'adn contient " + str(len(seq)) +" nucléotides"
    if isCorrectADN(seq):
        text += " et possède seulement des nucléotides ACTG. "
    else:
        text += " mais contient d'autres bases que ACTG. "
    text += " La séquence contient " + str(perc_CG(seq)) + "% de CG."
    return text


def seq_dic(seq):  # Return a dictionary containing every 4 nucleotides pattern of a sequence
    cods = defaultdict(list)
    i = 0
    while i < (len(seq) - 3):
        x = seq[i] + seq[i + 1] + seq[i + 2] + seq[i + 3]
        cods[str(x)].append(i+1)
        i += 1
    return cods


def freq_pattern(seq, pattern):  # return the frequency and the positions of a pattern in a sequence
    i = 0
    freq = 0
    position =[]
    while i < (len(seq) - (len(pattern)-1)):
        ii = 0
        x =''
        while ii < (len(pattern)):
            x += seq[i + ii]
            ii += 1
        if str(x) == str(pattern):
            freq += 1
            position.append(i+1)
        i += 1
    return freq, position


def random_seq(size):  # Generate a random nucleotides sequence with a given size
    Nucl= ['A', 'T', 'C', 'G']
    seq =''
    i = 0
    while i < size:
        seq += random.choice(Nucl)
        i +=1
    return seq


def translate(seq):  # Translate a nucleotides sequence into the corresponding amino acids
    codons = {
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "TGT": "C", "TGC": "C",
        "GAT": "D", "GAC": "D",
        "GAA": "E", "GAG": "E",
        "TTT": "F", "TTC": "F",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        "CAT": "H", "CAC": "H",
        "ATA": "I", "ATT": "I", "ATC": "I",
        "AAA": "K", "AAG": "K",
        "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "ATG": "/",
        "AAT": "N", "AAC": "N",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "TGG": "W",
        "TAT": "Y", "TAC": "Y",
        "TAA": "_", "TAG": "_", "TGA": "_"
    }
    cods = ""
    i = 0
    while i < (len(seq) - 2):
        x = seq[i] + seq[i + 1] + seq[i + 2]
        cods += codons.get(x)
        i += 1
    return cods


def align_seq(seq1,seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    return alignments
