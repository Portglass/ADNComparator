import pandas as pd

def read_file(path):
    f=open(path,'r')
    header = f.readline()
    ADN =""
    for line in f:
        ADN += line[0:-1]
    ADN = list(ADN)
    return header, ADN

def nbNucleo(ADN):
    return ADN.count("A"), ADN.count("T"), ADN.count("C"), ADN.count("G")

def pourcentNucleo(ADN):
    a,t,c,g = nbNucleo(ADN)
    list_pourcent = [a/len(ADN)*100, t/len(ADN)*100, c/len(ADN)*100, g/len(ADN)*100]
    source = pd.DataFrame({
        'Pourcent': list_pourcent,
        'Nucleotide': ['A', 'C', 'T','G']
     })
    return list_pourcent,source

def fromFileToDataframe(path):
    header, ADN = read_file(path)
    return pd.DataFrame(ADN,columns=['Nucleotide'])

def fromStringToDataframe(seq):
    return pd.DataFrame(seq,columns=['Nucleotide'])

def complementaire(seq):
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

def transposer(seq):
    return reversed(seq)

def nbMismatch(seq1,seq2):
    compt = 0
    if len(seq1) != len(seq2):
        print("Pas la même longueur")
        return -1
    else:
        for i in range(0,len(seq1)):
            if (seq1[i]=='A' and seq2[i]!='T') or (seq1[i]=='T' and seq2[i]!='A') or (seq1[i]=='C' and seq2[i]!='G') or (seq1[i]=='G' and seq2[i]!='C'):
                compt += 1
        return compt

