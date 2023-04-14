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