import sys, re

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
    return a/len(ADN), a/len(ADN), a/len(ADN), a/len(ADN)
