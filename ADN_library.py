import sys, re

def read_file(path):
    f=open(path,'r')
    print(f.readlines())

read_file('Data/sequence.fasta')
