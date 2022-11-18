#!/usr/bin/python3
import os
import sys
import subprocess
import pandas as pd

#step1 
#protein1 = input("what is your query protein: ")
#organism = input("what is your organism: ")
#ask user to input the query protein and organism
#subprocess.call(f"esearch -db protein -query '{protein1}[protein] AND {organism}[organism]' | efetch -format fasta > protein.txt" , shell=True)
#Use esearch then efetch to get the fasta seqs

#pro_count = open("protein.txt").read().count('>')
#count how many sequences are in protein.txt
#f = open("protein.txt")
#org = []
#for line in f:
#  if line.startswith('>'):
#    m = line.split("[")[1]
#    n = m.split("]")[0]
#    org.append(n)
#    set(org)
#    print(org)
#first create a new 

#org_num = len(set(org))
#if pro_count > 1000:
#  print('The number of sequences is too large to proceed to the next step, please exit or re-enter the protein and organism')
#else:
#  if org_num > 1:
#      print("The number of sequences is " + str(pro_count)+'\n' + "All the sequences are from " + str(org_num) + " different species")
#  else:
#      print("The number of sequences is " + str(pro_count)+'\n' + "All the sequences are from the same specie")
#
#step2
f = open("protein.txt")
acc = []
for line in f:
  if line.startswith('>'):
    a = line.split(" ")[0]
    acc.append(a)
    print(acc)
#Extract all accession numbers as a list
f = open("protein.txt")
spe = []
for line in f:
  if line.startswith('>'):
    b = line.split("[")[1]
    c = b.split("]")[0]
    spe.append(c)
    print(spe)
#Extract all species name as a list
f = open("protein.txt")
seqtotal = []
seq = []
for line in f:
    if line.startswith('>'):
        seqtm=''.join(seq)
        seqtotal.append(seqtm)
        seq=[]
        d = line.split("]")[1]
    else:
        e = line.replace('\n','')
        seq.append(e)

seqtm=''.join(seq)
seqtotal.append(seqtm)
seqtotal.pop(0)
print(seqtotal)

length = []
for i in seqtotal:
  g = len(i)
  length.append(g)

print(length)

s1 = pd.Series(acc)
s2 = pd.Series(spe)
s3 = pd.Series(seqtotal)
s4 = pd.Series(length)
pd.DataFrame( { 'Accession number' : s1, 'species' : s2, 'sequence' : s3, 'length' : s4} )


