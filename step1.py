#!/usr/bin/python3
import os, sys
import subprocess

#step1 
protein1 = input("what is your query protein: ")
organism = input("what is your organism: ")
#ask user to input the query protein and organism
subprocess.call(f"esearch -db protein -query '{protein1}[protein] AND {organism}[organism]' | efetch -format fasta > protein.txt" , shell=True)
#Use esearch then efetch to get the fasta seqs

pro_count = open("protein.txt").read().count('>')
#count how many sequences are in protein.txt
f = open("protein.txt")
org = []
for line in f:
  if line.startswith('>'):
    m = line.split("[")[1]
    n = m.split("]")[0]
    org.append(n)
    set(org)
    print(org)

org_num = len(set(org))
if pro_count > 1000:
  print('The number of sequences is too large to proceed to the next step')
else:
  if org_num > 1:
      print("The number of sequences is " + str(pro_count)+'\n' + "All the sequences are from " + str(org_num) + " different species")
  else:
      print("The number of sequences is " + str(pro_count)+'\n' + "All the sequences are from the same specie")
  




