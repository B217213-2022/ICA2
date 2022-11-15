#!/usr/bin/python3
import os, sys
import subprocess
 
protein1 = input("what is your query protein: ")
organism = input("what is your organism: ")
#ask user to input the query protein and organism
subprocess.call(f"esearch -db protein -query '{protein1}[protein] AND {organism}[organism]' | efetch -format fasta > protein.txt" , shell=True)
#Use esearch then efetch to get the fasta seqs

pro_count = open("protein.txt").read().count('>')
#count how many sequences are in protein.txt
if pro_count > 1000:
  print('The number of sequences is too large to proceed to the next step')
else:
  print("The number of sequences is " + str(pro_count))


