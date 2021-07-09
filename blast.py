from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
from difflib import SequenceMatcher
import Bio
import itertools

headers_blosum62 = np.array(['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'])

def get_combination(blosum, k=3):
  return [p for p in itertools.product(headers_blosum62, repeat=3)]


def division(letter,k = 3):
  array = []
  for i in range(len(letter)-k):
    str1 = letter[i:i+k]
    array.append(str1)
  return array

def find_neighbors(seeds, possibilities, blosum, threshold=13):
  local_neighbors = []
  index_seed = 0
  for i in seeds:
    for j in possibilities:
      score = 0
      for k in range(len(j)):
        score += blosum[i[k],j[k]] 
      if score>=threshold:
        local_neighbors.append((index_seed,''.join(j)))
    index_seed += 1
  return local_neighbors

def extend(sub_cad_db,  query, matrix, pos_cad_bd, pos_query, threshold=13):
  p_izq_q, p_der_q = pos_query[0],pos_query[1]
  p_izq_w, p_der_w = pos_cad_bd[0],pos_cad_bd[1]
  value = 0
  p_temp, q_temp = sub_cad_db[p_izq_w:p_der_w+1], query[p_izq_q:p_der_q+1]
  for i in range(len(p_temp)):
    value += matrix[p_temp[i], q_temp[i]]
  if value<threshold:
    return (None, None, -1)
  while p_izq_q >= 1 and p_izq_w >= 1 and p_der_w+1 < len(sub_cad_db) and p_der_q+1 < len(query):
    p_izq_q=p_izq_q-1
    p_izq_w=p_izq_w-1
    p_der_w=p_der_w+1
    p_der_q=p_der_q+1
    value = 0
    p_temp, q_temp = sub_cad_db[p_izq_w:p_der_w+1], query[p_izq_q:p_der_q+1]
    for i in range(len(p_temp)):
      value += matrix[p_temp[i], q_temp[i]]
    if value < threshold:
      p_izq_q=p_izq_q+1
      p_izq_w=p_izq_w+1
      p_der_w=p_der_w-1
      p_der_q=p_der_q-1
      break 

  value = 0
  p_temp, q_temp = sub_cad_db[p_izq_w:p_der_w+1], query[p_izq_q:p_der_q+1]
  for i in range(len(p_temp)):
    value += matrix[p_temp[i], q_temp[i]]

  return (sub_cad_db, sub_cad_db[p_izq_w:p_der_w+1], value)

def search_database(neighbors, DB, matrix, query, threshold = 22):
  salida = []
  for (index_seed, neighbor) in neighbors:
    for cadena in DB:
      #print(cadena)
      #print(neighbor)
      position_find = 0
      position_initial = 0
      while True:
        position_find = cadena.find(neighbor, position_initial)
        if position_find == -1:
          break
        valores = extend(cadena, query, matrix, (position_find, position_find+len(neighbor)-1), (index_seed, index_seed+len(neighbor)-1), threshold)
        if valores[2]!=-1:          
          salida.append(valores) 
        position_initial = position_find + 1
  salida.sort(key=lambda tup : tup[2],reverse=True)
  return salida

def main(LETRA = "GAGTT"):
  Datset_fasta = []
  sequences = SeqIO.parse("seqs2.fasta", "fasta")
  for record in sequences:
    Datset_fasta.append(record.seq)
  
  BLOSUM62 = Bio.Align.substitution_matrices.load("BLOSUM62")
  headers_blosum62 = np.array(['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'])
  print(headers_blosum62)
  a = division(LETRA)
  all_combinations = get_combination(headers_blosum62)
  neighbors = find_neighbors(a,all_combinations,BLOSUM62)

  return search_database(neighbors,Datset_fasta,BLOSUM62,LETRA,11)