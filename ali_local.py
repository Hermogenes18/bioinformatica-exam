import pandas as pd
import numpy as np


df = pd.DataFrame([[2,-7,-5,-7],[-7,2,-7,-5],[-5,-7,2,-7],[-7,-5,-7,2]], index=["A", "C", "G", "T"], columns=["A", "C", "G", "T"])

df.head()

def make_match_missmatch_score(A,B,missmatch,match):
  values = set(A+B)
  an_array = np.full((len(values), len(values)), missmatch)
  np.fill_diagonal(an_array,match)
  df = pd.DataFrame(an_array, index=values, columns=values)
  return df

import itertools
def cond_diagonal(str1,str2,A,B,d,s,i,j,matrix,previus_state=True):
  if matrix[j][i] == matrix[j-1][i-1] + s[A[i]][B[j]] or previus_state==False:
    str1+=A[i]
    str2+=B[j]
    i,j = i-1,j-1
    return True,str1,str2,i,j
  else:
    return False,str1,str2,i,j
def cond_up(str1,str2,A,B,d,s,i,j,matrix,previus_state=True):
  if matrix[j][i] == matrix[j-1][i]+d or previus_state==False:
    str1+="-"
    str2+=B[j]
    j = j-1
    return True,str1,str2,i,j
  else:
    return False,str1,str2,i,j
    
def cond_back(str1,str2,A,B,d,s,i,j,matrix,previus_state=True):
  if matrix[j][i] == matrix[j][i-1]+d or previus_state==False:
    str1+=A[i]
    str2+="-"
    i = i-1
    return True,str1,str2,i,j
  else:
    return False,str1,str2,i,j


def full_order(A,B,matrix,s,d):
  matrix_values = []
  array = [cond_diagonal,cond_back,cond_up]
  for subset in itertools.permutations(array, 3):
    #print(subset)
    str1 = ""
    str2 = ""
    i,j = len(A)-1,len(B)-1
    while i>0 and j>0:
      if subset[0](str1,str2,A,B,d,s,i,j,matrix)[0]:
        bool_,str1,str2,i,j = subset[0](str1,str2,A,B,d,s,i,j,matrix)
      elif subset[1](str1,str2,A,B,d,s,i,j,matrix)[0]:
        bool_,str1,str2,i,j = subset[1](str1,str2,A,B,d,s,i,j,matrix)
      else:
        bool_,str1,str2,i,j = subset[2](str1,str2,A,B,d,s,i,j,matrix)
    while i>0:
      str1+=A[i]
      str2+="-"
      i = i-1
    while j>0:
      str1+="-"
      str2+=B[j]
      j = j-1
    cont_sum = 0
    for i in range(len(str1)):
      if str1[i]==str2[i]:
        cont_sum+=s[str1[i]][str2[i]]
      elif str1[i]=="-" or str2[i]=="-":
        cont_sum+=d
      else:
        cont_sum+=s[str1[i]][str2[i]]
    matrix_values.append([str1[::-1],str2[::-1],cont_sum])
  return matrix_values

# d es gap
def funcion(A,B,d,s=df,missmatch = -2,match = 2,m_b = False):

  if m_b==True:
    s = make_match_missmatch_score(A,B,missmatch,match)
  A = "-"+A
  A = [i for i in A]
  B = "-"+B
  B = [i for i in B]
  #print(s.head())
  matrix = np.zeros((len(B),len(A)))
  
  matrix[0,:] = [d*i for i in range(len(A))]
  matrix[:,0] = [d*i for i in range(len(B))]

  for i in range(1,len(A)):
    for j in range(1,len(B)):
        v1 = matrix[j-1][i-1] + s[A[i]][B[j]]
        v3 = matrix[j-1][i] + d
        v2 = matrix[j][i-1] + d
        matrix[j][i] = max(v1,v2,v3)
  matrix_values = full_order(A,B,matrix,s,d)
  Values_DF = pd.DataFrame(np.transpose(np.array(matrix_values)),index=["sequ1","sequ2","score"])
  Values_DF = Values_DF.T.drop_duplicates().T
  values = [Values_DF[0]["sequ1"],Values_DF[0]["sequ2"],Values_DF[0]["score"]]
  return values,A,B

# values,A,B = funcion("AGC","AAG",d=-5,m_b=False)

# pdyn = pd.DataFrame(values, index=B, columns=A)
# pdyn