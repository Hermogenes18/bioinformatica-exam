#!/usr/bin/env python

def estimate_nucleotide_frequencies(seq):
    seq = seq.replace('-','').upper()
    A = seq.count('A')
    C = seq.count('C')
    G = seq.count('G')
    T = seq.count('T')
    length = float(len(seq))
    return [ x/length for x in [A,C,G,T] ]

def pdistance(seq1, seq2):
    p = 0
    pairs = []
    for x in zip(seq1,seq2):
        if '-' not in x: pairs.append(x)
    #for (x,y) in zip(seq1,seq2):
    for (x,y) in pairs:
        if x != y:
            p += 1
    #length = (len(seq1) + len(seq2)) / 2
    length = len(pairs)
    return float(p) / length

def JCdistance(seq1, seq2):
    from math import log
    b = 0.75
    p = pdistance(seq1,seq2)
    try: d = -b * log( 1-p/b )
    except ValueError: 
        print ("Intento tomar el registro de un número negativo")
        return None
    return d

def K2Pdistance(seq1,seq2):
    from math import log, sqrt
    pairs = []

    #collect ungapped pairs
    for x in zip(seq1,seq2):
        if '-' not in x: pairs.append(x)
        
    ts_count=0
    tv_count=0
    length = len(pairs)
    
    transitions = [ "AG", "GA", "CT", "TC"]
    transversions = [ "AC", "CA", "AT", "TA",
                      "GC", "CG", "GT", "TG" ]

    for (x,y) in pairs:
        if x+y in transitions: ts_count += 1 
        elif x+y in transversions: tv_count += 1
    
    p = float(ts_count) / length
    q = float(tv_count) / length
    try: d = -0.5 * log( (1 - 2*p - q) * sqrt( 1 - 2*q ) )
    except ValueError: 
        print ("Intento tomar el registro de un número negativo")
        return None
    return d    

# Quick test:
a = 'AATGAGCCGGTGGCTGATCCTTCCTCTGTACATGGGCCTCATAGCCTAGTCTGATGGTGATGAGCTGCCCAATTATAGCAGGGTGAAAATGCGCAACCCGCCGAAAACTTGCCACGGATTGGATCTTACATCTTGCGGGGGTACATCGACGATCAGGCACTAGTTGATCCTTCACACGAGTTCAACAGCTATTGACTCCCTCCAGAGGAACGGAGACTTACACTCTAGTTATTCTTTGACTGCGCCTGGTCGGCCCCCCCTGTACAGGTTGTACGAAAAAATAAGATCACCGCGCCAAGGCCCCCCTCCCGACCGCACACAATGCTACGGCCGTAATCTGGTTCTGTGCATATAGAGCAACTAGCATCACCCCGAAATGCAGCCGGTAGTTGGAAATATGTCGACGGCCTAAGGTCGGTAGCATTGCCCTTAGAGCGGAATACCTGGGCTACCTAAATCTTTCTCGCTCACTACCACTGGCCTAAGCCGATTAGAGGGCCAATTAGTTATTGTGCGTTAGTCGGACGGCAGTACGTGGAGACTTGGCGTGGAGATGCTAATCAGCTGCGCGCTGTCCAGTGGATCTATGCAGGGCATGGTGAGATCAAAATATCGGACGGAGTTGGTTACCCAATATCTCACAGAAGACCGGCTACCAAAAGGTAAGCAATCGCTCGGGAGGACCGAACCGGGGACTAGGTTACAATTTGTTCGGTACAACGCGCTATCCATGACAGGGTATTCCTTTGTAGTTACGCCCCATGCTTATGCCCCCGTTTTACCGTATTCGCATGAGAACGCTGAACACGGAGACGTACGTCTGAATCGCAGTGCTCTATCTGCCCCGCCGAGTTGGATACCTGGACCGCGACAAAGAATCGCTTAGAATTGTAGGTATTCGGGGATCGCCAAGGGATCGGTGCGTGCGCGCCTTAGCATGTACTCATCCCTAGCGCCCAGAGAAGACAACGGCCCCACTCGCAAACTTACGAAGACCTCCCAATCTCTCAAAACCTCACTAGTCCACTCCCAGTCACTGGAGGACAGCATTCGGTGAATAGCACTACTTACACAACGCCCGGACGCGTATCGGAGTATGCTGAGAATAATTTGTTATCGAAGCAGTAGGGTTGGGCTCTTCCGAACTGTTGACCTCGCGCGCTCGATTCTGCCCATGCAGCGGAAGGCACTCCGATGCGCATAAGTTCTGCAAAGCCACCCTGACTTCCAAGATCCGGGTGCTTACTACCAGTTGTAGCCTCTGCCAGAAGAGAAGTGGTGCTCAATGGTCAACCGTGTCCTCAGGCACCTTCAATTCCCCCCGTCCAAAGCGACGGGTATCATACCTACCCGCTAGGTTTCATGAGTGCAGAGTTTCATTTCCTGAAGCGAGGAGGCCAAATCTACACTTGCTTATATCTAGGTTAAATCTCCTGCTTGACATGTATCCTATTCGAAGTTGTTCGGTCACCAACGGTTAGTTACTGTTGCTCCGGTCTAGGCCCATTGAGGGACCACAAAGCACCGCGCAACTAGCACCCCGTTGTGCCGAAGTTACGAAGTTATGCCCCCTAATCATTAGTAAAAATCTTAACGGTCACCTTGGCTGGTATCTACTGCCTACGCCGAGCCTATCACGTGGCCTGGGCACAATATTCCGAGACACTGCGTATCTGACCTAGCTATCCATTTGACTATTCAGCGGTGTCCACTCTGCTTTATAGTGTTGATCTGGTAGGACAGCCATCGCGTTATCGCCCTGGCCAGCAAACCTCTACAGAGAGGTCACTGCAGAACCGTATAATTCTGGGGTACTACTGGGTAGGCCTTGGCTGGTCTTGTCCAACCCCAGGCCGTAGCGGATATGAGAGAACTTAGCTCAAACGTAATCGCCTGTTGCACGCACTTCTCGAGAGGATGACACTAATTGGCCTACGAAAATAGGACCAACTAATTGTCAAAAACCACCACTTCTTCACTATGAGGTGAGGAACGCGCCTTTACTGAGAACGTATCGTAAAATGAATAATCGACTCTCGGTCTCACGGCCGGAGATCTGCTGTGGGAGGTTATTTGTGAGGACGACAGTTGGCTCCAGCGCTGTAGAACGCCGGGTTTCTACTTTCAGAGACCACGGGTAAACAAATTCGTAGGATTTCGACTCTCTAAAATCATGACGACCGCGCTCTATAGTCGACTGTGTATTGTAGAGAAAAAGCCAAGTAAGTTGCGGGTTTCCCCCGGGCACGTGTACCGGGTTTACGTTTAAACGGATTGACGGGACTATAGATAACATTGTTTAGACCACGCATCTTAGTGAGTGAGTCGATTTTCTATTGTTGGATAGATGCCCGGACTAGAACGGTGATTAACCCCCGGCCATTCTGGCGTATGTGATGGCGCGGGGGAATTTTCACCAGATTCGCTTTCTGGTCCCCGCCCTTCGGGTTCCGACTTGTAGAGAGTCAACAGGTACCCCTAGGCGTAATTCTAAGTGGGGAACTAAGAGGTTTTCGTATTCCTCTCCATTACGTCATCGGATAAGCACAATGTAATGTGTAATAGCTCCCTTTAGTTGGAGTGTAGCCATATCAAAACTGAGCACGATAACTGACCTTCCTCTTTAATGATCGGGGCAACTACAATAACATGGGACCGAATTAAAGATAGGCTGAGCTAGCCTGATGTACCTACGCGCATTTGGAAATAGATCGCTCAACTTCCCGGCCCTAAGCAACTACCAGACGTTTGTGGTACGATGATCGGGAGTCCCCGGCGGAGGCATCGTCACGAAAACTGTCTGTTCATGCGGACAGTAAATCCCTATCTATGATGTTCGCGGGCTCATCTACCAGAAACGTCCCCCGGACCCGAACAGTTTCGGTGTAAATACTTGCTAGCCGGGTCCGAGGCTTCCGTGCCGCACAGTACTTTAGTCGGAGGCAGGTGACTGTTAGATTGTGCGAGCGCTCAAAGCGTGATATAGTTACGGGGGGTCTAATGTCAGCTGTGACAAAGTCCTGTCTGTCCGTTTTTGCCCAGACACGCTCGTCTCTCCATGCTCTTTTACAGTGCTGGAAATGGGAGCTTCTGGTGCATTGAAACGAACCTTTTAGCAGATGTGGCATACCTCTTGGTACAGACGAAGTATAGACTACACAAAAGCGCTCGAGGCAAGTTTTGCGTAATCGGCCCGTTTCGCGCAGGTTGGAAACCTATACTGCACACTCCAATTAGTTAATTCACACTAGATTTTGGATCTTACGTCTCAGCGTTCCACTTGCGAATGAACATCCAGTCGATTAACTACCAGCTATACTTCCCTTAGATTAATTGTGACCAAATTAACGCATGTTGATCCGAACATCTGCTGTTGACTCGCTGTAGTCTTACGTGCTAAAATTTCGTCACCCCAACAGCATCCTTATGGTTTAAAGGTAGCCCCATCAACGGAATATTGCAAACTACTAACGTGGTTGGCTCGAAGATTTTTACCTTTGGTTATGCCTCGGACCTGCCGTGCTGTATATAGGTCCCAGAGTCACACATGGATGTACTACTTTCATAATACCATGGTGAACCACGATTATACTAGGTATGCGTGGAGTGATCAGTTTGTGCTACGACCGAGCGGAGCGGGCGTGACATGGGCGGTGTACGACGACCTTCTGGAAAGGACCCCCTGGTTCAGTTCAGCAACAGGCGGACTTGTCTATAAACACGTCTTAAGGAGTAAAGTAGAGGCGGTTTGCTAGGTGCAGATTGTTTGGCCACGCTGCCTGAGCGTAATGTTAACGAATTATTGGACAGTTAGGGCGTTCGAGTGGGATGTTAATTTGGCCTCTTGCACCGGAACCGTTAGCCGGACCACCGGCTTCTGTGGGGGTCCGTAGGAACACAGGTGCCTAGAGCGCAGCAGGCTGTGCCTGTGGGCCCGCCAATGAAGCCAAATCAGATTTCGATGAGCATGATCAGCACACGCCAGTGCGATACCGTTCTCTCGGTAGGGCGTCCCAATGGGAAGAGCGCAGTGAATCTTGTGGGTGCGGACTATAACGGATCCGTGATTGCCTCCGGCGGAAATGAAGTCTGTGGTAGTGTTTTAATAAATATAATGCGCCAACGTGGACTCTCGGATAACGTCGGTGTATCCGGCATATCCAGACACCGTTGAACAATTTAAATCGTCTGCCCGTGCTCTAGAGTCCTCATAGTCTAAGTGTAGTCTATTCCCTCTAAAATTTTTGAGATTAGTAAGTCTTGACTCTCAGGTCTTGTTCTTAAGCGGGGCAATGCTCGACTAAATAACACGATCACGCCCACGCCCGTAGCCAGTGCCTCAATCTATGCGCACGTGAGGTAACTGCGGACGGCGTAACGCGTTCTGGATATTATATCCTTGCCGTGGACCCGGTCAGAAGCAGGGTCGGAAGACCCCAGTACTGCCGGCCCGGGGCGTACAAAGTGTCTGTGCTTTGCGCCGGGGGACGGTAAGGCTTTGTGATCCACCGTCAGTTCCGCTAATAGTGTTCACTGGTAATACCTAGGGTAGCAACGCACAGATTCGAGAATTAAGCACATAGGTTTGACCGACATGGTGGGCGCACACGACCGAAGTAGATGTACGATGTTGACGTTTACGAACACTATCTAATTCTCTGCAATATGTAGGAATGCAGCCCAAAATTTACATTGAGTAGCCGGACGAAGCTCTGTGATCTAGGGCATAATATTGATAGCACTCCAAGACTGTCAACTTATAGATGAGCGAGAACACTGGTCCTGAGCCCCAAGGGAAAGCTCTGCAAGAGGCACCTCCTTGACCAAGGAGCCCATTGGTAGCCGGTTACATAGAGTCCGCAAATCTCGAAACCGGTGTCCAGATGAGAGATAAAGAACTTGTTGTTGGCTGAGTGGGCAGAAAACACCGATGCTGCGCTTGCGGACGGCATGCAAGCCACTCAGGCTCCCGTTGAAACACGCGCCCCTTCTGTCAAGCC'
b = 'ATTCCCCGGTAAAAGCAAAATTCTGTAGTCAGACGCCCGTACCACGTAGATGAGGATCTGCCCCCGTCATGAGGGCGACCACAATAGTAGGCGGCCACTGTCGGCAGTTACAAGAGCAGCGGCGGATTTCATCGGTACAATCGCGGTTGCGCCGCCCGCCTAGTAAGCATTGCGCTCGAGTCCCTGTAATTCTCAGGCTAGGCAAGCACACAGATATGTGAGTGTATCGTCATAATGTGCTCGTTTAAGAGCTCATAATTTGTTAACTGTATGGCCGCATTTTCAAGATGGATCGTGCGCCCAATAACAATCCGGCAGTAAGAGTTCGCTCCAGGAGATTGCCTGGCGCGCTAGTACTACGGTAGAGACGTCAAACACACGAGCGGGGCGAGGATTTTGGCTTGGGTTCTAGCTACCCTAAACCAGCGGGTATCGGCCCAAGGGCTTTTCTCTCCTGTTTCCAGCGAAATTAGGGGTACGGTGGCAATCACTGAAATCACCCCCGATGTGTAGTGACGCTATACCTAAGTCTAAGAGAATCTGGGGCAAGTGTCTAAGGACTGCTGGCTCCAAAGAGCCTGGATACACCGCTGTACCTTTGAGACGTTGGTGGTATCTAATGCATAACGCGCTCGGAGCGAACGCCTCCACAACGGTGCGTTTCCCGTATCGTATAGAGGATAAGGATATCGTCTCTATTAAAGCCCCAATGGGGTGCAACAGCAACGAATAAGTCTGTATTTAGAAGAGCGAAAAAATTCGCAACCTTTAGAAATCCTTAGGCCTACGGCTATCCACGCATTGGTCGGGGCACTAAAGTCATAGATTCACCACGATCCCCCAGTTTCGGATAGGAATTGCGTCGTGGGCGTGATGAGTCTACGTCAACCAAGCCAGCTACATGCTCAGTGCAGCGATAATTTGCCTGGTGAGGATAATACAAATCGCATTGCGGAATGCTTTCTAACCTACGTTAGGTCGGTAAGTTGCCATGCTTGTCGGTGTTGTACTACTCTCAGACTTAAGATGTTCAAATGGACGGTATCCCATGGGTAGCCAGTGATGTGGATGTCGTCGGGCCTATAAAGGGCGTGGCGTGGCCCCTAACCAATCCACTAGTCCGGCGTAAGCAGCACCCTCACGGGGGAAAGACCGAGCATCTCGAGCAGTTCGTTCCACGATATTATCAGGTGACTTCCTATCTAGGTAACAGTCTCCAGCTGCAGACGCCTATTGGTACCCCAACGGCGCGGAGTTTAATAAGTGACCTAGCATCCTGGGTCTCCCGGTATACAAAGTCGAACGTCGACGTAGTCCCCCCCGTGACGCTATGAGGCTGGTCGGGAGGTTTCCGTAGACCGAGAGTGATACACACCAAAGGTTCACAATGTTGTTTGATATGCAAACGTCACCTTAATGGATAAGAACAGTGCTCATCAGAAGTCAGAATACAGATATCAATAAGAACGGCCTCACGGCCGTTAGAGCTATTGAACGAAAGAGCGCCTCCAACACATAGGATCGGCCAATCGACAGTCTAACCTGGCACATTGTGGGCCGCGGCAAAGCGTCCGCGACAGCGTTGTAGGTCCGGAGCCTTATCTAGAAGTGAGATCAGCACAACATCTCATAGTATCCCATGGTTTCCCCTTCGGAAGTATATGGAGAGGTGACATGTAACCGAGGGTCATCAACGGACTAGGTGATGATTCTTTCCCGTGAAACCCCATCTGATGCTGTGTCTGCGGGACTTGGATTGTACGGAACGACTCGGCTTGTTGAATCATCGAATTTGGGGAGGGGTGCGTATGTGGTGGGTTTGACAACTAGCTCCTTCCCTCGATCTTGAAGTCAGTTGTCTGAACGTTTTAGACATATGACATCCCGACCGCAGGGAAGAATCGTCTGATAGGGGATTGGAGCTTGATAGTAAGTGGCCCCTAATACGGCTGTATAGGTCGTGAAGGTCAAGGTAACGGTGACATATAAACACTCCCCTGATTGTAGTATAGCACATCGAGTTGTACACGAGTACACCTCTTATTGACCAATGGGCATCGCTTCCGCTACGGCCCCAAAGGGTGGCCGCCTGAGCCCTGGGCTAACGACGATGAAATAATCAGCTTGTCTCGTGTGCAGGCATAGATTAGGTGACCCAAGGAGCAATAGCCATTGGCTCATCTTAGGGGTTCAACGGACAAATTTTGTCGAGAAGACACAATCCCGGATCAGTTGAGTGGCAATGCGAGATCTTCCTCCTCTTAATTACGTCCCTTATCGTTTGGCTTCTGGAAGACGCGATAACCACGGAATGCAGTGTAGCATTCTTTGCGCGCAAAGGTTCAGATGGAACTACAAAGCCGCATGGTTATTCCCCGTCGAAGAAGTCGCACCTTACTTTTGCACTCTTCGTGCTCGAGACCAAAAGTGCGAATTTAGGGGCTCCTAATGGGTCGCCTAGACTTGTATCGACTACACACATGCCTTCATAACTGTGTGATTTATAGAGCTAACACACTTGCGAGTTGCACGTTGAGCTCCCGCGTAAGTGTTCTCTACAGCCCGATGCCTGCTTGTCTTCTGTCCGCGGCTTGTTGCAGGTCTGACCGCCATGTTCGTTCCCTACTAACACACAGACGGTGTATAGCCCTATCATATCGATCACTTAGGACTGTGTTTTATCGGAACCGAGGTGACCATGCGTAGTCGACTAGTTGAGCGAGAATGCCCAGACGTTAAAAGGACAACGTCGCGAATGGCGTTGCGCGGGTATCCAATCGCGTAGCCAGAGGTAACCAAGTGACCCGATTAACTCCGACATAACACTAGATTGAGCGGTCAGGACTTTCGCCAAGAAGCGTCAATGGGTGGGCGATGGCGTATACTCGCTCTAACCTCGGAGAGAAGCCATTCAGCATTACCAGAGAACTGTGGCGATAGGTAACCGCTTTCTCAACTATTATACATTCGGTGCACCCTGCTTTAGGCATGTTCGAATCGCTGGGAACACAAGCCATCCAACTTCGCAGTCTGGGTAGCCGCGGACGTTTTAACGGTCGTACAGTTGAGCGCCGCGGTATTGGTACGTTCAGGCTAATTCGATCTCCGTTGGCGCCCGGTGCCGTGAGACTTGAATGAATGCCCTATGGAATGTGGTAGACATTGGTAGCTGAGTATTGCGAAGTCCGCGAATAGGGCCATACCAGGTTAGTAAATGTAACCAAAGTAAGGCTTGGTGTTGCTGACCGAGTTTTTCGTTTTAGTAGATATCGTAACACGACTGACGATGCCCTGCTTACCCAATTTAGTTGACGTTAGATAAGATTGATTCGATAGAAACCTTTCCTGTTTCTTAATTTCTTATCATACTATATTGCGGACGGGTGTAACATAGCGGTTTAAAGGTCGACAGCTGAAGCCTTGTTAGGCTGATTATGCCTCCTCACCACTCAATCTAACGAATTTTCCAAAGCGGTCGCTTTAAACTAACGGCGAGACTATGCGTGCTCGGTCATTTCCGCCTACGCGAAATCACGTCGGAGCGCGGACCCGGCGGGTAGTGACGGCGAACTCTACTCCCCCGAAATCAACCAAGGCACGTTGCGACCACTTTTGGGGTGGAGCCGCCCACGCAATCCCTCAACCGATATTTGGAGTAGCGGACGCCCCGGGTATCCTAGGTCTGGAGCACCTATCCGTACAGCGATGGCTGCCTACCAGTACTGCAAAGAGACCTACATCGTTTCTAAGAACCTCTTACACGTCCCGATGCACAGCATGACTCCCCTGGCTGATTCCCACCCTTGTGGGACCGCTTAGATGTTAACGGTTCTCAATGTGCATACTAGACCCAGGGGCCGGTACAAGAACCCATAAGGAGAACACATCTATCAACGGTCGTGTATAACAAGGCAGCGTTGTTTTGGATATCAGATGTAGATACCATCGTACTGTACCACCCTAGAGGTGATTCCTGGAGCGAGTGTCGACGGACGGACCACGCTATGCGCCTTCCCCATCAAGATAAAGTGGGTAGCTTCGACTGACTGTGATTGAACAGATTTTTGGAGGTCCAAAGACAGCTTCGGCGCGTGCGCACGCCGAAGCACAACAGTACCTACTTGATGAGCTTGTTTCGGTTGGTAAGGTGGACCAGTACGATCCGACGTTTGTGAAATTTTACGGACAACCCCAAGAAGTCGGGAGCATGACCTTGATCTTGGGCTTGACCACGCATTGTCAGAGAACAGTCTACCTATACTAGCCACCTGTACGGTCCTTAATTATCACTATTACCATACGCTATTTTGAAACCCATCTCATTTCCTTCCGTTTCCACTCATCATGTGAGATATAAACGTACCAGGACCGCGTGCCGATGCGACGGATTCGGTGGTGTCCGTTGGCCAAACGGCCTCTCTACACAGTAGGGATACATCATCCGCCGTGTAATAGCAGGGGTGGCCGGACTAGCTTATGCACTGATAGCAAACTTTAGACCGGTCCGAATGGTCGTGGCTCGGTCTCGGTCCGCTGATATTGCTTGAGTCCCCGCTGACGTGATTCTAACTGAATGTATCGCGGTTAATGTACGCGCTTAGCTACCCCAAAGGCCATAGGCCGCGGCGACCGGTATTGCAGATTTTATAGCTCTAGTTCGTTAGTCAAGGGTCCGGATTGTTCCCATTTGCGACTGAAAAGCGCTTTAAGTGGCTACTACCCTTGTGACTAGCTAGCAGACTCTCGGCACTACGTGTATAATAACGCAACAAGTCTGTTCCACGAGAAAGCTGGAGGCAGGAATTTGTGTGAAAAGTCGTCTTTAGAACTTAGACCATGATACCGAACCTGGGGCTGTCGCTCCGTCAGCGGCGCAATTACCGTCACGGCCAGACGGTCCACCCCTATTGTAAATATAGCTTCTTGCCTGCTGATGATCGTTAGATGAGCCACCTGGCATTAAAGCCAACGGCACAACTTACCATATGAATTCGCAGGTT'
#a = 'AATGAGCCGG'
#b = 'ATTCCCCGGT'

print ("p-distance                  = ",pdistance(a,b))
print ("Jukes-Cantor distance       = ",JCdistance(a,b))
print ("Kimura 2-Parameter distance = ",K2Pdistance(a,b))
