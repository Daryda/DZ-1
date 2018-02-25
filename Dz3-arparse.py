# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 16:40:53 2017

@author: Darida
"""
from Bio.SubsMat import MatrixInfo
from Bio import SeqIO
import numpy as np

import matplotlib.pyplot as plt

import argparse

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Needleman–Wunsch algorithm')
    parser.add_argument('-m', '--mistake', help='Penalty for gap or insert', metavar="float", type=float, default=-5)
    parser.add_argument('-W', '--W_output', help='Print W', action='store_true')
    parser.add_argument('-l', '--label', help='The Needleman–Wunsch algorithm to align protein sequences and build heat-map.', metavar='Str', default='ALL')
    parser.add_argument('-i', '--input', help='Input filename for analise', metavar='FILE', required=True)
    
    args=parser.parse_args()
    
    mis = args.mistake
    W_flag = args.W_output
    in_file = args.input
    zag=args.label
          
    # Open FASTA files
    with open(in_file, "rU") as handle:
        seq_list=[]   
        for record in SeqIO.parse(handle, "fasta"):
            seq_list.append(str(record.seq))
    nseq=len(seq_list)        
    # Work with substitution matrix 
    subs_matrix = MatrixInfo.blosum62
    #Так как blosum62 матрица аминокислотных замен, взяты ак последовательности
    all_letters = set( [key[0] for key in subs_matrix.keys()] )
    subs_matrix.update(  {(letter,'-'):mis for letter in all_letters} )
    subs_matrix.update(  {('-',letter):mis for letter in all_letters} )
    #Матрица для оценок выравниваний
    Z=np.zeros((nseq, nseq))
    #Заполним матрицы с ценой для каждого выравнивания 
    for ali1 in range (0, nseq):
        seq1=seq_list[ali1]
        n=len(seq1)
        for ali2 in range (0, nseq):
            seq2=seq_list[ali2]
            m=len(seq2)
            Wmax=np.zeros((n,m))
    # -5 - штраф за вставку/удаление      
            for j in range(1,m):
                Wmax[0][j]= Wmax[0][j - 1] +mis
                   
            for i in range(1,n):
                Wmax[i][0] = Wmax[i - 1][0] +mis 
                
            for i in range(1, n):
                for j in range(1, m):
                    a1=Wmax[i,j-1]+mis
                    a2=Wmax[i-1,j]+mis
                    key=seq1[i-1],seq2[j-1]
    #Проверка наличия ключа в словаре, если в "треугольной матрице" ключ отсутствует
    #вызовем ключ с обратным следованием элементов
                    if key in subs_matrix:
                        a3=Wmax[i-1,j-1]+subs_matrix[seq1[i-1],seq2[j-1]]
                    else:
                        a3=Wmax[i-1,j-1]+subs_matrix[seq2[j-1],seq1[i-1]]
                    Wmax[i,j]=max([a1,a2,a3])
    #Проверка заполнения каждой матрицы выравнивания и записи в итоговую матрицу 
    #этот блок можно удалить 
            if W_flag:
                print('Wmax',Wmax)
            Z[ali1][ali2]=Wmax[n-1][ m-1]
            print(ali1, ali2, Z[ali1][ali2])
    #Итоговая матрица 
    print(zag)       
    print('Z',Z)
    # Plot heatmaps
    plt.imshow(Z,interpolation="nearest")
    plt.colorbar(orientation='vertical')