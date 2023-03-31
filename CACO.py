import time
import numpy as np
import math
from collections import defaultdict
from tqdm.auto import tqdm
import os

from sklearn.metrics.pairwise import cosine_similarity
from sklearn import preprocessing

from model.clustermethod import *

import argparse
def parse_args():
    parser = argparse.ArgumentParser(description='CACO')
    parser.add_argument('--PPIfile', default='data/string_dataset.txt')
    parser.add_argument('--GOfile', default='data/human_mus_go_reviewed_count.txt')
    parser.add_argument('--Core_threshold', default=0.4)
    parser.add_argument('--output', default='results/temp.txt')
    
    args = parser.parse_args()
    return args

def prepare_data(PPI_file, GO_file):
    relations, protein_id, id_protein, edge_num = load_data(PPI_file)
    protein_num = len(protein_id)
    
    go_id = {}
    with open(GO_file, 'r') as f:
        content = f.readlines()
        for line in tqdm(content, desc='Encode GO_ID'):
            line = line.strip('\n').strip().split()
            protein = line[0]
            go = line[1]
            if go in go_id.keys():
                continue
            else:
                newid = len(go_id)
                go_id[go] = newid
        go_number = len(go_id)
        gocnt_matrix = np.zeros((protein_num, go_number))
        for line in tqdm(content, desc = 'GO count Matrix'):
            line = line.strip('\n').strip().split()
            protein, go, cnt = line[0], line[1], float(line[2])
            if protein in protein_id.keys():
                row = protein_id[protein]
                col = go_id[go]
                gocnt_matrix[row][col] = cnt
    return gocnt_matrix

def main(ppifile, emds, resultfile, lamda):
    relations, protein_id, id_protein, edge_num = load_data(ppifile)

    protein_num = len(protein_id)
    print(f'Protein Number: {protein_num}')
    edge_num = np.sum([len(x) for x in relations.values()])
    print(f'PPI Adj Edges: {edge_num/2}')

    w = np.zeros((protein_num, protein_num))
    w_jcs = cal_second(id_protein, relations, w)
    relations = relation_redu(id_protein, relations, w_jcs)
    
    w = np.zeros((protein_num, protein_num))
    w, relations = construct_net_from_emd(id_protein, relations, emds, w)

    edge_num = np.sum([len(x) for x in relations.values()])
    print(f'Final Edges: {edge_num/2}')

    theta = 1

    Core_PCs = Core_algorithm(id_protein, relations, lamda)

    complexes = Complex_algorithm(Core_PCs, relations, w, theta)

    complexes = Complex_redundancy(complexes)

    save_result(complexes, id_protein, result_file = resultfile)

if __name__=='__main__':
    args = parse_args()
    PPI_file = args.PPIfile
    GO_file = args.GOfile
    lamda = float(args.Core_threshold)
    result_file = args.output
    emds = prepare_data(PPI_file, GO_file)
    main(PPI_file, emds, result_file, lamda)
    