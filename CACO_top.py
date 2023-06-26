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
    parser.add_argument('--Core_threshold', default=0.4)
    parser.add_argument('--output', default='results/temp.txt')
    
    args = parser.parse_args()
    return args

def main(ppifile, resultfile, lamda):
    relations, protein_id, id_protein, edge_num = load_data(ppifile)

    protein_num = len(protein_id)
    print(f'Protein Number: {protein_num}')
    edge_num = np.sum([len(x) for x in relations.values()])
    print(f'PPI Adj Edges: {edge_num/2}')

    w = np.zeros((protein_num, protein_num))
    w_jcs = cal_second(id_protein, relations, w)
    relations = relation_redu(id_protein, relations, w_jcs)
    
    w = np.zeros((protein_num, protein_num))
    w = cal_jcs(id_protein, relations, w) 

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
    lamda = float(args.Core_threshold)
    result_file = args.output
    result_file = 'results/{0}'.format(PPI_file.split('/')[-1])
    main(PPI_file, result_file, lamda)
    