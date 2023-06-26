import time
import numpy as np
import math
from collections import defaultdict
from tqdm.auto import tqdm,trange
import os
from sklearn.metrics.pairwise import cosine_similarity

def load_data(PPIfile):
    edge_num = 0
    relations = defaultdict(list)
    protein_id, id_protein = {}, {}
    
    with open(PPIfile, 'r') as f:
        content = f.readlines()
        for line in content:
            line = line.strip('\n').split('\t')
            protein_source = line[0]
            protein_destination = line[1]
            
            if protein_source not in protein_id.keys():
                newid = len(protein_id)
                protein_id[protein_source] = newid
                id_protein[newid] = protein_source
            if protein_destination not in protein_id.keys():
                newid = len(protein_id)
                protein_id[protein_destination] = newid
                id_protein[newid] = protein_destination
            
            protein_source_id = protein_id[protein_source]
            protein_destination_id = protein_id[protein_destination]
            if protein_source_id == protein_destination_id:
                continue
            if protein_destination_id not in relations[protein_source_id]:
                relations[protein_source_id].append(protein_destination_id)
                relations[protein_destination_id].append(protein_source_id)
                edge_num += 1
    return relations, protein_id, id_protein, edge_num


def construct_net_from_emd(id_protein, relations, emds, w):
    protein_num = len(id_protein)
    cosine_w = cosine_similarity(emds)
    for u in range(protein_num):
        N_u = set(relations[u])
        for v in N_u:
            w[u][v] = cosine_w[u][v]
    return w, relations

def cal_jcs(id_protein, relations, w):
    protein_num = len(id_protein)
    for u in trange(protein_num, desc='Jaccard'):
        N_u = relations[u]
        for v in range(u+1, protein_num):
            N_v = relations[v]
            if len(N_u)>1 or len(N_v)>1:
                inter = set(N_u) & set(N_v)
                union = set(N_u) | set(N_v)
                if len(inter) < 1:
                    w[u][v] = 0.0
                    w[v][u] = 0.0
                else:
                    w[u][v] = float(len(inter)) / len(union)
                    w[v][u] = float(len(inter)) / len(union)
            else:
                w[u][v] = 0.0
                w[v][u] = 0.0
    return w

def cal_second(id_protein, relations, w):
    protein_num = len(id_protein)
    for u in range(protein_num):
        N_u = relations[u]
        for v in N_u:
            w[u][v] = 1
    second_w = np.dot(w, w)
    
    return second_w

def relation_redu(id_protein, relations, w_hocn):
    protein_num = len(id_protein)
    for u in trange(protein_num, desc='Delete Edges'):
        N_u = relations[u]
        
        new_relation = []
        for v in N_u:
            if w_hocn[u][v] == 0.0:
                continue
            else:
                new_relation.append(v)
        relations[u] = new_relation
    return relations

def cal_similarity(v, u, relations):
    N_v = relations[v]
    N_u = relations[u]
    SN_v = set(N_v)
    SN_v.add(v)
    SN_u = set(N_u)
    SN_u.add(u)
    inter = SN_v & SN_u
    SS_vu = float(len(inter)) / math.sqrt(len(SN_v)*len(SN_u))
    return SS_vu


def Core_algorithm(id_protein, relations, lamda):
    core_pcs = {}
    protein_num = len(id_protein)
    for v in trange(protein_num, desc='Core'):
        N_v = relations[v]
        Core_v = set()
        Core_v.add(v)
        for u in N_v:
            # my 
            SS_vu = cal_similarity(v, u, relations)
            if SS_vu > lamda:
                Core_v.add(u)
        if len(Core_v) >=2 :
            core_pcs[v] = list(Core_v)
    return core_pcs

def cal_sumEcore(core_pc, relations, w):
    sum_Ecore = 0.0
    proteins = sorted(core_pc)
    proteins_set = set(proteins)
    for v in proteins:
        N_v = set(relations[v])
        inter = N_v & proteins_set
        inter = list(inter)
        inter = sorted(inter)
        for u in inter:
            if u > v:
                sum_Ecore += w[v][u]
    return sum_Ecore

def new_attachment_algorithm(core_pc, relations, w):
    attachment = set()
    
    sum_Ecore = cal_sumEcore(core_pc, relations, w)
    w_avg_core = 2.0*sum_Ecore / len(core_pc)
    
    for u in core_pc:
        N_u = set(relations[u])
        attachment = attachment | N_u
    attachment = attachment - set(core_pc)
    
    attachment_nodes_dict_in = {}
    attachment_nodes_dict_out = {}
    
    candidate_attachment_proteins = []
    for p in attachment:
        count = 0
        N_p = relations[p]
        sum_p_core_in = 0.0
        sum_p_core_out = 0.0
        for t in N_p:
            if t in core_pc:
                sum_p_core_in += w[p][t]
                count += 1
            else:
                sum_p_core_out += w[p][t]
        if count>=2:
            candidate_attachment_proteins.append(p)
            attachment_nodes_dict_in[p] = sum_p_core_in
            attachment_nodes_dict_out[p] = sum_p_core_out
    # save overlapping node
    attachment = set()
    if len(candidate_attachment_proteins) > 0:
        for p in candidate_attachment_proteins:
            if (attachment_nodes_dict_in[p] >= w_avg_core):
                attachment.add(p)
    attachment_pc = list(attachment)
    return attachment_pc


def Complex_algorithm(Core_PCs, relations, w):
    complexes = {}
    for v in tqdm(Core_PCs.keys(), desc='Attachment'):
        core_pc = Core_PCs[v]
        #附属蛋白质
        attachment_pc = new_attachment_algorithm(core_pc, relations, w)
        complex_v = list(set(core_pc + attachment_pc))
        complex_v.sort()
        complexes[v] = complex_v
    return complexes

def cal_complex_similarity(cp1, cp2):
    '''
    cp1=cp2 时 sm = 1
    其它情况 sm<1
    '''
    inter = cp1 & cp2
    union = cp1 | cp2
    sm = len(inter) / len(union)
    return sm


def Complex_redundancy(complexes):
    cp_redu = set()
    a = set(complexes.keys())
    b = set(complexes.keys())
    for u in tqdm(a, desc='Removing'):
        if u not in cp_redu:
            for v in b:
                if v not in cp_redu:
                    if v > u:
                        set_u = set(complexes[u])
                        set_v = set(complexes[v])
                        sm = cal_complex_similarity(set_u, set_v)
                        if sm>=1:
                            cp_redu.add(v)
    complexes_new = {}
    for u in a:
        if u in cp_redu:
            continue
        else:
            complexes_new[u] = complexes[u]
    
    temp_dic = {}
    for v,cp in complexes_new.items():
        if len(cp) in temp_dic.keys():
            temp_dic[len(cp)] += 1
        else:
            temp_dic[len(cp)] = 1
    
    return complexes_new


def save_result(complexes, id_protein, result_file):
    with open(result_file, 'w') as f:
        rowid = 0
        for u in complexes.keys():
            cps = complexes[u]
            if len(cps) >= 3:
                rowid += 1
                content = ''
                for pid in cps:
                    content = content + ' ' + id_protein[pid]
                content +='\n'
                f.writelines(content)
#     print('Save Results to file: ', result_file)
