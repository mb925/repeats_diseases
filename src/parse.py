import sys
import config as cfg
import pandas as pd
import math


def main():

    table_count()


def how_many(elements):
    list_pdb = []
    for uniprot in elements:
        if isinstance(uniprot, str):
            for pdb in uniprot.split(';'):
                if pdb not in list_pdb and pdb != '':
                    list_pdb.append(pdb)
    return len(list_pdb)

def how_many_merge_classes(trp_class1, trp_class2):
    trp_unique_1 = trp_class1['Entry'].drop_duplicates().to_frame()
    trp_unique_2 = trp_class2['Entry'].drop_duplicates().to_frame()
    trp_unique_all = trp_unique_1.merge(trp_unique_2, how='inner', on='Entry')
    return len(trp_unique_all.index)


def how_many_trp(trp_classes):
    trp_unique_all = trp_classes.drop_duplicates(subset='Entry', keep="first")
    trp_all = len(trp_unique_all.index)
    trp_unique_human = trp_classes.loc[(trp_classes.Organism == 'Homo sapiens (Human)')].drop_duplicates(subset='Entry', keep="first")
    trp_human = len(trp_unique_human.index)
    trp_pdb_all = trp_unique_all['Cross-reference (PDB)'].tolist()
    trp_pdb_human = trp_unique_human['Cross-reference (PDB)'].tolist()
    trp_pdb_all = how_many(trp_pdb_all)
    trp_pdb_human = how_many(trp_pdb_human)
    trp_omim = len(trp_unique_all['Involvement in disease'].dropna().index)
    return [trp_all, trp_human, trp_pdb_all, trp_pdb_human, trp_omim]

def table_count():
    df_uniprots = pd.read_csv(cfg.data['trp_diseases'] + '/uniprot_referenceproteome.tsv', sep='\t')
    df_pdbs = pd.read_csv(cfg.data['trp_diseases'] + '/pdb_chain_uniprot.tsv', sep='\t')
    df_repeats_classification = pd.read_csv(cfg.data['trp_diseases'] + '/repeatsdb_type_classification.tsv', sep='\t',  index_col=0).drop_duplicates()
    df_repeats = pd.read_csv(cfg.data['trp_diseases'] + '/repeatsdb_type_disease.tsv', sep='\t')
    df_trp_dis = df_repeats[df_repeats['Involvement in disease'].notna()]

    uniprot_all = 225578953
    uniprot_human = 79038 # fonte: uniprot search by human organism

    # TRP
    # filter
    df_pdb_trp = df_repeats['Cross-reference (PDB)'].tolist()
    df_pdb_trp_human = df_repeats.loc[(df_repeats.Organism == 'Homo sapiens (Human)')]['Cross-reference (PDB)'].tolist()
    # count
    pdb_trp = how_many(df_pdb_trp)
    pdb_human_trp = how_many(df_pdb_trp_human)
    uniprot_trp = len(df_repeats['Entry'].index) ## duplicated lines because of the different classification
    uniprot_human_trp = len(df_repeats.loc[(df_repeats.Organism == 'Homo sapiens (Human)')]['Entry'].index)



    # DISEASE
    # filter
    pdb_human_disease = df_uniprots.loc[(df_uniprots.Organism == 'Homo sapiens (Human)')]['Cross-reference (PDB)'].tolist()
    # count
    pdb_human_disease = how_many(pdb_human_disease)
    uniprot_disease_human = len(df_uniprots.loc[(df_uniprots.Organism == 'Homo sapiens (Human)')]['Involvement in disease'].dropna().index)  # info matching with: uniprot search by disease and organism

    uniprot_trp_dis = len(df_trp_dis['Entry'].index)
    uniprot_all_omim = 5147 # fonte: uniprot search by disease

    pdb_all = len(df_pdbs.drop_duplicates(subset=['PDB', 'CHAIN'], keep="first").index)

    # TRP CLASSES

    trp_ii_classes = df_repeats_classification.loc[(df_repeats_classification['Classification'] >= 2) & (df_repeats_classification['Classification'] <= 2.9)]
    trps_classes = how_many_trp(trp_ii_classes)
    trp_ii_all = trps_classes[0]
    trp_ii_human = trps_classes[1]
    trp_ii_pdb_all = trps_classes[2]
    trp_ii_pdb_human = trps_classes[3]
    trp_ii_omim = trps_classes[4]

    trp_iii_classes = df_repeats_classification.loc[(df_repeats_classification['Classification'] >= 3) & (df_repeats_classification['Classification'] <= 3.9)]
    trps_classes = how_many_trp(trp_iii_classes)
    trp_iii_all = trps_classes[0]
    trp_iii_human = trps_classes[1]
    trp_iii_pdb_all = trps_classes[2]
    trp_iii_pdb_human = trps_classes[3]
    trp_iii_omim = trps_classes[4]


    trp_iv_classes = df_repeats_classification.loc[(df_repeats_classification['Classification'] >= 3) & (df_repeats_classification['Classification'] <= 4.9)]
    trps_classes = how_many_trp(trp_iv_classes)
    trp_iv_all = trps_classes[0]
    trp_iv_human = trps_classes[1]
    trp_iv_pdb_all = trps_classes[2]
    trp_iv_pdb_human = trps_classes[3]
    trp_iv_omim = trps_classes[4]


    trp_v_classes = df_repeats_classification.loc[(df_repeats_classification['Classification'] >= 5) & (df_repeats_classification['Classification'] <= 5.9)]
    trps_classes = how_many_trp(trp_v_classes)
    trp_v_all = trps_classes[0]
    trp_v_human = trps_classes[1]
    trp_v_pdb_all = trps_classes[2]
    trp_v_pdb_human = trps_classes[3]
    trp_v_omim = trps_classes[4]

    # TRP CLASSES INTERSECTIONS

    trp_ii_iii = how_many_merge_classes(trp_ii_classes, trp_iii_classes)
    trp_ii_iv = how_many_merge_classes(trp_ii_classes, trp_iv_classes)
    trp_ii_v = how_many_merge_classes(trp_ii_classes, trp_v_classes)
    trp_iii_iv = how_many_merge_classes(trp_iii_classes, trp_iv_classes)
    trp_iii_v = how_many_merge_classes(trp_iii_classes, trp_v_classes)
    trp_iv_v = how_many_merge_classes(trp_iv_classes, trp_v_classes)

    w, h = 15, 15
    matrix = [[0 for x in range(w)] for y in range(h)]
    # PDB all
    matrix[0][0] = pdb_all # PDB all
    matrix[0][1] = 0  # PDB Human
    matrix[0][2] = pdb_all # UniProt all
    matrix[0][3] = 0  # 'UniProt Human'
    matrix[0][4] = pdb_trp # TRP UniProt
    matrix[0][5] = 0 # OMIM UniProt
    matrix[0][6] = trp_ii_pdb_all # TRP II
    matrix[0][7] = trp_iii_pdb_all # TRP III
    matrix[0][8] = trp_iv_pdb_all # TRP IV
    matrix[0][9] = trp_v_pdb_all # TRP V
    matrix[0][10] = 0 # TRP PDB
    matrix[0][11] = 0 # TRP II PDB
    matrix[0][12] = 0 # TRP III PDB
    matrix[0][13] = 0 # TRP IV PDB
    matrix[0][14] = 0 # TRP V PDB

    # PDB Human
    matrix[1][1] = 0 # PDB Human
    matrix[1][2] = 0 # UniProt all
    matrix[1][3] = 0 # 'UniProt Human'
    matrix[1][4] = pdb_human_trp # TRP UniProt
    matrix[1][5] = pdb_human_disease # OMIM UniProt
    matrix[1][6] = trp_ii_pdb_human # TRP II
    matrix[1][7] = trp_iii_pdb_human # TRP III
    matrix[1][8] = trp_iv_pdb_human # TRP IV
    matrix[1][9] = trp_v_pdb_human # TRP V
    matrix[1][10] = 0 # TRP PDB
    matrix[1][11] = 0 # TRP II PDB
    matrix[1][12] = 0 # TRP III PDB
    matrix[1][13] = 0 # TRP IV PDB
    matrix[1][14] = 0 # TRP V PDB

    # UniProt all
    matrix[2][2] = uniprot_all # UniProt all
    matrix[2][3] = uniprot_human # 'UniProt Human'
    matrix[2][4] = uniprot_trp # TRP UniProt
    matrix[2][5] = uniprot_all_omim # OMIM UniProt
    matrix[2][6] = trp_ii_all # TRP II
    matrix[2][7] = trp_iii_all # TRP III
    matrix[2][8] = trp_iv_all # TRP IV
    matrix[2][9] = trp_v_all # TRP V
    matrix[2][10] = 0 # TRP PDB
    matrix[2][11] = 0 # TRP II PDB
    matrix[2][12] = 0 # TRP III PDB
    matrix[2][13] = 0 # TRP IV PDB
    matrix[2][14] = 0 # TRP V PDB

    # 'UniProt Human'
    matrix[3][3] = uniprot_human # 'UniProt Human'
    matrix[3][4] = uniprot_human_trp # TRP UniProt
    matrix[3][5] = uniprot_disease_human # OMIM UniProt
    matrix[3][6] = trp_ii_human # TRP II
    matrix[3][7] = trp_iii_human # TRP III
    matrix[3][8] = trp_v_human # TRP IV
    matrix[3][9] = trp_v_human # TRP V
    matrix[3][10] = 0 # TRP PDB
    matrix[3][11] = 0 # TRP II PDB
    matrix[3][12] = 0 # TRP III PDB
    matrix[3][13] = 0 # TRP IV PDB
    matrix[3][14] = 0 # TRP V PDB

    # TRP UniProt
    matrix[4][4] = uniprot_trp # TRP UniProt
    matrix[4][5] = uniprot_trp_dis # OMIM UniProt
    matrix[4][6] = trp_ii_all # TRP II
    matrix[4][7] = trp_iii_all
    matrix[4][8] = trp_iv_all
    matrix[4][9] = trp_v_all
    matrix[4][10] = 0 # TRP PDB
    matrix[4][11] = 0 # TRP II PDB
    matrix[4][12] = 0 # TRP III PDB
    matrix[4][13] = 0 # TRP IV PDB
    matrix[4][14] = 0 # TRP V PDB

    # OMIM UniProt
    matrix[5][5] = 0
    matrix[5][6] = trp_ii_omim
    matrix[5][7] = trp_iii_omim
    matrix[5][8] = trp_iv_omim
    matrix[5][9] = trp_v_omim
    matrix[5][10] = 0 # TRP PDB
    matrix[5][11] = 0 # TRP II PDB
    matrix[5][12] = 0 # TRP III PDB
    matrix[5][13] = 0 # TRP IV PDB
    matrix[5][14] = 0 # TRP V PDB
    # TRP II
    matrix[6][6] = trp_ii_all
    matrix[6][7] = trp_ii_iii
    matrix[6][8] = trp_ii_iv
    matrix[6][9] = trp_ii_v
    matrix[6][10] = 0  # TRP PDB
    matrix[6][11] = 0  # TRP II PDB
    matrix[6][12] = 0  # TRP III PDB
    matrix[6][13] = 0  # TRP IV PDB
    matrix[6][14] = 0  # TRP V PDB

    # TRP III
    matrix[7][7] = trp_iii_all
    matrix[7][8] = trp_iii_iv
    matrix[7][9] = trp_iii_v
    matrix[7][10] = 0  # TRP PDB
    matrix[7][11] = 0  # TRP II PDB
    matrix[7][12] = 0  # TRP III PDB
    matrix[7][13] = 0  # TRP IV PDB
    matrix[7][14] = 0  # TRP V PDB

    # TRP IV
    matrix[8][8] = trp_iv_all
    matrix[8][9] = trp_iv_v
    matrix[8][10] = 0  # TRP PDB
    matrix[8][11] = 0  # TRP II PDB
    matrix[8][12] = 0  # TRP III PDB
    matrix[8][13] = 0  # TRP IV PDB
    matrix[8][14] = 0  # TRP V PDB

    # TRP V
    matrix[9][9] = trp_v_all
    matrix[9][10] = 0  # TRP PDB
    matrix[9][11] = 0  # TRP II PDB
    matrix[9][12] = 0  # TRP III PDB
    matrix[9][13] = 0  # TRP IV PDB
    matrix[9][14] = 0  # TRP V PDB

    # TRP PDB
    matrix[10][10] = 0  # TRP PDB
    matrix[10][11] = 0  # TRP II PDB
    matrix[10][12] = 0  # TRP III PDB
    matrix[10][13] = 0  # TRP IV PDB
    matrix[10][14] = 0  # TRP V PDB

    # TRP I PDB
    matrix[11][11] = 0  # TRP II PDB
    matrix[11][12] = 0  # TRP III PDB
    matrix[11][13] = 0  # TRP IV PDB
    matrix[11][14] = 0  # TRP V PDB

    # TRP II PDB
    matrix[12][12] = 0  # TRP III PDB
    matrix[12][13] = 0  # TRP IV PDB
    matrix[12][14] = 0  # TRP V PDB

    # TRP III PDB
    matrix[13][13] = 0  # TRP IV PDB
    matrix[13][14] = 0  # TRP V PDB

    # TRP IV PDB
    matrix[14][14] = 0  # TRP V PDB

    # filling the lower triangle of the matrix (copying upper)
    for i in range(10):
        for j in range(i, 10):
            matrix[j][i] = matrix[i][j]

    label = ['PDB all', 'PDB Human', 'UniProt all', 'UniProt Human', 'TRP UniProt', 'OMIM UniProt', 'TRP II', 'TRP III', 'TRP IV', 'TRP V', 'TRP PDB', 'TRP II PDB', 'TRP III PDB', 'TRP IV PDB', 'TRP V PDB']
    df_count_disease = pd.DataFrame(matrix, columns=label)
    df_count_disease.index = label
    print(df_count_disease)

    df_count_disease.to_csv(cfg.data['trp_diseases'] + '/count.tsv', sep='\t')


if __name__ == '__main__':
    sys.exit(main())
