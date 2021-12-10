import sys
import config as cfg
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats


def main():

    table_count()
    # disease_by_class()
    # oddsratio, pvalue = stats.fisher_exact([[105, 4613], [175, 74425]])
    # print(oddsratio)
    # print(pvalue)


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


def how_many_trp(trp_classes, df_human_reference, df_swissprot, df_human):
    trp_unique_all = trp_classes.drop_duplicates(subset='Entry', keep="first")
    trp_all = len(trp_unique_all.index)
    trp_unique_human = trp_classes.loc[(trp_classes.Organism == 'Homo sapiens (Human)')].drop_duplicates(subset='Entry', keep="first")
    trp_human_reference = trp_unique_human.merge(df_human_reference, on=['Entry'], how='inner')
    trp_human_swissprot = trp_unique_human.merge(df_swissprot, on=['Entry'], how='inner')
    trp_human = trp_unique_human.merge(df_human, on=['Entry'], how='inner')
    trp_human_reference = len(trp_human_reference.index)
    trp_human_swissprot = len(trp_human_swissprot.index)
    trp_human = len(trp_human.index)
    trp_pdb_all = trp_unique_all['Cross-reference (PDB)'].tolist()
    trp_pdb_human = trp_unique_human['Cross-reference (PDB)'].tolist()
    trp_pdb_all = how_many(trp_pdb_all)
    trp_pdb_human = how_many(trp_pdb_human)
    trp_omim = len(trp_unique_all['Involvement in disease'].dropna().index)
    return [trp_all, trp_human_reference, trp_pdb_all, trp_pdb_human, trp_omim, trp_human_swissprot, trp_human]

def merge_left_only(df_uniprot, df_repeats):
    df = pd.merge(df_uniprot, df_repeats, on=['Entry'], how="outer", indicator=True).query(
        '_merge=="left_only"')
    cols = [17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35]
    df.drop(df.columns[cols], axis=1, inplace=True)
    df.columns = df_uniprot.columns
    df_uniprot.update(df)
    return df

def table_count():
    # filtrare i 3 dataset da repeats e repeats_uniprot
    df_uniprots_reference = pd.read_csv(cfg.data['trp_diseases'] + '/uniprot_referenceproteome.tsv', sep='\t')
    df_uniprots_swissprot = pd.read_csv(cfg.data['trp_diseases'] + '/uniprot-homosapiens-swissprot.tsv', sep='\t')
    df_uniprots = pd.read_csv(cfg.data['trp_diseases'] + '/uniprot-advancedseach-homosapiens.tsv', sep='\t')


    df_pdbs = pd.read_csv(cfg.data['trp_diseases'] + '/pdb_chain_uniprot.tsv', sep='\t')
    df_repeats_classification = pd.read_csv(cfg.data['trp_diseases'] + '/repeatsdb_type_classification.tsv', sep='\t', index_col=0).drop_duplicates()
    df_repeats = pd.read_csv(cfg.data['trp_diseases'] + '/repeatsdb_type_disease.tsv', sep='\t')
    df_trp_dis = df_repeats[df_repeats['Involvement in disease'].notna()]

    # remove repeats from uniprot databases
    df_uniprots_reference = merge_left_only(df_uniprots_reference, df_repeats)
    df_uniprots_swissprot = merge_left_only(df_uniprots_swissprot, df_repeats)
    df_uniprots = merge_left_only(df_uniprots, df_repeats)

    # UniProt
    uniprot_pdb = len(df_pdbs['Entry'].drop_duplicates().index)
    df_uniprot_human_reference = df_uniprots_reference.loc[(df_uniprots_reference.Organism == 'Homo sapiens (Human)')]
    uniprot_pdb_human = len(df_uniprot_human_reference.merge(df_pdbs, on=['Entry'], how='inner')['Entry'].drop_duplicates().index)
    df_uniprot_human_swissprot = df_uniprots_swissprot.loc[(df_uniprots_swissprot.Organism == 'Homo sapiens (Human)')]
    df_uniprot_human = df_uniprots.loc[(df_uniprots.Organism == 'Homo sapiens (Human)')]
    uniprot_human_reference = len(df_uniprot_human_reference.index)
    uniprot_human_swissprot = len(df_uniprot_human_swissprot.index)
    uniprot_human = len(df_uniprot_human.index)
    uniprot_pdb_merge = df_uniprot_human_reference.merge(df_pdbs, on=['Entry'], how='inner')
    pdb_human = len(uniprot_pdb_merge.drop_duplicates(subset=['PDB', 'CHAIN'], keep="first").index)

    # TRP
    uniprot_human_trp = len(df_repeats.loc[(df_repeats.Organism == 'Homo sapiens (Human)')]['Entry'].index)

    # DISEASE
    uniprot_all_omim = 5147  # fonte: uniprot search by disease
    uniprot_disease_reference = len(df_uniprots_reference.loc[(df_uniprots_reference.Organism == 'Homo sapiens (Human)')][
                                    'Involvement in disease'].dropna().index)  # info matching with: uniprot search by disease and organism
    uniprot_disease_human = len(df_uniprots.loc[(df_uniprots.Organism == 'Homo sapiens (Human)')][
                                    'Involvement in disease'].dropna().index)
    uniprot_disease_swissprot = len(df_uniprots_swissprot.loc[(df_uniprots_swissprot.Organism == 'Homo sapiens (Human)')][
                                    'Involvement in disease'].dropna().index)

    uniprot_trp_dis = len(df_trp_dis['Entry'].index)
    uniprot_trp_dis_human = len(df_repeats.loc[(df_uniprots_reference.Organism == 'Homo sapiens (Human)')]['Involvement in disease'].dropna().index)


    # TRP CLASSES

    trp_ii_classes = df_repeats_classification.loc[(df_repeats_classification['Classification'] >= 2) & (df_repeats_classification['Classification'] <= 2.9)]
    [trp_ii_all, trp_ii_reference, trp_ii_pdb_all, trp_ii_pdb_human, trp_ii_omim, trp_ii_swissprot, trp_ii_human] = how_many_trp(trp_ii_classes, df_uniprots_reference, df_uniprots_swissprot, df_uniprots)


    trp_iii_classes = df_repeats_classification.loc[(df_repeats_classification['Classification'] >= 3) & (df_repeats_classification['Classification'] <= 3.9)]
    [trp_iii_all, trp_iii_reference, trp_iii_pdb_all, trp_iii_pdb_human, trp_iii_omim, trp_iii_swissprot, trp_iii_human] = how_many_trp(trp_iii_classes, df_uniprots_reference, df_uniprots_swissprot, df_uniprots)


    trp_iv_classes = df_repeats_classification.loc[(df_repeats_classification['Classification'] >= 3) & (df_repeats_classification['Classification'] <= 4.9)]
    [trp_iv_all, trp_iv_reference, trp_iv_pdb_all, trp_iv_pdb_human, trp_iv_omim, trp_iv_swissprot, trp_iv_human] = how_many_trp(trp_iv_classes, df_uniprots_reference, df_uniprots_swissprot, df_uniprots)

    trp_v_classes = df_repeats_classification.loc[(df_repeats_classification['Classification'] >= 5) & (df_repeats_classification['Classification'] <= 5.9)]
    [trp_v_all, trp_v_reference, trp_v_pdb_all, trp_v_pdb_human, trp_v_omim, trp_v_swissprot, trp_v_human] = how_many_trp(trp_v_classes, df_uniprots_reference, df_uniprots_swissprot, df_uniprots)

    pdb_all = len(df_pdbs.drop_duplicates(subset=['PDB', 'CHAIN'], keep="first").index)

    trp_pdb_omim = df_repeats['Cross-reference (PDB)'].tolist()
    trp_pdb_omim = how_many(trp_pdb_omim)

    trp_ii_pdb_omim = trp_ii_classes.drop_duplicates(subset='Entry', keep="first")['Cross-reference (PDB)'].tolist()
    trp_ii_pdb_omim = how_many(trp_ii_pdb_omim)
    trp_iii_pdb_omim = trp_iii_classes.drop_duplicates(subset='Entry', keep="first")['Cross-reference (PDB)'].tolist()
    trp_iii_pdb_omim = how_many(trp_iii_pdb_omim)
    trp_iv_pdb_omim = trp_iv_classes.drop_duplicates(subset='Entry', keep="first")['Cross-reference (PDB)'].tolist()
    trp_iv_pdb_omim = how_many(trp_iv_pdb_omim)
    trp_v_pdb_omim = trp_v_classes.drop_duplicates(subset='Entry', keep="first")['Cross-reference (PDB)'].tolist()
    trp_v_pdb_omim = how_many(trp_v_pdb_omim)


    # TRP CLASSES INTERSECTIONS

    trp_ii_iii = how_many_merge_classes(trp_ii_classes, trp_iii_classes)
    trp_ii_iv = how_many_merge_classes(trp_ii_classes, trp_iv_classes)
    trp_ii_v = how_many_merge_classes(trp_ii_classes, trp_v_classes)
    trp_iii_iv = how_many_merge_classes(trp_iii_classes, trp_iv_classes)
    trp_iii_v = how_many_merge_classes(trp_iii_classes, trp_v_classes)
    trp_iv_v = how_many_merge_classes(trp_iv_classes, trp_v_classes)

    w, h = 17, 17
    matrix = [[0 for x in range(w)] for y in range(h)]
    # PDB all
    matrix[0][0] = pdb_all # PDB all
    matrix[0][1] = pdb_human  # PDB Human
    matrix[0][2] = uniprot_pdb # UniProt all
    matrix[0][3] = uniprot_pdb_human  # UniProt Human reference
    matrix[0][4] = 0 # TRP UniProt
    matrix[0][5] = 0 # OMIM UniProt
    matrix[0][6] = 0 # TRP II
    matrix[0][7] = 0 # TRP III
    matrix[0][8] = 0 # TRP IV
    matrix[0][9] = 0 # TRP V
    matrix[0][10] = trp_ii_pdb_all + trp_iii_pdb_all + trp_iv_pdb_all + trp_v_pdb_all # TRP PDB
    matrix[0][11] = trp_ii_pdb_all # TRP II PDB
    matrix[0][12] = trp_iii_pdb_all # TRP III PDB
    matrix[0][13] = trp_iv_pdb_all # TRP IV PDB
    matrix[0][14] = trp_iv_pdb_all # TRP V PDB
    matrix[0][15] = 0 # UniProt Human SwissProt
    matrix[0][16] = 0 # UniProt Human

    # PDB Human
    matrix[1][1] = 0 # PDB Human
    matrix[1][2] = 0 # UniProt all
    matrix[1][3] = 0 # UniProt Human reference
    matrix[1][4] = 0 # TRP UniProt
    matrix[1][5] = 0 # OMIM UniProt
    matrix[1][6] = 0 # TRP II
    matrix[1][7] = 0 # TRP III
    matrix[1][8] = 0 # TRP IV
    matrix[1][9] = 0 # TRP V
    matrix[1][10] = trp_ii_pdb_human + trp_iii_pdb_human + trp_iv_pdb_human + trp_v_pdb_human # TRP PDB
    matrix[1][11] = trp_ii_pdb_human # TRP II PDB
    matrix[1][12] = trp_iii_pdb_human # TRP III PDB
    matrix[1][13] = trp_iv_pdb_human # TRP IV PDB
    matrix[1][14] = trp_v_pdb_human # TRP V PDB
    matrix[1][15] = 0  # UniProt Human SwissProt
    matrix[1][16] = 0  # UniProt Human

    # UniProt all
    matrix[2][2] = 0 # UniProt all
    matrix[2][3] = uniprot_human_reference # UniProt Human reference
    matrix[2][4] = trp_ii_all + trp_iii_all + trp_iv_all + trp_v_all # TRP UniProt
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
    matrix[2][15] = uniprot_human_swissprot  # UniProt Human SwissProt
    matrix[2][16] = uniprot_human  # UniProt Human

    # 'UniProt Human'
    matrix[3][3] = uniprot_human_reference # UniProt Human reference
    matrix[3][4] = uniprot_human_trp # TRP UniProt
    matrix[3][5] = uniprot_disease_reference # OMIM UniProt
    matrix[3][6] = trp_ii_reference # TRP II UniProt
    matrix[3][7] = trp_iii_reference # TRP III UniProt
    matrix[3][8] = trp_iv_reference # TRP IV UniProt
    matrix[3][9] = trp_v_reference # TRP V UniProt
    matrix[3][10] = 0 # TRP PDB
    matrix[3][11] = 0 # TRP II PDB
    matrix[3][12] = 0 # TRP III PDB
    matrix[3][13] = 0 # TRP IV PDB
    matrix[3][14] = 0 # TRP V PDB
    matrix[3][15] = 0  # UniProt Human SwissProt
    matrix[3][16] = 0  # UniProt Human

    # TRP UniProt
    matrix[4][4] = trp_ii_all + trp_iii_all + trp_iv_all + trp_v_all # TRP UniProt
    matrix[4][5] = uniprot_trp_dis # OMIM UniProt
    matrix[4][6] = trp_ii_all # TRP II
    matrix[4][7] = trp_iii_all # TRP III
    matrix[4][8] = trp_iv_all # TRP IV
    matrix[4][9] = trp_v_all # TRP V
    matrix[4][10] = trp_ii_all + trp_iii_all + trp_iv_all + trp_v_all # TRP PDB
    matrix[4][11] = 0 # TRP II PDB
    matrix[4][12] = 0 # TRP III PDB
    matrix[4][13] = 0 # TRP IV PDB
    matrix[4][14] = 0 # TRP V PDB
    matrix[4][15] = 0  # UniProt Human SwissProt
    matrix[4][16] = 0  # UniProt Human

    # OMIM UniProt
    matrix[5][5] = 0
    matrix[5][6] = trp_ii_omim
    matrix[5][7] = trp_iii_omim
    matrix[5][8] = trp_iv_omim
    matrix[5][9] = trp_v_omim
    matrix[5][10] = trp_pdb_omim # TRP PDB
    matrix[5][11] = trp_ii_pdb_omim # TRP II PDB
    matrix[5][12] = trp_iii_pdb_omim # TRP III PDB
    matrix[5][13] = trp_iv_pdb_omim # TRP IV PDB
    matrix[5][14] = trp_v_pdb_omim # TRP V PDB
    matrix[5][15] = uniprot_disease_swissprot  # UniProt Human SwissProt
    matrix[5][16] = uniprot_disease_human  # UniProt Human

    # TRP II UniProt
    matrix[6][6] = trp_ii_all
    matrix[6][7] = trp_ii_iii
    matrix[6][8] = trp_ii_iv
    matrix[6][9] = trp_ii_v
    matrix[6][10] = 0  # TRP PDB
    matrix[6][11] = 0  # TRP II PDB
    matrix[6][12] = 0  # TRP III PDB
    matrix[6][13] = 0  # TRP IV PDB
    matrix[6][14] = 0  # TRP V PDB
    matrix[6][15] = 0  # UniProt Human SwissProt
    matrix[6][16] = 0  # UniProt Human

    # TRP III UniProt
    matrix[7][7] = trp_iii_all
    matrix[7][8] = trp_iii_iv
    matrix[7][9] = trp_iii_v
    matrix[7][10] = 0  # TRP PDB
    matrix[7][11] = 0  # TRP II PDB
    matrix[7][12] = 0  # TRP III PDB
    matrix[7][13] = 0  # TRP IV PDB
    matrix[7][14] = 0  # TRP V PDB
    matrix[7][15] = 0  # UniProt Human SwissProt
    matrix[7][16] = 0  # UniProt Human

    # TRP IV UniProt
    matrix[8][8] = trp_iv_all
    matrix[8][9] = trp_iv_v
    matrix[8][10] = 0  # TRP PDB
    matrix[8][11] = 0  # TRP II PDB
    matrix[8][12] = 0  # TRP III PDB
    matrix[8][13] = 0  # TRP IV PDB
    matrix[8][14] = 0  # TRP V PDB
    matrix[8][15] = 0  # UniProt Human SwissProt
    matrix[8][16] = 0  # UniProt Human

    # TRP V UniProt
    matrix[9][9] = trp_v_all
    matrix[9][10] = 0  # TRP PDB
    matrix[9][11] = 0  # TRP II PDB
    matrix[9][12] = 0  # TRP III PDB
    matrix[9][13] = 0  # TRP IV PDB
    matrix[9][14] = 0  # TRP V PDB
    matrix[9][15] = 0  # UniProt Human SwissProt
    matrix[9][16] = 0  # UniProt Human

    # TRP PDB
    matrix[10][10] = trp_ii_pdb_all + trp_iii_pdb_all + trp_iv_pdb_all + trp_v_pdb_all  # TRP PDB
    matrix[10][11] = trp_ii_pdb_all  # TRP II PDB
    matrix[10][12] = trp_iii_pdb_all  # TRP III PDB
    matrix[10][13] = trp_iv_pdb_all  # TRP IV PDB
    matrix[10][14] = trp_v_pdb_all  # TRP V PDB
    matrix[10][15] = 0  # UniProt Human SwissProt
    matrix[10][16] = 0  # UniProt Human

    # TRP II PDB
    matrix[11][11] = trp_ii_pdb_all  # TRP II PDB
    matrix[11][12] = trp_iii_pdb_all  # TRP III PDB
    matrix[11][13] = trp_iv_pdb_all  # TRP IV PDB
    matrix[11][14] = trp_v_pdb_all  # TRP V PDB
    matrix[11][15] = trp_ii_swissprot  # UniProt Human SwissProt
    matrix[11][16] = trp_ii_human  # UniProt Human

    # TRP III PDB
    matrix[12][12] = trp_iii_pdb_all  # TRP III PDB
    matrix[12][13] = trp_iv_pdb_all  # TRP IV PDB
    matrix[12][14] = trp_v_pdb_all  # TRP V PDB
    matrix[12][15] = trp_iii_swissprot  # UniProt Human SwissProt
    matrix[12][16] = trp_iii_human  # UniProt Human

    # TRP IV PDB
    matrix[13][13] = trp_iv_pdb_all  # TRP IV PDB
    matrix[13][14] = trp_v_pdb_all  # TRP V PDB
    matrix[13][15] = trp_iv_swissprot  # UniProt Human SwissProt
    matrix[13][16] = trp_iv_human  # UniProt Human

    # TRP V PDB
    matrix[14][14] = trp_v_pdb_all  # TRP V PDB
    matrix[14][15] = trp_v_swissprot  # UniProt Human SwissProt
    matrix[14][16] = trp_v_human  # UniProt Human

    # UniProt Human Swissprot
    matrix[15][15] = uniprot_human_swissprot  # UniProt Human SwissProt
    matrix[15][16] = 0  # UniProt Human

    # # UniProt Human
    matrix[16][16] = uniprot_human  # UniProt Human

    # # filling the lower triangle of the matrix (copying upper)
    # for i in range(17):
    #     for j in range(i, 17):
    #         matrix[j][i] = matrix[i][j]

    label = ['PDB all', 'PDB Human', 'UniProt all', 'UniProt Human reference', 'TRP UniProt', 'OMIM UniProt', 'TRP II', 'TRP III', 'TRP IV', 'TRP V', 'TRP PDB', 'TRP II PDB', 'TRP III PDB', 'TRP IV PDB', 'TRP V PDB', 'UniProt Swissprot', 'Uniprot Human']
    df_count_disease = pd.DataFrame(matrix, columns=label)
    df_count_disease.index = label
    # print(df_count_disease)

    df_count_disease.to_csv(cfg.data['trp_diseases'] + '/contingency_table.tsv', sep='\t')


def disease_by_class():
    # Create dataset
    height = [1, 59, 84, 26]
    bars = ('II', 'III', 'IV', 'V')
    x_pos = np.arange(len(bars))

    # Create bars
    plt.bar(x_pos, height, color=['#6DA6C5'])

    # Create names on the x-axis
    plt.xticks(x_pos, bars)



    # Show graphic
    plt.savefig(cfg.data['trp_diseases'] + '/{}.png'.format('disease_by_class'))


if __name__ == '__main__':
    sys.exit(main())
