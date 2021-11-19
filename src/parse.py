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

def table_count():
    df_uniprots = pd.read_csv(cfg.data['trp_diseases'] + '/uniprot_type_disease.tsv', sep='\t')
    df_repeats = pd.read_csv(cfg.data['trp_diseases'] + '/repeatsdb_type_disease.tsv', sep='\t')
    df_trp_dis = df_repeats[df_repeats['Involvement in disease'].notna()]


    # TRP
    # filter
    df_pdb_trp = df_repeats['Cross-reference (PDB)'].tolist()
    df_pdb_trp_human = df_repeats.loc[(df_repeats.Organism == 'Homo sapiens (Human)')]['Cross-reference (PDB)'].tolist()
    # count
    pdb_trp = how_many(df_pdb_trp)
    pdb_human_trp = how_many(df_pdb_trp_human)
    uniprot_trp = len(df_repeats['Entry'].index)
    uniprot_human_trp = len(df_repeats.loc[(df_repeats.Organism == 'Homo sapiens (Human)')]['Entry'].index)


    # DISEASE
    # filter
    df_pdb_disease = df_uniprots['Cross-reference (PDB)'].tolist()
    pdb_human_disease = df_uniprots.loc[(df_uniprots.Organism == 'Homo sapiens (Human)')]['Cross-reference (PDB)'].tolist()
    # count
    pdb_disease = how_many(df_pdb_disease)
    pdb_human_disease = how_many(pdb_human_disease)
    uniprot_disease = len(df_uniprots['Entry'].index)
    uniprot_human_disease = len(df_uniprots.loc[(df_uniprots.Organism == 'Homo sapiens (Human)')]['Entry'].index)

    # TRP & DISEASE
    # filter

    pdb_trp_dis = df_trp_dis['Cross-reference (PDB)'].tolist()
    pdb_human_trp_dis = df_trp_dis.loc[(df_repeats.Organism == 'Homo sapiens (Human)')]['Cross-reference (PDB)'].tolist()
    # count
    pdb_trp_dis = how_many(pdb_trp_dis)
    pdb_human_trp_dis = how_many(pdb_human_trp_dis)
    uniprot_trp_dis = len(df_trp_dis['Entry'].index)
    uniprot_human_trp_dis = len(df_trp_dis.loc[(df_trp_dis.Organism == 'Homo sapiens (Human)')]['Entry'].index)

    # DISEASE & NON TRP

    pdb_trp_ndis = pdb_disease - pdb_trp_dis
    pdb_human_trp_ndis = pdb_human_disease - pdb_human_trp_dis
    uniprot_trp_ndis = uniprot_disease - uniprot_trp_dis
    uniprot_human_trp_ndis = uniprot_human_disease - uniprot_human_trp_dis
    rows = [(pdb_trp, pdb_human_trp, uniprot_trp, uniprot_human_trp),
            (pdb_disease, pdb_human_disease, uniprot_disease, uniprot_human_disease),
            (pdb_trp_dis, pdb_human_trp_dis, uniprot_trp_dis, uniprot_human_trp_dis),
            (pdb_trp_ndis, pdb_human_trp_ndis, uniprot_trp_ndis, uniprot_human_trp_ndis)]
    df_count_disease = pd.DataFrame(rows, columns=['PDB', 'PDB Human', 'UniProt', 'UniProt Human'], index=['TRP', 'disease', 'disease & TRP ', 'disease non TRP'])
    print(df_count_disease)

    # todo save the row index
    df_count_disease.to_csv(cfg.data['trp_diseases'] + '/count_disease.tsv', sep='\t', index=False)


if __name__ == '__main__':
    sys.exit(main())
