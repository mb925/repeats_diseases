import sys
import config as cfg
import pandas as pd
import numpy
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas
import numpy
import json
import subprocess
import logging
from sklearn.metrics import classification_report


def main():

    # modified_residues()
    # reducefile()
    report()


def reducefile():
        content = pd.read_csv(cfg.data['interactors'] + '/intact.txt', sep='\t')
        columns = ['#ID(s) interactor A', 'ID(s) interactor B',
                   'Feature(s) interactor A', 'Feature(s) interactor B']
        print(len(content), 'total interactions')
        df = [d for d in content[columns].to_numpy() if 'binding' in d[2] or 'binding' in d[3]]
        print(len(df), 'interactions annotated with region')
        fl = open(cfg.data['interactors'] + '/intact_reduced.txt', 'w+')
        fl.write('\t'.join(columns) + '\n')
        for row in df:
            fl.write('\t'.join(row))
            fl.write('\n')

def report():

    def parselocation(loc):
        resrange = []
        try:
            key, val = loc.split(':')
            # TODO: deal with n-nand c-c notations
            if '?' not in val and 'n' not in val and 'c' not in val:
                for value in val.split(','):
                    if '..' in value:
                        resrange += range(int(value.split('-')[0].split('..')[1]) - 1, int(value.split('-')[1].split('..')[0]))
                    else:
                        resrange += range(int(value.split('-')[0])-1, int(value.split('-')[1].split('(')[0]))
        except:
            return resrange
        return resrange

    def preload():

        def splitid(key):
            return [d.split(':')[1].split('-')[0] for d in intact[key] if d != '-']

        intact = pd.read_csv(cfg.data['interactors'] + '/intact_reduced.txt', sep='\t')
        repeatsdb = json.loads(open(cfg.data['repeats'] + '/uniprots.json').read())
        intact_ids = set(splitid('#ID(s) interactor A')).union(
            set(splitid('ID(s) interactor B'))
        )
        repeatsdb_ids = set([d['uniprot_id'] for d in repeatsdb])
        dataset = repeatsdb_ids & intact_ids
        print(len(dataset), 'dataset length')
        # fastafl = '\n'.join(['>' + d['uniprot_id'] + '\n' + d['uniprot_sequence']
        #                      for d in repeatsdb if d['uniprot_id'] in dataset])
        # with open(cfg.data['repeats'] + '/intact_dataset.mfas', 'w+') as fl:
        #     fl.write(fastafl)
        #     fl.close()
        # # print(len(dataset))
        # ## process info in the dataset ##
        # # CD-hit
        # # TODO: to check compatibility with environment
        # cmd = 'cdhit -i ' + cfg.data['repeats'] + '/intact_dataset.mfas -o ' + \
        #       cfg.data['repeats'] + '/intact_dataset_cdhit_100.fa -c 0.7'
        # with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8") as proc:
        #     stdout = proc.stdout.readlines()
        #     stderr = proc.stderr.readlines()
        # if stderr:
        #     logging.error(f'{stderr}. Check again...')
        #     stdout = []

        reduced_dataset = [ln.strip('>') for ln in
                           open(cfg.data['repeats'] + '/intact_dataset_cdhit_100.fa').read().split('\n') if '>' in ln]
        intact = pandas.read_csv(cfg.data['interactors'] + '/intact_reduced.txt', sep='\t')
        finaldic = {u: [] for u in reduced_dataset}
        for index, row in intact.iterrows():
            a = row['#ID(s) interactor A'].split(':')[1].split('-')[0]
            b = row['ID(s) interactor B'].split(':')[1].split('-')[0] if row['ID(s) interactor B'] != '-' else None
            if a in reduced_dataset:
                if 'binding' in row['Feature(s) interactor A']:
                    finaldic[a].append({
                        'interactor': b,
                        'location': row['Feature(s) interactor A']
                    })
            if b and b in reduced_dataset:
                if 'binding' in row['Feature(s) interactor B']:
                    finaldic[b].append({
                        'interactor': a,
                        'location': row['Feature(s) interactor B']
                    })
        repdic = {u['uniprot_id']: u for u in json.loads(open(cfg.data['repeats'] + '/uniprots.json').read())}
        reducedfinaldic = {
            u: {
                'interactions': finaldic[u],
                'repeats': repdic[u]
            } for u in finaldic if finaldic[u]}
        open(cfg.data['repeats'] + '/intact_dataset.json', 'w+').write(json.dumps(reducedfinaldic))

        dataset = json.loads(open(cfg.data['repeats'] + '/intact_dataset.json').read())
        interestingtypes = set(['sufficient binding region', 'necessary binding region',
                                'direct binding region', 'binding-associated region'])
        filtered_dataset = [d for d in dataset if interestingtypes & set(
            [loc.split(':')[0] for inte in dataset[d]['interactions'] for loc in inte['location'].split('|')]
        )]
        for protein in filtered_dataset:
            prange = range(len(dataset[protein]['repeats']['uniprot_sequence']))
            # print(protein, len(prange))
            repeat_array = [1 if index in set([ii
                                               for reg in dataset[protein]['repeats']['repeatsdb_consensus_one'] for ii
                                               in range(
                    reg['start'] - 1, reg['end'])]) else 0
                            for index in prange]
            interaction_array = [1 if index in set(
                [ii for inte in dataset[protein]['interactions'] for loc in inte['location'].split('|')
                 for ii in parselocation(loc)]
            ) else 0 for index in prange]
            dataset[protein]['repeat_array'] = repeat_array
            dataset[protein]['interaction_array'] = interaction_array
        open(cfg.data['repeats'] + '/intact_dataset_arrays.json', 'w+').write(json.dumps(dataset))

    # preload() # comment to run only statistical analysis on precalculate files
    dataset = json.loads(open(cfg.data['repeats'] + '/intact_dataset_arrays.json').read())
    overlaprepeat, overlapinteraction = [], []
    correct_flat, predicted_flat = [], []
    for protein in dataset:
        if 'repeat_array' in dataset[protein] and 'interaction_array' in dataset[protein]:
            repeat_array, interaction_array = dataset[protein]['repeat_array'], dataset[protein]['interaction_array']
            prange = range(len(dataset[protein]['repeats']['uniprot_sequence']))
            # overlap with respect to repeats
            sumoverlap = float(sum([1 for i in prange if repeat_array[i] == 1 and interaction_array[i] == 1]))
            overlaprepeat.append(sumoverlap/sum(repeat_array))
            correct_flat += repeat_array
            if sum(interaction_array):
                overlapinteraction.append(sumoverlap / sum(interaction_array))
                if 1 > (sumoverlap / sum(interaction_array)) > 0.9:
                    print(protein)
            else:
                overlapinteraction.append(0)
            predicted_flat += interaction_array
            # print(classification_report(repeat_array, interaction_array))
    print(classification_report(correct_flat, predicted_flat))
    print(numpy.mean(overlaprepeat), 'repeat')
    print(numpy.mean(overlapinteraction), 'interaction')

def interactors_to_file(file):
    df_interactors = []

    with open(file) as file:
        for line in file:
            if 'IntAct;' in line:
                entry = line.split(';')[1].strip()
                interactors = line.split(';')[2].split('.')[0].strip()

                print(interactors)
                df_interactors.append(
                    {
                        'Entry': entry,
                        'IntAct': interactors
                    }
                )

    df_interactors = pd.DataFrame(df_interactors)
    return df_interactors

def execute_t_test(repeat_distribution, non_repeat_distribution, folder, name, title):
    # plot
    plt.figure(figsize=(10, 7), dpi=300)
    kwargs = dict(hist_kws={'alpha': .6}, kde_kws={'linewidth': 2})
    sns.distplot(repeat_distribution, bins=1000, color="dodgerblue",
                 label="Repeat containing proteins in SwissProt, mean: " + str(
                     round(numpy.mean(repeat_distribution), 2)), **kwargs)
    sns.distplot(non_repeat_distribution, bins=1000, color="orange",
                 label="Proteins in SwissProt without repeats, mean: " + str(
                     round(numpy.mean(non_repeat_distribution), 2)), **kwargs)
    plt.xlim(0, 30)
    extraticks = [13, 18, 21]
    plt.xticks(list(plt.xticks()[0]) + extraticks)
    plt.ylabel("Density")
    plt.xlabel("Number of ".format(title))
    plt.legend()
    plt.savefig(cfg.data[folder] + '/{}.png'.format(name))
    # statistical test
    sttest = stats.ttest_ind(repeat_distribution, non_repeat_distribution)
    # for one tailed, divide the p-value by 2
    pvalue = sttest.pvalue / 2
    outfile = open(cfg.data[folder] + '/{}.txt'.format(name), 'w+')
    outfile.write("Test difference between repeats and background through Student's t test, p-value: " + str(pvalue))

def interactors_mean():
    df_uniprots_reference = cfg.data['interactors'] + '/uniprot_referenceproteome_intact.txt'
    df_repeats = cfg.data['interactors'] + '/repeatsdb_interactors.txt'

    df_repeats_interactors = interactors_to_file(df_repeats)
    df_uniprots_reference_interactors = interactors_to_file(df_uniprots_reference)

    # mean_repeat = df_repeats_interactors['IntAct'].astype(float).mean()
    # mean_non_repeat = df_uniprots_reference_interactors['IntAct'].astype(int).mean()

    repeat_distribution = df_repeats_interactors['IntAct'].to_numpy().astype(np.float)
    non_repeat_distribution = df_uniprots_reference_interactors['IntAct'].to_numpy().astype(np.float)
    execute_t_test(repeat_distribution, non_repeat_distribution, 'interactors', 'interactors_distributions', 'interactors')
    print(df_repeats_interactors)

def PTMs_mean():
    df_uniprots_reference = pd.read_csv(cfg.data['trp_diseases'] + '/uniprot_referenceproteome.tsv', sep='\t')['Post-translational modification']
    df_repeats = pd.read_csv(cfg.data['trp_diseases'] + '/repeatsdb_type_disease.tsv', sep='\t')['Post-translational modification']


    non_repeat_distribution = []
    for el in df_uniprots_reference:
        if type(el) != str:
            non_repeat_distribution.append(0)
            print(el)
        else:
            non_repeat_distribution.append(len(el.split(';')))

    repeat_distribution = []
    for el in df_repeats:
        if type(el) != str:
            repeat_distribution.append(0)
            print(el)
        else:
            repeat_distribution.append(len(el.split(';')))

    # plot
    plt.figure(figsize=(10, 7), dpi=300)
    kwargs = dict(hist_kws={'alpha': .6}, kde_kws={'linewidth': 2})
    sns.distplot(repeat_distribution, bins=10, color="dodgerblue",
                 label="Repeat containing proteins in SwissProt, mean: " + str(
                     round(numpy.mean(repeat_distribution), 2)), **kwargs)
    sns.distplot(non_repeat_distribution, bins=10, color="orange",
                 label="Proteins in SwissProt without repeats, mean: " + str(
                     round(numpy.mean(non_repeat_distribution), 2)), **kwargs)
    plt.xlim(0, 5)
    extraticks = [4, 6, 8]
    plt.xticks(list(plt.xticks()[0]) + extraticks)
    plt.ylabel("Density")
    plt.xlabel("Number of PTMs")
    plt.legend()
    plt.savefig(cfg.data['PTMs'] + '/{}.png'.format('PTMs_distributions'))
    # statistical test
    sttest = stats.ttest_ind(repeat_distribution, non_repeat_distribution)
    # for one tailed, divide the p-value by 2
    pvalue = sttest.pvalue / 2
    outfile = open(cfg.data['PTMs'] + '/{}.txt'.format('PTMs_distributions'), 'w+')
    outfile.write("Test difference between repeats and background through Student's t test, p-value: " + str(pvalue))

def modified_residues():
    df_uniprots_reference = pd.read_csv(cfg.data['trp_diseases'] + '/uniprot_referenceproteome.tsv', sep='\t')[
        'Modified residue']
    df_repeats = pd.read_csv(cfg.data['trp_diseases'] + '/repeatsdb_type_disease.tsv', sep='\t')[
        'Modified residue']

    non_repeat_distribution = []
    for el in df_uniprots_reference:
        if type(el) != str:
            non_repeat_distribution.append(0)
            print(el)
        else:
            non_repeat_distribution.append(el.count("MOD_RES"))
    repeat_distribution = []
    for el in df_repeats:
        if type(el) != str:
            repeat_distribution.append(0)
            print(el)
        else:
            repeat_distribution.append(el.count("MOD_RES"))
    print(non_repeat_distribution)

    execute_t_test(repeat_distribution, non_repeat_distribution, 'modified_residues', 'modified_residues_distribution', 'modified residues')


if __name__ == '__main__':
    sys.exit(main())
