import sys
import config as cfg
import matplotlib.pyplot as plt

def main():

    draw_pie()


def draw_pie():
    # UNIPROT HUMAN ---
    # Pie chart, TRP & disease
    labels = 'disease & TRP human', 'TRP human'
    sizes = [103, 177] # 105 = 208 - 103
    explode = (0.1, 0.0)  # only "explode" the 2nd slice (i.e. 'Hogs')

    fig1, ax1 = plt.subplots()
    ax1.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%',
            shadow=True, startangle=90)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

    plt.savefig(cfg.data['figures'] + '/trp&disease-uniprot_human.png')

    # Pie chart, disease non TRP
    labels = 'disease non TRP human', 'UniProt human'
    sizes = [4510, 1660824]  # 1660824 = 1665334 - 4510
    explode = (0.1, 0.0)  # only "explode" the 2nd slice

    fig1, ax1 = plt.subplots()
    ax1.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%',
            shadow=True, startangle=90)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

    plt.savefig(cfg.data['figures'] + '/disease_nontrp-uniprot_human.png')

    # UNIPROT ALL ORGANISMS ---
    # Pie chart, TRP & disease
    labels = 'disease & TRP', 'TRP'
    sizes = [105, 1406]  # 1406 = 1511 - 105
    explode = (0.1, 0.0)  # only "explode" the 2nd slice

    fig1, ax1 = plt.subplots()
    ax1.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%',
            shadow=True, startangle=90)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

    plt.savefig(cfg.data['figures'] + '/trp&disease-uniprot_all.png')

    # Pie chart, disease non TRP
    labels = 'disease non TRP', 'UniProt'
    sizes = [5042, 225573911]  # 225573911 = 225578953 - 5042
    explode = (0.1, 0.0)  # only "explode" the 2nd slice (i.e. 'Hogs')

    fig1, ax1 = plt.subplots()
    ax1.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%',
            shadow=True, startangle=90)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

    plt.savefig(cfg.data['figures'] + '/disease_nontrp-uniprot_all.png')





if __name__ == '__main__':
    sys.exit(main())
