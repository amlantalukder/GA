import pdb, re, os
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
from utils import *

# -------------------------------------------
def plot(x, y, legends, title='', x_label='', y_label='', file_name='plot', errorbar=True):
    fig1, ax1 = plt.subplots(figsize=(12, 8))
    ax1.set_title(title, fontsize=28)

    line_styles = ['-', ':']
    colors = ['k', 'b', 'r', 'g', 'y']

    for i in range(len(y)):
        if errorbar:
            ax1.errorbar(x, y[i][0], yerr=y[i][1], c=colors[i], ls=line_styles[0], lw=2)
        else:
            ax1.plot(x, y[i][0], c=colors[i], ls=line_styles[0], lw=2)

    ax1.legend(legends, fontsize=20)
    ax1.set_xlabel(x_label, fontsize=25)
    ax1.set_ylabel(y_label, fontsize=25)
    ax1.tick_params(axis='both', which='major', labelsize=20)

    plt.savefig(file_name)
    plt.close()

# -------------------------------------------
def plotComparison(x, y, legends, title='', x_label='', y_label='', file_name='plot'):
    fig1, ax1 = plt.subplots(figsize=(20, 10))
    ax1.set_title(title, fontsize=28)

    colors = ['g', 'b']
    line_styles = ['-', ':']

    for i in range(len(y)):
        ax1.errorbar(x, y[i][0], yerr=y[i][1], c=colors[i], ls=line_styles[0], lw=2)
        ax1.errorbar(x, y[i][2], yerr=y[i][3], c=colors[i], ls=line_styles[1], lw=2)

    ax1.margins(0)
    ax1.legend(legends, bbox_to_anchor=(1, 0.5), loc='center left', fontsize=20)
    ax1.set_xlabel(x_label, fontsize=25)
    ax1.set_ylabel(y_label, fontsize=25)
    ax1.tick_params(axis='both', which='major', labelsize=20)
    plt.subplots_adjust(right=0.7)
    plt.savefig(file_name)
    plt.close()

# -------------------------------------------
def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, std, se = np.mean(a), np.std(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n - 1)
    return m, std, m - h, m + h


# ------------------------------------
def printDec(msg):
    horizontal_border = '_' * 50
    vertical_border = '|'

    l = len(horizontal_border)

    print(horizontal_border)
    print(' ' * l)
    msg_part = msg
    while len(msg_part) >= l - 4:
        print(vertical_border + ' ' + msg_part[:l - 4] + ' ' + vertical_border)
        msg_part = msg_part[l - 4:]
    print(vertical_border + ' ' + msg_part + ' ' * (l - 3 - len(msg_part)) + vertical_border)
    print(horizontal_border)
    print("")


# -----------------------------------------------
def analyze(file_name, best_fitness = None):
    data = readFile(file_name)
    data = (''.join(data)).split('\n\n')

    stats = [[[float(item2) for item2 in item1.strip().split(',')[2:]] for item1 in item.split('\n')] for item in data[:-1]]
    best_stats = [[float(item1) for item1 in item.strip().split(',')[1:]] for item in data[-1].split('\n')]

    num_runs = len(stats)
    num_gens = len(stats[0])

    # -----------------------------------------------
    # Average best fitness and standard deviation
    # -----------------------------------------------
    x = range(1, num_gens + 1)
    y = [[[np.average([stats[i][j][0] for i in range(num_runs)]) for j in range(num_gens)], \
         [np.std([stats[i][j][0] for i in range(num_runs)]) for j in range(num_gens)]], \
         [[np.average([stats[i][j][1] for i in range(num_runs)]) for j in range(num_gens)], \
         [np.std([stats[i][j][1] for i in range(num_runs)]) for j in range(num_gens)]]]
    legends = ['Avg best fitness', 'Std dev (best fitness)', 'Avg avg fitness', 'Std dev (avg fitness)']
    plot_file_name = file_name.split('/')[-1][:-3]
    plot(x, y, legends, title='', x_label='Number of generations', y_label='Fitness', file_name='{}/{}'.format(results_dir, plot_file_name))

    # -----------------------------------------------
    # Average best fitness and standard deviation
    # over all runs
    # -----------------------------------------------
    best_fitness_all_runs = [max([stats[i][j][0] for j in range(num_gens)]) for i in range(num_runs)]
    print(mean_confidence_interval(best_fitness_all_runs))

    # print(best_individual)
    # print('Best individual {}, best fitness {}'.format(best_individual, best_fitness))
    # print('Best fitness: ', max(best_fitness_all_runs))

    if best_fitness != None:
        best_indices = []
        for stats_gen in stats:
            for i in range(len(stats_gen)):
                if stats_gen[i][0] == best_fitness:
                    best_indices.append(i+1)
                    break

        if len(best_indices) > 0:
            print('Earliest generation to achieve best fitness: ', best_indices[0])
            print('Best fitness achieving generation stats: ', mean_confidence_interval(best_indices))

    return x, y, legends


# -----------------------------------------------
# Results path
# -----------------------------------------------
results_dir = 'results'

ops_credit_probs = readDataTable('{}/operator_probabilities.csv'.format(results_dir))

data = []
for i in range(1, len(ops_credit_probs)):
    op_wise_score = {}
    for j in range(len(ops_credit_probs[i])):
        optype = ops_credit_probs[0][j]
        if optype not in op_wise_score:
            op_wise_score[optype] = [float(ops_credit_probs[i][j])]
        else:
            op_wise_score[optype].append(float(ops_credit_probs[i][j]))

    op_wise_score = dict([[optype, [np.average(op_wise_score[optype]), np.std(op_wise_score[optype])]] for optype in op_wise_score])
    data.append(op_wise_score)

x = range(1, len(data)+1)
legends = sorted(data[0].keys())
y = [[[item[optype][0] for item in data], [item[optype][1] for item in data]] for optype in legends]
plot(x, y, legends, title='', x_label='New chromosomes', y_label='Probabilities', file_name='{}/operator_evolution'.format(results_dir), errorbar=False)


# -----------------------------------------------
# Compare results of different GAs defined by
# different param files
# -----------------------------------------------
y_all, legends_all = [], []

summary_file_paths = [f for f in os.listdir(results_dir) \
                 if os.path.isfile(os.path.join(results_dir, f)) and f.split('_')[-1] == 'summary.csv']

for summary_file_path in summary_file_paths:
    info = [item.split('-') for item in summary_file_path.split('_')[1:-1]]
    x, y, l = analyze('{}/{}'.format(results_dir, summary_file_path))
    y_all += y
    legends_all += ['{} ({})'.format(item, formatDataTable(info, ' ', ',')) for item in l]

# -----------------------------------------------
# Plot the results of different GAs defined by
# different param files
# -----------------------------------------------
#plotComparison(x, y_all, legends_all, title='', \
 #              x_label='Number of generations', y_label='Fitness', \
  #             file_name='{}/comparison'.format(results_dir))