import pdb, os
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
from utils import *

#-------------------------------------------
def plot(x, y, legends, title='', x_label='', y_label='', file_name='plot'):

    fig1, ax1 = plt.subplots(figsize=(12,8))
    ax1.set_title(title, fontsize=28)

    line_styles = ['-', '--', '-.', ':']

    for i in range(len(y)):
        ax1.plot(x, y[i], c='k', ls=line_styles[i], lw=2)

    ax1.legend(legends, fontsize=20)
    ax1.set_xlabel(x_label, fontsize=25)
    ax1.set_ylabel(y_label, fontsize=25)
    ax1.tick_params(axis='both', which='major', labelsize=20)
    
    plt.savefig(file_name)
    plt.close()

#-------------------------------------------
def plotComparison(x, y, legends, title='', x_label='', y_label='', file_name='plot'):

    fig1, ax1 = plt.subplots(figsize=(20,10))
    ax1.set_title(title, fontsize=28)

    colors = ['r', 'g', 'b', 'k', 'c', 'y']

    for i in range(0, len(y), 2):
        ax1.plot(x, y[i], c=colors[i//2], ls='-', lw=2)
        ax1.plot(x, y[i+1], c=colors[(i+1)//2], ls='--', lw=2)

    ax1.margins(0)
    ax1.legend(legends, bbox_to_anchor=(1, 0.5), loc='center left', fontsize=20)
    ax1.set_xlabel(x_label, fontsize=25)
    ax1.set_ylabel(y_label, fontsize=25)
    ax1.tick_params(axis='both', which='major', labelsize=20)
    plt.subplots_adjust(right=0.7)
    plt.savefig(file_name)
    plt.close()

#-------------------------------------------
def mean_confidence_interval(data, confidence=0.95):
    
    a = 1.0 * np.array(data)
    n = len(a)
    m, std, se = np.mean(a), np.std(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, std, m-h, m+h

# ------------------------------------
def printDec(msg):

    horizontal_border = '_'*50
    vertical_border = '|'

    l = len(horizontal_border)

    print(horizontal_border)
    print(' '*l)
    msg_part = msg
    while len(msg_part) >= l-4:
        print(vertical_border + ' ' + msg_part[:l-4] + ' ' + vertical_border)
        msg_part = msg_part[l-4:]
    print(vertical_border + ' ' + msg_part + ' '*(l-3-len(msg_part)) + vertical_border)
    print(horizontal_border)
    print("")

#-----------------------------------------------
def compareSettings(settings):

    changes = []

    for key in settings_default:
        key1 = key.split(' (')[0]

        assert key1 in settings, '{} not found in parameter settings !!'.format(key1)

        if settings_default[key] != settings[key1]:
            changes.append('{} changed from {} (default) to {}'.format(key1, settings_default[key], settings[key1]))
            
    if len(changes) > 0:
        printDec('\n'.join(changes))
    else:
        printDec('Nothing was changed from default')

#-----------------------------------------------
def analyze(file_name):

    data = [re.sub(r'[ ]+', ' ', item.strip()) for item in readFile(file_name)]
    
    #-----------------------------------------------
    # Get GA settings
    #-----------------------------------------------
    settings = dict([item.split(' : ') for item in data[:16]])

    num_runs = int(settings['Number of Runs'])
    num_gens = int(settings['Generations per Run'])
    size_gene = int(settings['Size of Genes'])

    #-----------------------------------------------
    # Compare the settings of the param file with
    # the default settings, to report what
    # parameters were changed.
    #-----------------------------------------------
    #compareSettings(settings)

    #-----------------------------------------------
    # Get GA data for all gens and all runs
    #-----------------------------------------------
    data = [item.split(' ') for item in data[18:] if item != '']

    # -----------------------------------------------
    # Separate best results from all results
    # -----------------------------------------------
    best_fitness = size_gene
    best_path = data[-(num_gens + 2)][0]
    best_raw_fitness = float(data[-(num_gens + 2)][1])

    # -----------------------------------------------
    # Separate results for all runs
    # -----------------------------------------------
    overall_stats = [[float(item[1]), float(item[2])] for item in data[-num_gens:]]

    assert len(overall_stats) == num_gens, 'Overall stats size ({}) is not equal to num gens ({})'.format(len(overall_stats), num_gens)

    # -----------------------------------------------
    # Get stats for all runs and record the indices
    # that had the best fitness
    # -----------------------------------------------
    best_indices = []

    stats = [[] for i in range(num_runs)]
    
    for counter in range(len(data)):
        if data[counter][0] == 'B':
            best_individual = data[counter+1][0]
            best_fitness = data[counter+1][1]
        item = data[counter]
        if item[0] != 'R':
            continue
        i = int(item[1])-1
        stats[i].append([float(item[4]), float(item[5]), float(item[6])])
        if stats[i][-1][0] == best_fitness:
            best_indices.append(len(stats[i])-1)

    #-----------------------------------------------
    # Average best fitness and standard deviation
    #-----------------------------------------------
    x = range(1, num_gens+1)
    y = [[np.average([stats[i][j][0] for i in range(num_runs)]) for j in range(num_gens)], \
         [np.std([stats[i][j][0] for i in range(num_runs)]) for j in range(num_gens)], \
         [np.average([stats[i][j][1] for i in range(num_runs)]) for j in range(num_gens)], \
         [np.std([stats[i][j][1] for i in range(num_runs)]) for j in range(num_gens)]]
    legends = ['Avg best fitness', 'Std dev (best fitness)', 'Avg avg fitness', 'Std dev (avg fitness)']

    plot_file_name = file_name.split('/')[-1][:-3]
    #plot(x, y, legends, title='', x_label='Number of generations', y_label='Fitness', file_name='{}/{}'.format(results_dir, plot_file_name))
    
    #-----------------------------------------------
    # Average best fitness and standard deviation
    # over all runs
    #-----------------------------------------------
    best_fitness_all_runs = [max([stats[i][j][0] for j in range(num_gens)]) for i in range(num_runs)]
    #print(mean_confidence_interval(best_fitness_all_runs))

    print(best_individual)
    #print('Best individual {}, best fitness {}'.format(best_individual, best_fitness))
    #print('Best fitness: ', max(best_fitness_all_runs))

    if len(best_indices) > 0:
        print('Earliest generation to achieve best fitness: ', best_indices[0])
        print('Best fitness achieving generation stats: ', mean_confidence_interval(best_indices))

    return x, [y[0], y[2]], [legends[0], legends[2]], stats

#-----------------------------------------------
# Parameter file and results path
#-----------------------------------------------
param_file_dir = 'IPD_GA/cap5512.code'
results_dir = 'results'

#-----------------------------------------------
# Get the settings from the default param file
# to compare the settings of other param files,
# to report what parameters were changed.
#-----------------------------------------------
settings_default = getSettings('{}/ipd.params'.format(param_file_dir))

#-----------------------------------------------
# Compare results of different GAs defined by
# different param files
#-----------------------------------------------
y_all, legends_all = [], []

param_file_paths = [f for f in os.listdir('IPD_GA/') if os.path.isfile(os.path.join('IPD_GA', f)) and f.split('_')[0] == 'exp-2' and f.split('_')[-1] == 'summary.txt']

for param_file_path in param_file_paths:

    info = dict([item.split('-') for item in param_file_path.split('_')[:-1]])
    #printDec(param_file_path)
    x, y, l, _ = analyze('IPD_GA/{}'.format(param_file_path))
    y_all += y
    legends_all += ['{} (mem {}, xover {}, mut {}, games {})'.format(item, info['mem'], info['xrate'], info['mrate'], info['ng']) for item in l]

#-----------------------------------------------
# Plot the results of different GAs defined by
# different param files
#-----------------------------------------------
#plotComparison(x, y_all, legends_all, title='', x_label='Number of generations', y_label='Fitness', file_name='{}/comparison'.format(results_dir))
