import pdb, os
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
from utils import *

#-------------------------------------------
def plot(x, y, legends, title='', x_label='', y_label='', file_name='plot'):

    fig1, ax1 = plt.subplots(figsize=(12,8))
    ax1.set_title(title, fontsize=28)

    line_styles = ['-', ':']

    ax1.errorbar(x, y[0], yerr=y[1], c='k', ls=line_styles[0], lw=2)
    ax1.errorbar(x, y[2], yerr=y[3], c='k', ls=line_styles[1], lw=2)

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

#-------------------------------------------
def mean_confidence_interval(data, confidence=0.95):
    
    a = 1.0 * np.array(data)
    n = len(a)
    m, std, se = np.mean(a), np.std(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, std, m-h, m+h

# ------------------------------------
def printDec(msg):
    horizontal_border = '_' * 50
    vertical_border = '|'

    l = len(horizontal_border)

    print(horizontal_border)
    print(' ' * l)
    msg_part = msg.strip()
    while len(msg_part) >= l - 4:
        print(vertical_border + ' ' + msg_part[:l - 4] + ' ' + vertical_border)
        msg_part = msg_part[l - 4:].strip()
    print(vertical_border + ' ' + msg_part + ' ' * (l - 3 - len(msg_part)) + vertical_border)
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
    optimum_fitness = None
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

    gen_counter = 0
    for item in data:
        if item[0] != 'R':
            continue
        i = int(item[1])-1
        try:
            stats[i].append([float(item[4]), float(item[5]), float(item[6])])
        except:
            pdb.set_trace()

        if optimum_fitness and stats[i][-1][0] == optimum_fitness:
            best_indices.append(len(stats[i])-1)

    #-----------------------------------------------
    # Average best fitness and standard deviation
    #-----------------------------------------------
    x = range(1, num_gens+1)
    y = [[np.average([stats[i][j][0] for i in range(num_runs)]) for j in range(num_gens)], \
         [np.std([stats[i][j][0] for i in range(num_runs)]) for j in range(num_gens)], \
         [np.average([stats[i][j][1] for i in range(num_runs)]) for j in range(num_gens)], \
         [np.std([stats[i][j][1] for i in range(num_runs)]) for j in range(num_gens)]]
    legends = ['Avg best fitness', 'Avg avg fitness']

    #-----------------------------------------------
    # Average best fitness and standard deviation
    # over all runs
    #-----------------------------------------------
    best_fitness_all_runs = [min([stats[i][j][0] for j in range(num_gens)]) for i in range(num_runs)]
    conf_int = mean_confidence_interval(best_fitness_all_runs)
    best_fitness = min(best_fitness_all_runs)
    best_fitness_stat = [best_fitness] + list(conf_int)

    #print('Best fitness: ', best_fitness)

    if len(best_indices) > 0:
        print('Earliest generation to achieve best fitness: ', best_indices[0])
        print('Best fitness achieving generation stats: ', mean_confidence_interval(best_indices))

    return x, y, legends, best_fitness_stat

#-----------------------------------------------
header = ['Rep', 'Data set', 'Best Fitness', 'Average Best Fitness', 'Std Dev Best Fitness', '95% confidence interval', 'Param File']

reprs = ['Representation_1', 'Representation_2']
repr_suff = ['Rep 1', 'Rep 2']

plot_compare_info = {}

for i in range(len(reprs)):

    #-----------------------------------------------
    # Parameter file and results path
    #-----------------------------------------------
    param_file_dir = 'codes/{}'.format(reprs[i])
    results_dir = 'results'

    #-----------------------------------------------
    # Get the settings from the default param file
    # to compare the settings of other param files,
    # to report what parameters were changed.
    #-----------------------------------------------
    #settings_default = getDefaultSettings('{}/TSP.params'.format(param_file_dir))

    #-----------------------------------------------
    # Run analysis on results of GA run by different
    # settings.
    #-----------------------------------------------
    param_files_all = [item for item in os.listdir(param_file_dir) if os.path.isfile(os.path.join(param_file_dir, item)) and item.split('.')[-2] == 'params_summary']
    data_sets = ['att48', 'berlin52', 'rl1323']

    for d in data_sets:

        printDec(d)

        perc = 10
        counter = 0

        stats = []
        plot_info = {}

        for param_file in param_files_all:
            counter += 1
            perc = showPercBar(counter, len(param_files_all), perc)
            if d not in param_file:
                continue
            x, y, l, b = analyze(os.path.join(param_file_dir, param_file))

            plot_info[param_file] = [x, y, l]

            stats.append([reprs[i], d] + b + [param_file])

        stats = sorted(stats)
        print(formatDataTable([header] + stats[:10]))
        best_param_file = stats[0][-1]
        x, y, l = plot_info[stats[0][-1]]
        plot_file_name = best_param_file.split('/')[-1][:-3]
        plot(x, y, l, title='', x_label='Number of generations', y_label='Fitness', file_name='{}/{}'.format(results_dir, plot_file_name))

        if d not in plot_compare_info:
            plot_compare_info[d] = [x, [], []]
        plot_compare_info[d][1].append(y)
        plot_compare_info[d][2] += ['{} ({})'.format(item, repr_suff[i]) for item in l]

#-----------------------------------------------
# Plot the results of different GAs defined by
# different param files
#-----------------------------------------------
for d in data_sets:
    x, y_all, legends_all = plot_compare_info[d]
    plotComparison(x, y_all, legends_all, title='', x_label='Number of generations', y_label='Fitness', file_name='{}/comparison_{}'.format(results_dir, d))
