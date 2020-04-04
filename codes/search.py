from utils import *
import pdb

# -----------------------------------------------
class Chromosome:

    # -----------------------------------------------
    def __init__(self):

        self.raw_fitness = -1
        self.scl_fitness = -1
        self.pro_fitness = -1
        
        self.chromo = ''
        for i in range(num_genes):
            for j in range(gene_size):
                self.chromo += '0' if random.random() > 0.5 else '1'

    # -----------------------------------------------
    def mutate(self):

        assert mut_type == 1, 'Error: Invalid mutation type {}'.format(mut_type)

        x = ''
        for i in range(num_genes * gene_size):
            if random.random() < mut_rate:
                x += '0' if self.chromo[i] == '1' else '0'
            else:
                x += self.chromo[i]
        self.chromo = x

    # -----------------------------------------------
    def calcFitness(self):

        self.raw_fitness = sum([int(i) for i in self.chromo])

    # -----------------------------------------------
    def clone(self):

        a = Chromosome()
        a.chromo = self.chromo
        a.raw_fitness = self.raw_fitness
        a.scl_fitness = self.scl_fitness
        a.pro_fitness = self.pro_fitness
        return a

    # -----------------------------------------------
    def show(self):

        print(self.chromo)
        print(self.raw_fitness)
        print(self.scl_fitness)
        print(self.pro_fitness)

# -----------------------------------------------
def selectParent():

    assert sel_type in [1, 2, 3], 'Error: Invalid selection type {}'.format(sel_type)

    if sel_type == 1:
        r_wheel = 0
        r = random.random()
        for i in range(pop_size):
            r_wheel += members[i].pro_fitness
            if r < r_wheel:
                return i

    elif sel_type == 2:    
        tour_size = 0.5
        tour_prob = 0.5

        tournament = [random.randint(0, pop_size-1) for i in range(int(tour_size * pop_size))]
        tournament = sorted(tournament, key=lambda i: members[i].profitness, reverse=True)

        for i in range(len(tournament)):
            if random.random() <= tour_prob:
                return i
        return len(tournament) - 1

    elif sel_type == 3:
        return random.randint(0, pop_size-1)

# -----------------------------------------------
def xover(p1, p2):

    assert xover_type in [1, 2, 3], 'Error: Invalid cross over type {}'.format(xover_type)

    c1 = Chromosome()
    c2 = Chromosome()

    if xover_type == 1:
        # -----------------------------------------------
        # Select crossover point
        # -----------------------------------------------
        xp = 1 + (int)(random.random() * (num_genes * gene_size-1))

        # -----------------------------------------------
        # Create child chromosome from parental material
        # -----------------------------------------------
        c1.chromo = p1.chromo[:xp] + p2.chromo[xp:]
        c2.chromo = p2.chromo[:xp] + p1.chromo[xp:]

    elif xover_type == 2:
        # -----------------------------------------------
        # Select crossover points
        # -----------------------------------------------
        xp1 = 1 + (int)(random.random() * (num_genes * gene_size-1))
        xp2 = 1 + (int)(random.random() * (num_genes * gene_size-1))

        if xp1 > xp2:
            xp1, xp2 = xp2, xp1

        # -----------------------------------------------
        # Create child chromosome from parental material
        # -----------------------------------------------
        c1.chromo = p1.chromo[:xp1] + p2.chromo[xp1:xp2] + p1.chromo[xp2:]
        c2.chromo = p2.chromo[:xp1] + p1.chromo[xp1:xp2] + p2.chromo[xp2:]

    return c1, c2

# -----------------------------------------------
def scaleFitness():

    assert scale_type in [0, 1, 2, 3], 'Error: Invalid fitness scaling type'

    sum_sf = 0
    if scale_type == 0:
        for i in range(pop_size):
            members[i].scl_fitness = members[i].raw_fitness + 0.000001
            sum_sf += members[i].scl_fitness

    elif scale_type == 1:
        for i in range(pop_size):
            members[i].scl_fitness = 1/(members[i].raw_fitness + 0.000001)
            sum_sf += members[i].scl_fitness

    elif scale_type == 2:
        member_indices = sorted(range(pop_size), key=lambda i: members[i].raw_fitness, reverse=True)
        for i in range(pop_size):
            members[member_indices[i]].scl_fitness = i
            sum_sf += members[member_indices[i]].scl_fitness

    elif scale_type == 3:
        member_indices = sorted(range(pop_size), key=lambda i: members[i].raw_fitness)
        for i in range(pop_size):
            members[member_indices[i]].scl_fitness = i
            sum_sf += members[member_indices[i]].scl_fitness

    for i in range(pop_size):
        members[i].pro_fitness = members[i].scl_fitness/sum_sf

# -----------------------------------------------
def evolveGeneration(members):

    new_members = []

    for i in range(pop_size):

        p_index1 = selectParent()
        p_index2 = p_index1
        while p_index2 == p_index1:
            p_index2 = selectParent()

        if random.random() <= xover_rate:
            c1, c2 = xover(members[p_index1], members[p_index2])
        else:
            c1, c2 = members[p_index1].clone(), members[p_index2].clone()

        c1.mutate()
        c2.mutate()

        new_members += [c1, c2]

    return new_members

# -----------------------------------------------
param_file = '../onemax.params' #sys.argv[1]
printDec('Parameter filename: {}'.format(param_file))
params = getSettings(param_file)

printDec('Problem name: {}'.format(params['Problem Type']))
exp_id = params['Experiment ID']
random_seed = int(params['Random Number Seed'])
num_runs = int(params['Number of Runs'])
pop_size = int(params['Population Size'])
num_gens = int(params['Generations per Run'])
sel_type = int(params['Selection Method (1)'])
scale_type = int(params['Fitness Scaling Type (2)'])
xover_type = int(params['Crossover Type (3)'])
xover_rate = float(params['Crossover Rate (4)'])
mut_type = int(params['Mutation Type (5)'])
mut_rate = float(params['Mutation Rate (6)'])
num_genes = int(params['Number of Genes/Points (7)'])
gene_size = int(params['Size of Genes (18)'])

out_file = '../results/{}_summary.csv'.format(exp_id)
writeFile(out_file, '')
random.seed(random_seed)
min_or_max  = 'max' if scale_type in [0, 2] else 'min'



best_overall, best_overall_r, best_overall_g = None, -1, -1

stats_overall = []

for r in range(1, num_runs):

    members = [Chromosome() for i in range(pop_size)]
    
    best_of_run, best_of_run_g = None, -1

    stats_all_gen = []

    for g in range(num_gens):

        sum_rf = 0
        sum_rf_2 = 0

        best_of_gen = None

        for i in range(pop_size):

            members[i].calcFitness()

            sum_rf += members[i].raw_fitness
            sum_rf_2 += members[i].raw_fitness ** 2
                
            if best_of_gen == None or \
               (min_or_max == 'max' and best_of_gen.raw_fitness < members[i].raw_fitness) or \
               (min_or_max == 'min' and best_of_gen.raw_fitness > members[i].raw_fitness):

                best_of_gen = members[i].clone()

        scaleFitness()
        members = evolveGeneration(members)

        avg_rf = sum_rf/pop_size
        std_dev_rf = math.sqrt(abs(sum_rf_2-sum_rf**2/pop_size)/(pop_size-1))

        print('{}\t{}\t{}\t{}\t{}'.format(r, g, best_of_gen.raw_fitness, avg_rf, std_dev_rf))

        stats_all_gen.append([r, g, best_of_gen.raw_fitness, avg_rf, std_dev_rf])

        if best_of_gen != None and \
           (best_of_run == None or \
           (min_or_max == 'max' and best_of_run.raw_fitness < best_of_gen.raw_fitness) or \
           (min_or_max == 'min' and best_of_run.raw_fitness > best_of_gen.raw_fitness)):
            
            best_of_run, best_of_run_g = best_of_gen.clone(), g

    stats_all_gen.append([])
    stats_all_gen.append([])
    
    writeDataTable(stats_all_gen, out_file, mode='a')

    printDec('Run: {}, Best gen: {}, Best fitness: {}'.format(r, best_of_run_g, best_of_run.raw_fitness))
    stats_overall.append([r, best_of_run_g, best_of_run.raw_fitness])

    if best_of_run != None and \
       (best_overall == None or \
       (min_or_max == 'max' and best_overall.raw_fitness < best_of_run.raw_fitness) or \
       (min_or_max == 'min' and best_overall.raw_fitness > best_of_run.raw_fitness)):

        best_overall, best_overall_g, best_overall_r = best_of_run.clone(), g, r

writeDataTable(stats_overall, out_file, mode='a')









    
