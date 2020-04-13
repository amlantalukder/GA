from utils import *
import pdb

class Parameters:

    def __init__(self, params):

        self.problem_type = self.getDictVal(params, 'Problem Type')
        self.exp_id = self.getDictVal(params, 'Experiment ID')
        self.random_seed = int(self.getDictVal(params, 'Random Number Seed'))
        self.num_runs = int(self.getDictVal(params, 'Number of Runs'))
        self.pop_size = int(self.getDictVal(params, 'Population Size'))
        self.num_gens = int(self.getDictVal(params, 'Generations per Run'))
        self.sel_type = int(self.getDictVal(params, 'Selection Method (1)'))
        self.scale_type = int(self.getDictVal(params, 'Fitness Scaling Type (2)'))
        self.xover_type = int(self.getDictVal(params, 'Crossover Type (3)'))
        self.xover_rate = float(self.getDictVal(params, 'Crossover Rate (4)'))
        self.mut_type = int(self.getDictVal(params, 'Mutation Type (5)'))
        self.mut_rate = float(self.getDictVal(params, 'Mutation Rate (6)'))
        self.num_genes = int(self.getDictVal(params, 'Number of Genes/Points (7)'))
        self.gene_size = int(self.getDictVal(params, 'Size of Genes (18)'))
        self.decay = float(self.getDictVal(params, 'Decay'))
        self.qlen = int(self.getDictVal(params, 'Queue Length'))

        # Rank selection
        if self.sel_type == 4:
            self.sel_type = 1
            if self.scale_type == 0:
                self.scale_type = 2
            elif self.scale_type == 1:
                self.scale_type = 3

    def getDictVal(self, d, key):
        return d[key] if key in d else None


# -----------------------------------------------
class Chromosome:

    # -----------------------------------------------
    def __init__(self):

        self.raw_fitness = -1
        self.scl_fitness = -1
        self.pro_fitness = -1

        self.optype = None
        self.p1 = None
        self.p2 = None
        self.credit_xover = None
        self.credit_mut = None
        
        self.chromo = ''
        for i in range(params.num_genes):
            for j in range(params.gene_size):
                self.chromo += '0' if random.random() < 0.5 else '1'

    # -----------------------------------------------
    def mutate(self, p):

        assert params.mut_type == 1, 'Error: Invalid mutation type {}'.format(params.mut_type)

        x = ''
        for i in range(params.num_genes * params.gene_size):
            if random.random() < params.mut_rate:
                x += '0' if self.chromo[i] == '1' else '0'
            else:
                x += self.chromo[i]
        self.chromo = x

        self.setOperatorHierarchy('mut', p1=p)

    # -----------------------------------------------
    def calcFitness(self):

        self.raw_fitness = sum([int(i) for i in self.chromo])

    # -----------------------------------------------
    def setOperatorHierarchy(self, optype, p1=None, p2=None):
        self.optype = optype
        self.p1 = p1
        self.p2 = p2
        self.setOpertorCredit()

    # -----------------------------------------------
    def setOpertorCredit(self):
        credit_xover, credit_mut = 0, 0

        level = 0
        count_xover, count_mut = (1, 0) if self.optype == 'xover' else (0, 1)
        credit_xover += count_xover * params.decay ** level
        credit_mut += count_mut * params.decay ** level

        parents = self.p1, self.p2
        while len(parents) > 0:
            count_xover, count_mut = 0, 0
            newparents = []
            for p in parents:
                if not p: continue
                if p.optype:
                    if p.optype == 'xover':
                        count_xover += 1
                    else:
                        count_mut += 1
                if p.p1: newparents.append(p.p1)
                if p.p2: newparents.append(p.p2)
            parents = newparents
            level += 1
            credit_xover += count_xover * params.decay ** level
            credit_mut += count_mut * params.decay ** level

        self.credit_xover = credit_xover
        self.credit_mut = credit_mut

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

    assert params.sel_type in [1, 2, 3], 'Error: Invalid selection type {}'.format(params.sel_type)

    if params.sel_type == 1:
        r_wheel = 0
        r = random.random()
        for i in range(params.pop_size):
            r_wheel += members[i].pro_fitness
            if r < r_wheel:
                return i
        return -1

    elif params.sel_type == 2:
        tour_size = 0.5
        tour_prob = 0.5

        tournament = random.sample(range(params.pop_size), int(tour_size*params.pop_size))
        tournament = sorted(tournament, key=lambda i: members[i].pro_fitness, reverse=True)
        
        for i in range(len(tournament)):
            if random.random() < tour_prob:
                return tournament[i]
        return tournament[-1]

    elif params.sel_type == 3:
        return random.randint(0, params.pop_size-1)

# -----------------------------------------------
def xover(p1, p2):

    assert params.xover_type in [1, 2, 3], 'Error: Invalid cross over type {}'.format(params.xover_type)

    c1 = Chromosome()
    c2 = Chromosome()

    if params.xover_type == 1:
        # -----------------------------------------------
        # Select crossover point
        # -----------------------------------------------
        xp = random.randint(0, params.num_genes * params.gene_size-1)

        # -----------------------------------------------
        # Create child chromosome from parental material
        # -----------------------------------------------
        c1.chromo = p1.chromo[:xp] + p2.chromo[xp:]
        c2.chromo = p2.chromo[:xp] + p1.chromo[xp:]

        c1.setOperatorHierarchy('xover', p1, p2)
        c2.setOperatorHierarchy('xover', p2, p1)

    elif xover_type == 2:
        # -----------------------------------------------
        # Select crossover points
        # -----------------------------------------------
        xp1 = random.randint(0, params.num_genes * params.gene_size-1)
        xp2 = random.randint(0, params.num_genes * params.gene_size-1)

        if xp1 > xp2:
            xp1, xp2 = xp2, xp1

        # -----------------------------------------------
        # Create child chromosome from parental material
        # -----------------------------------------------
        c1.chromo = p1.chromo[:xp1] + p2.chromo[xp1:xp2] + p1.chromo[xp2:]
        c2.chromo = p2.chromo[:xp1] + p1.chromo[xp1:xp2] + p2.chromo[xp2:]

        c1.setOperatorHierarchy('xover', p1, p2)
        c2.setOperatorHierarchy('xover', p2, p1)

    return c1, c2

# -----------------------------------------------
def scaleFitness():

    assert params.scale_type in [0, 1, 2, 3], 'Error: Invalid fitness scaling type {}'.format(params.scale_type)

    sum_sf = 0
    if params.scale_type == 0:
        for i in range(params.pop_size):
            members[i].scl_fitness = members[i].raw_fitness + 0.000001
            sum_sf += members[i].scl_fitness

    elif params.scale_type == 1:
        for i in range(params.pop_size):
            members[i].scl_fitness = 1/(members[i].raw_fitness + 0.000001)
            sum_sf += members[i].scl_fitness

    elif params.scale_type == 2:
        member_indices = sorted(range(params.pop_size), key=lambda i: members[i].raw_fitness)
        for i in range(params.pop_size):
            members[member_indices[i]].scl_fitness = i
            sum_sf += members[member_indices[i]].scl_fitness

    elif params.scale_type == 3:
        member_indices = sorted(range(params.pop_size), key=lambda i: members[i].raw_fitness, reverse=True)
        for i in range(params.pop_size):
            members[member_indices[i]].scl_fitness = i
            sum_sf += members[member_indices[i]].scl_fitness
    
    for i in range(params.pop_size):
        members[i].pro_fitness = members[i].scl_fitness/sum_sf

# -----------------------------------------------
def evolveGeneration(members):

    worst_two_members_info = []

    for i in range(params.pop_size):
        # When the list is empty
        if len(worst_two_members_info) == 0:
            worst_two_members_info.append((members[i].pro_fitness, i))
        # When the list has only one entry
        elif len(worst_two_members_info) == 1:
            # The new fitness is worse than the entry, enlist it
            if worst_two_members_info[0][0] > members[i].pro_fitness:
                worst_two_members_info.append((members[i].pro_fitness, i))
            # The new fitness is better than the entry,
            # move the entry as the worst fitness,
            # enlist the new one as the second worst
            else:
                worst_two_members_info.append(worst_two_members_info[0])
                worst_two_members_info[0] = (members[i].pro_fitness, i)
        # When the list has both entries,
        # if the new fitness is worst than the enlisted worse,
        # then move the enlisted worst to second worst and
        # enlist the new fitness as the worst
        elif worst_two_members_info[1][0] > members[i].pro_fitness:
            worst_two_members_info[0] = worst_two_members_info[1]
            worst_two_members_info[1] = (members[i].pro_fitness, i)
        # When the list has both entries,
        # if the new fitness is better than the enlisted worse but worse than the second worst
        # then enlist the new fitness as the second worst
        elif worst_two_members_info[0][0] > members[i].pro_fitness:
            worst_two_members_info[0] = (members[i].pro_fitness, i)

    p_index1 = selectParent()
    p_index2 = p_index1
    while p_index2 == p_index1:
        p_index2 = selectParent()

    p1, p2 = members[p_index1], members[p_index2]

    global queue, total_credit_xover, total_credit_mut, num_xover, num_mut

    #print(total_credit_xover, total_credit_mut, num_xover, num_mut)

    if queue.size() >= params.qlen:
        xover_ratio = (total_credit_xover/num_xover) if num_xover > 0 else 0
        mut_ratio = (total_credit_mut/num_mut) if num_mut > 0 else 0
        xover_rate = xover_ratio/(xover_ratio + mut_ratio)
    else:
        xover_rate = params.xover_rate

    if random.random() < xover_rate:
        c1, c2 = xover(p1, p2)
        num_xover += 2
    else:
        c1, c2 = p1.clone(), p2.clone()
        c1.mutate(p1)
        c2.mutate(p2)
        num_mut += 2

    for c in [c1, c2]:

        if queue.size() > params.qlen:
            optype, credit_xover, credit_mut = queue.dequeue()
            total_credit_xover -= credit_xover
            total_credit_mut -= credit_mut
            if optype == 'xover':
                num_xover -= 1
            elif optype == 'mut':
                num_mut -= 1

        total_credit_xover += c.credit_xover
        total_credit_mut += c.credit_mut

        queue.enqueue((c.optype, c.credit_xover, c.credit_mut))

    members[worst_two_members_info[0][1]] = c1
    members[worst_two_members_info[1][1]] = c2

    return members

# -----------------------------------------------
# Load parameters
# -----------------------------------------------
try:
    param_file = sys.argv[1]
except:
    print("Error: invalid parameter file")
    sys.exit(1)

verbose = True if len(sys.argv) > 2 and sys.argv[2] == '1' else False
    
printDec('Parameter filename: {}'.format(param_file))
params = Parameters(getSettings(param_file))

printDec('Problem name: {}'.format(params.problem_type))

out_file = 'results/{}_summary.csv'.format(params.exp_id)
writeFile(out_file, '')
random.seed(params.random_seed)
min_or_max  = 'max' if params.scale_type in [0, 2] else 'min'

# -----------------------------------------------
# Run GA
# -----------------------------------------------
best_overall, best_overall_r, best_overall_g = None, -1, -1

stats_overall = []

for r in range(1, params.num_runs+1):

    members = [Chromosome() for i in range(params.pop_size)]
    
    best_of_run, best_of_run_g = None, -1

    num_xover, num_mut = 0, 0
    total_credit_xover, total_credit_mut = 0, 0
    queue = Queue()

    stats_all_gen = []

    perc = 10

    for g in range(1, params.num_gens+1):

        if not verbose: perc = showPercBar(g, params.num_gens, perc)

        sum_rf = 0
        sum_rf_2 = 0

        best_of_gen = None

        for i in range(params.pop_size):

            members[i].calcFitness()

            sum_rf += members[i].raw_fitness
            sum_rf_2 += members[i].raw_fitness ** 2
                
            if best_of_gen == None or \
               (min_or_max == 'max' and best_of_gen.raw_fitness < members[i].raw_fitness) or \
               (min_or_max == 'min' and best_of_gen.raw_fitness > members[i].raw_fitness):

                best_of_gen = members[i].clone()

        scaleFitness()
        members = evolveGeneration(members)

        avg_rf = sum_rf/params.pop_size
        std_dev_rf = math.sqrt(abs(sum_rf_2-sum_rf**2/params.pop_size)/(params.pop_size-1))

        if verbose: print('{}\t{}\t{}\t{}\t{}'.format(r, g, best_of_gen.raw_fitness, avg_rf, std_dev_rf))

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

        best_overall, best_overall_g, best_overall_r = best_of_run.clone(), best_of_run_g, r

writeDataTable(stats_overall, out_file, mode='a')
printDec('Best Run: {}, Best gen: {}, Best fitness: {}'.format(best_overall_r, best_overall_g, best_overall.raw_fitness))








    
