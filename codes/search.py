from utils import *
import pdb
import math

# -----------------------------------------------
class Parameters:

    # -----------------------------------------------
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
        self.depth = int(self.getDictVal(params, 'Depth'))

        # Rank selection
        if self.sel_type == 4:
            self.sel_type = 1
            if self.scale_type == 0:
                self.scale_type = 2
            elif self.scale_type == 1:
                self.scale_type = 3

    # -----------------------------------------------
    def getDictVal(self, d, key):
        return d[key] if key in d else None


# -----------------------------------------------
class Chromosome:

    # -----------------------------------------------
    def __init__(self, chromo=None):

        self.raw_fitness = None
        self.scl_fitness = None
        self.pro_fitness = None

        self.optype = None
        self.p1 = None
        self.p2 = None
        self.credit = {}
        self.chromo = ''

        if chromo == None:
            for i in range(params.num_genes):
                for j in range(params.gene_size):
                    self.chromo += '0' if random.random() < 0.5 else '1'
        else:
            self.chromo = chromo

    # -----------------------------------------------
    def mutate(self):

        assert params.mut_type in [1, 2], 'Error: Invalid mutation type {}'.format(params.mut_type)

        # Flip
        if params.mut_type == 1:
            x = ''
            for i in range(len(self.chromo)):
                if random.random() < params.mut_rate:
                    x += '0' if self.chromo[i] == '1' else '0'
                else:
                    x += self.chromo[i]
            self.chromo = x

        # Swap
        elif params.mut_type == 2:

            x = [i for i in self.chromo]
            random.shuffle(x)
            self.chromo = ''.join(x)

    # -----------------------------------------------
    def calcFitness(self):
        if params.problem_type == 'OM':
            self.raw_fitness = sum([int(i) for i in self.chromo])
        if params.problem_type == 'BF6':
            self.raw_fitness = Chromosome.calcBinaryF6(self.chromo)

    # -----------------------------------------------
    def setOperatorHierarchy(self, optype, p1=None, p2=None):
        self.optype = optype
        self.p1 = p1
        self.p2 = p2
        self.setOpertorCredit()

    # -----------------------------------------------
    def setOpertorCredit(self):

        self.credit = dict([[optype, 0] for optype in operator_credit_info])

        level = 1
        nodes = [self]
        while len(nodes) > 0 and level <= params.depth:

            new_nodes = []

            for p in nodes:
                if not p: continue
                if p.optype: self.credit[p.optype] += params.decay ** (level - 1)
                if p.p1: new_nodes.append(p.p1)
                if p.p2: new_nodes.append(p.p2)

            nodes = new_nodes
            level += 1

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

    @staticmethod
    def calcBinaryF6(bitstr):
        l = len(bitstr) // 2
        ratioNum = 200 / 4194303
        x = int(bitstr[:l], 2) * ratioNum - 100
        y = int(bitstr[l:], 2) * ratioNum - 100
        sq = x ** 2 + y ** 2
        return (0.5 - (((math.sin(math.sqrt(sq)) ** 2) - 0.5) / (1.0 + 0.001 * sq ** 2)))


# -----------------------------------------------
def selectParent(members):
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
        tour_size = 0.1
        tour_prob = 1

        tournament = random.sample(range(params.pop_size), int(tour_size * params.pop_size))
        tournament = sorted(tournament, key=lambda i: members[i].pro_fitness, reverse=True)

        for i in range(len(tournament)):
            if random.random() < tour_prob:
                return tournament[i]
        return tournament[-1]

    elif params.sel_type == 3:
        return random.randint(0, params.pop_size - 1)


# -----------------------------------------------
def xover(p1, p2):
    assert params.xover_type in [1, 2, 3], 'Error: Invalid cross over type {}'.format(params.xover_type)

    c1 = Chromosome()
    c2 = Chromosome()

    # One point
    if params.xover_type == 1:
        # -----------------------------------------------
        # Select crossover point
        # -----------------------------------------------
        xp = random.randint(0, params.num_genes * params.gene_size - 1)

        # -----------------------------------------------
        # Create child chromosome from parental material
        # -----------------------------------------------
        c1.chromo = p1.chromo[:xp] + p2.chromo[xp:]
        c2.chromo = p2.chromo[:xp] + p1.chromo[xp:]

    # Two point
    elif params.xover_type == 2:
        # -----------------------------------------------
        # Select crossover points
        # -----------------------------------------------
        xp1 = random.randint(0, params.num_genes * params.gene_size - 1)
        xp2 = random.randint(0, params.num_genes * params.gene_size - 1)

        if xp1 > xp2:
            xp1, xp2 = xp2, xp1

        # -----------------------------------------------
        # Create child chromosome from parental material
        # -----------------------------------------------
        c1.chromo = p1.chromo[:xp1] + p2.chromo[xp1:xp2] + p1.chromo[xp2:]
        c2.chromo = p2.chromo[:xp1] + p1.chromo[xp1:xp2] + p2.chromo[xp2:]

    # Uniform
    elif params.xover_type == 2:

        x1, x2 = '', ''

        for i in range(len(p1.chromo)):

            # -----------------------------------------------
            # Create child chromosome from parental material
            # -----------------------------------------------
            if random.random() < 0.5:
                x1 += p1.chromo[i]
                x2 += p2.chromo[i]
            else:
                x1 += p2.chromo[i]
                x2 += p1.chromo[i]

        c1 = Chromosome(x1)
        c2 = Chromosome(x2)

    return c1, c2


# -----------------------------------------------
def scaleFitness(members):
    assert params.scale_type in [0, 1, 2, 3], 'Error: Invalid fitness scaling type {}'.format(params.scale_type)

    sum_sf = 0
    if params.scale_type == 0:
        for i in range(params.pop_size):
            members[i].scl_fitness = members[i].raw_fitness + 0.000001
            sum_sf += members[i].scl_fitness

    elif params.scale_type == 1:
        for i in range(params.pop_size):
            members[i].scl_fitness = 1 / (members[i].raw_fitness + 0.000001)
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
        members[i].pro_fitness = members[i].scl_fitness / sum_sf

    return members


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

    p_index1 = selectParent(members)
    p_index2 = selectParent(members)
    while p_index2 == p_index1:
        p_index2 = selectParent(members)

    p1, p2 = members[p_index1], members[p_index2]

    global queue, operator_credit_info

    # print(queue.size(), operator_credit_info)

    s = 0
    if queue.isFull():
        for optype in operator_credit_info:
            total_credit, num = operator_credit_info[optype][:2]
            operator_credit_info[optype][2] = (total_credit / num) if num > 0 else 0
            s += operator_credit_info[optype][2]

    # -----------------------------------------------
    # Adaptive probabilities for the operators
    # -----------------------------------------------
    if s > 0:
        for optype in operator_credit_info:
            operator_credit_info[optype][2] /= s
    # -----------------------------------------------
    # Default probabilities for the operators
    # -----------------------------------------------
    else:
        for optype in operator_credit_info:
            operator_credit_info[optype][2] = 1 / len(operator_credit_info)

    toss, r = random.random(), 0
    for optype in operator_credit_info:
        if toss < r + operator_credit_info[optype][2]:
            if optype[:3] == 'xov':
                params.xover_type = int(optype.split('_')[1])
                c1, c2 = xover(p1, p2)
                c1.setOperatorHierarchy(optype, p1, p2)
                c2.setOperatorHierarchy(optype, p2, p1)
                operator_credit_info[optype][0] += (c1.credit[optype] + c2.credit[optype])
                operator_credit_info[optype][1] += 2
                break
            elif optype[:3] == 'mut':
                params.mut_type = int(optype.split('_')[1])
                c1, c2 = p1.clone(), p2.clone()
                c1.mutate()
                c1.setOperatorHierarchy(optype, p1)
                c2.mutate()
                c2.setOperatorHierarchy(optype, p2)
                operator_credit_info[optype][0] += (c1.credit[optype] + c2.credit[optype])
                operator_credit_info[optype][1] += 2
                break
        r += operator_credit_info[optype][2]

    for c in [c1, c2]:

        if queue.isFull():
            optype, credit = queue.dequeue()
            operator_credit_info[optype][0] -= credit[optype]
            operator_credit_info[optype][1] -= 1
            if operator_credit_info[optype][1] < 0:
                pdb.set_trace()

        queue.enqueue((c.optype, c.credit))

    members[worst_two_members_info[0][1]] = c1
    members[worst_two_members_info[1][1]] = c2

    return members


# -----------------------------------------------
def runGA(params, verbose=False):
    global queue, operator_credit_info

    printDec('Problem name: {}'.format(params.problem_type))

    out_file = 'results/{}_summary.csv'.format(params.exp_id)
    writeFile(out_file, '')
    ops_probs_out_file = 'results/{}_operator_probabilities.csv'.format(params.exp_id)
    random.seed(params.random_seed)
    min_or_max = 'max' if params.scale_type in [0, 2] else 'min'

    # -----------------------------------------------
    # Run GA
    # -----------------------------------------------
    best_overall, best_overall_r, best_overall_g = None, -1, -1

    stats_overall = []
    ops_credit_probs = [[] for i in range(params.num_gens + 1)]

    for r in range(1, params.num_runs + 1):

        ops_credit_probs[0] += [optype for optype in operator_credit_info]

        members = [Chromosome() for i in range(params.pop_size)]

        best_of_run, best_of_run_g = None, -1

        operator_credit_info = {'xover_1': [0, 0, 0], 'mut_1': [0, 0, 0]}
        queue = Queue(params.qlen)

        stats_all_gen = []

        perc = 10

        for g in range(1, params.num_gens + 1):

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

            members = scaleFitness(members)
            members = evolveGeneration(members)

            ops_credit_probs[g] += [operator_credit_info[optype][2] for optype in operator_credit_info]

            avg_rf = sum_rf / params.pop_size
            std_dev_rf = math.sqrt(abs(sum_rf_2 - sum_rf ** 2 / params.pop_size) / (params.pop_size - 1))

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
    writeDataTable(ops_credit_probs, ops_probs_out_file)
    printDec(
        'Best Run: {}, Best gen: {}, Best fitness: {}'.format(best_overall_r, best_overall_g, best_overall.raw_fitness))


operator_credit_info = {'xover_1': [0, 0, 0], 'mut_1': [0, 0, 0]}
queue = None

if __name__ == '__main__':

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

    runGA(params, verbose)