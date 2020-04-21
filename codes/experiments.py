import search as s

params = {\
 'Experiment ID'                :'BF6', \
 'Problem Type'                 :'BF6', \
 'Data Input File Name'         :'NA', \
 'Number of Runs'               :30, \
 'Generations per Run'          :1000, \
 'Population Size'              :300, \
 'Selection Method (1)'         :2, \
 'Fitness Scaling Type (2)'     :0, \
 'Crossover Type (3)'           :1, \
 'Crossover Rate (4)'           :0.8, \
 'Mutation Type (5)'            :1, \
 'Mutation Rate (6)'            :0.05, \
 'Random Number Seed'           :75982, \
 'Number of Genes/Points (7)'   :2, \
 'Size of Genes (18)'           :22, \
 'Decay'                        :0.8, \
 'Queue Length'                 :100, \
 'Depth'                        :5
}

gens_per_run = [5000]
pop_sizes = [100] #[300]
qlens = [10] #[10, 500]
decay_values = [0.5] #[0.5, 0.95]
depths = [1, 5] # [1, 8]

for g in gens_per_run:
 params['Generations per Run'] = g
 for p in pop_sizes:
  params['Population Size'] = p
  for qlen in qlens:
   params['Queue Length'] = qlen
   for decay in decay_values:
    params['Decay'] = decay
    for depth in depths:
     params['Depth'] = depth
     params['Experiment ID'] = 'BF6_g-{}_p-{}_q-{}_d-{}_dp-{}'.format(g, p, qlen, decay, depth)
     s.params = s.Parameters(params)
     s.runGA(s.params)