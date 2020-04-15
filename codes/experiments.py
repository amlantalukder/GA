import search as s

params = {\
 'Experiment ID'                :'onemax', \
 'Problem Type'                 :'OM', \
 'Data Input File Name'         :'NA', \
 'Number of Runs'               :3, \
 'Generations per Run'          :1000, \
 'Population Size'              :300, \
 'Selection Method (1)'         :4, \
 'Fitness Scaling Type (2)'     :2, \
 'Crossover Type (3)'           :1, \
 'Crossover Rate (4)'           :0.8, \
 'Mutation Type (5)'            :1, \
 'Mutation Rate (6)'            :0.001, \
 'Random Number Seed'           :75982, \
 'Number of Genes/Points (7)'   :1, \
 'Size of Genes (18)'           :200, \
 'Decay'                        :0.8, \
 'Queue Length'                 :100, \
 'Depth'                        :5
}

gens_per_run = [100]
pop_sizes = [100]
qlens = [10, 100]
decay_values = [0.8]
depths = [2, 5, 10]

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
     params['Experiment ID'] = 'onemax_g-{}_p-{}_q-{}_d-{}_dp-{}'.format(g, p, qlen, decay, depth)
     s.params = s.Parameters(params)
     s.runGA(s.params)