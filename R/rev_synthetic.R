rev_synth_data = function(cohort, n, tree = NULL, replace = FALSE,
                          onlyWithMultipleSolutions = TRUE,
                          n.truncalMuts = 1,
                          n.subclonalMuts = 30,
                          n.childrenPerMut = 3,
                          n.driverPerCloneMuts = 3)
{
  print(cohort)
  
  if(onlyWithMultipleSolutions){
    cat(cyan('All patients in the cohort:'), length(cohort$patients), '\n')
    which.patients = sapply(
      cohort$patients,
      function(w) {
        # keys = lapply(
        #   cohort$phylogenies[[w]],
        #   function(p) paste(sort(DataFrameToEdges(p$transfer$driver)), collapse = ' '))
        # keys = Reduce(rbind, keys)
        # print(duplicated(keys))
        # 
        
        length(cohort$phylogenies[[w]]) > 1
        })
    
    cohort$patients = cohort$patients[which.patients]
    cat(cyan('Patients that can transfer different orderings:'), length(cohort$patients), '\n')
  }
  
  patients = sample(cohort$patients, n, replace)
  dataset = cohort$dataset[cohort$patients %in% patients, , drop = FALSE]

  # Sample randomly a subset of the patients from the template cohort -- and decide their true models
  clonal.trees = cohort$phylogenies[patients]
  true.trees.idx = lapply(clonal.trees, 
                      function(w) sample(1:length(w), 1))
  true.trees = lapply(
    sapply(1:length(true.trees.idx), list), 
    function(w) 
      clonal.trees[[names(true.trees.idx)[w]]][[true.trees.idx[[w]]]] 
    )
  names(true.trees) = names(clonal.trees)
  
  #  Pick a tree for the evolutionary relations among synthetic drivers
  if(is.null(tree)) {
    subcl.tree = make_tree(n.subclonalMuts, children = n.childrenPerMut) 
    subcl.tree = set.vertex.attribute(subcl.tree, 'name', value = paste('S', 1:n.subclonalMuts, sep =''))
    
    cl.trunk = make_tree(n.truncalMuts, children = 1)
    cl.trunk = set.vertex.attribute(cl.trunk, 'name', value = paste('C', 1:n.truncalMuts, sep =''))
    cl.trunk = set.vertex.attribute(cl.trunk, 'name', value = paste('S', 1, sep =''), index = n.truncalMuts)
    
    tree = cl.trunk + subcl.tree
    
    random_code = paste(LETTERS[sample(1:15, 3, replace = T)], collapse = '')
    tree = set.vertex.attribute(tree, "name", value = paste(random_code, 1:vcount(tree), sep =''))
  }
  
  pdf('syth.pdf')
  plot(tree, layout = layout.reingold.tilford(tree))
  dev.off()
  
  morph = function(clones.tree, trajectories) {

    # All paths in the trajectory "space"    
    asp = all_simple_paths(trajectories, from = V(trajectories)$name[1], to = leaves(as_adj(tree, sparse = FALSE)))
    asp = asp[[sample(1:length(asp), 1)]]
    path = names(asp)
    
    
    cat('Morphing: ', paste(path, collapse = ' -> '), '\n')
    adj_mat = clones.tree$adj_mat
    # cat('\t@ ', colnames(adj_mat), '\n')
    
    mapping = NULL
    idx = children(adj_mat, 'GL')
    
    repeat{
      # We map a random number of entries from the path, in linear order
      mapped = min(length(path), sample(1:n.driverPerCloneMuts, 1))
      mapping = append(mapping, list(path[1:mapped]))
      names(mapping)[length(mapping)] = idx
      
      cat('\t Mapping', path[1:mapped], '@', idx, '-- ', mapped, ' hops mapped\n')
      
      # we stop if we do not have more to sample from the path, or if we do not have more clones
      
      ch = children(adj_mat, idx)
      
      # print(ch)
      # print(mapped)
      # print(length(path))
      # 
      if(mapped == length(path)) break;
      if(length(ch) == 0) break;
      
      idx = sample(ch, 1)
      path = path[(mapped + 1):length(path)]
    }
    
    return(mapping)
  }
  
  new.dataset = NULL
  new.patients = NULL
  for(p in 1:length(patients))
  {
    cat(cyan('\nMapping drivers to clone tree for', patients[p]), '\n')
    mapping = morph(true.trees[[p]], tree)
    # dataset[dataset$patientID == patients[p], ]
    
    # patient.dataset = true.trees[[p]]$dataset
    
    patient.dataset = cohort$dataset[cohort$dataset$patientID == patients[p], ]
    
    data  = lapply(
      sapply(1:length(mapping), list), 
      function(w){
        w = names(mapping)[w]  
        
        template = patient.dataset[patient.dataset$cluster == as.character(w), ]
        nentries = length(mapping[[w]])
        
        # cat(red('template'), w)
        # print(template)
        # print(dataset[dataset$patientID == patients[p], ])
        
        new.entries = NULL
        for(t in 1:nentries) new.entries = rbind(new.entries, template[1, ])
        new.entries[,  'variantID'] = mapping[[w]]
        new.entries[, 'is.driver'] = TRUE
        
        # print(new.entries)
        return(new.entries)
    })
    data = Reduce(rbind, data)
    # data.thispatient = dataset[dataset$patientID == patients[p] & !dataset$is.driver, ]
    
    data.thispatient = 
      rbind(
        patient.dataset[!patient.dataset$is.driver, ],
        # cohort$dataset[cohort$dataset$patientID == patients[p] & !cohort$dataset$is.driver, ],
        data)
    data.thispatient$patientID = paste('SynthPatient', p, sep = '-')
    
    new.dataset = rbind(new.dataset, data.thispatient)
  }  

  # we copy a dataset now...
  cat(cyan('* Refactoring a new cohort\n'))

  y = revolver_cohort(
    dataset = new.dataset,
    CCF.parser = cohort$CCF.parser,
    annotation = 'Syntetic cohort cohort of 100 patients',
    options = list(
      ONLY.DRIVER = FALSE, 
      MIN.CLUSTER.SIZE = 0)
  )
  
  names(true.trees.idx) = y$patients
  
  cat(cyan('* Refactoring phylogienes\n'))
  
  phylos = cohort$phylogenies[patients]
  names(phylos) = y$patients
  # 
  # for(p in names(phylos))
  # {
  #   cat('\n', p, '@ phylogeny num:')
  #   
  #   for(ph in 1:length(phylos[[p]]))
  #   {
  #     cat(ph , ' ')
  #     phylos[[p]][[ph]]$dataset = y$dataset[y$dataset$patientID == p, ]
  #     phylos[[p]][[ph]]$patient_ID = p 
  #     phylos[[p]][[ph]]$CCF = clusters.table(
  #       phylos[[p]][[ph]]$dataset,
  #       phylos[[p]][[ph]]$samples_ID)
  #     
  #     phylos[[p]][[ph]]$transfer = information.transfer(phylos[[p]][[ph]])
  #   }
  # }
  # 
  
  for(p in names(phylos))
  {
    matrices = lapply(phylos[[p]], function(w) {
      M = w$adj_mat
      M[rownames(M) != 'GL', colnames(M) != 'GL', drop = FALSE]})
    
    scores = unlist(lapply(phylos[[p]], function(w) w$score))
    
    y = compute_rev_phylogenies(
      y, 
      p, 
      precomputed.trees = matrices,
      precomputed.scores = scores,
      options = c(overwrite = TRUE, store.max = 200)
    )
    
  }
  

  cat(cyan('TL identifiability:\n'))
  w = sapply(
      y$patients,
      function(w) {
        keys = lapply(
           y$phylogenies[[w]],
           function(p) paste(sort(DataFrameToEdges(p$transfer$driver)), collapse = ' '))
         keys = Reduce(rbind, keys)
         # print(keys)
         cat(w, ' --> # Phylo', length(y$phylogenies[[w]]), ' # Transf', length(unique(keys)), '\n')
      })
    
  

  return(list(cohort = y, true.models = unlist(true.trees.idx), true.orderings = tree))
}


rev_analyze_experiment = function(fit, dataset)
{
  performance.fit = data.frame(
    fitID = fit$fit$solutionID,
    trueID = dataset$true.models)
  
  cohort = dataset$cohort
  
  w = lapply(
    cohort$patients,
    function(w) {
      keys = lapply(
        cohort$phylogenies[[w]],
        function(p) paste(sort(DataFrameToEdges(p$transfer$driver)), collapse = ' '))
      keys = Reduce(rbind, keys)
      return(
        data.frame(
          NumPhylo = length(cohort$phylogenies[[w]]), 
          NumTR = length(unique(keys)))
        )
    })
  performance.fit = cbind(performance.fit, Reduce(rbind, w))
  
  w =  lapply(1:nrow(performance.fit),
          function(w){
            
            fT = fit$fit$phylogenies[[rownames(performance.fit)[w]]]$transfer$driver
            tT = cohort$phylogenies[[rownames(performance.fit)[w]]][[performance.fit[w, 'trueID']]]$transfer$driver
            
            # print(fT)
            # print(tT)
            # print(performance.fit[w, 'trueID'])
            # print(rownames(performance.fit)[w])
            # print(cohort$phylogenies[[rownames(performance.fit)[w]]])
            
            found = length(intersect(DataFrameToEdges(fT), DataFrameToEdges(tT)))
            
            return(data.frame(
              TLdrvTP =  found/nrow(fT),
              TLdrvFP = 1 - found/nrow(fT)
            ))
    })
    
    performance.fit = cbind(performance.fit, Reduce(rbind, w))
    
    w =  lapply(1:nrow(performance.fit),
                function(w){
                  fT = fit$phylogenies[[rownames(performance.fit)[w]]][[1]]$transfer$driver
                  tT = cohort$phylogenies[[rownames(performance.fit)[w]]][[performance.fit[w, 'trueID']]]$transfer$driver
                  
                  found = length(intersect(DataFrameToEdges(fT), DataFrameToEdges(tT)))
                  
                  return(data.frame(
                    NTLdrvTP =  found/nrow(fT),
                    NTLdrvFP = 1 - found/nrow(fT)
                  ))
                })
 
    performance.fit = cbind(performance.fit, Reduce(rbind, w))
    
    list(performance.fit,
         overall = list(
           TLdrvTP = mean(performance.fit$TLdrvTP),
           TLdrvFP = mean(performance.fit$TLdrvFP),
           NTLdrvTP = mean(performance.fit$NTLdrvTP),
           NTLdrvFP = mean(performance.fit$NTLdrvFP),
           NumPhylo = mean(performance.fit$NumPhylo),
           NumTR = mean(performance.fit$NumTR)
           )
    )
           
}


#####################################################################
##################################################################### DENOVO
#####################################################################

# rev_generate_true_model = function(n, r, noise = c(0, 0.02, 0.05, 0.1), identify = 0)
# {
# 
#   hops = NULL
#   # for(i in 1:n){
#   #   ch = paste(i, 1:sample(1:2, 1), sep = '-') 
#   #   hops[[i]] = expand.grid(from = i, to = ch, stringsAsFactors = FALSE)
#   # }
#   
#   # for(i in 2:n){
#   #   hops[[i-1]] = rbind(
#   #     hops[[i-1]],
#   #     data.frame(from = sample(hops[[i-1]]$to, 1), to = i, stringsAsFactors = FALSE)
#   #   )
#   # }
#   # tree = Reduce(rbind, hops)
#   
#   hops = expand.grid(from = 'GL', to = 1, stringsAsFactors = FALSE)
#   
#   for(i in 2:n){
#     hops = rbind(hops,
#       data.frame(from = sample(hops$to, 1), to = i, stringsAsFactors = FALSE)
#     )
#   }
#   plot(graph_from_adjacency_matrix(DataFrameToMatrix(hops)))
#   hops = hops[hops$from != 'GL', , drop = FALSE]
#   
#   tree = hops
#   clones = unique(unlist(tree))
# 
#   data = matrix(0, nrow = length(clones), ncol = r)
#   colnames(data) = paste('R', 1:r, sep = '')
#   rownames(data) = clones
#   
#   # G = graph_from_adjacency_matrix(DataFrameToMatrix(tree))
#   # plot(G, layout = layout.reingold.tilford(G, root = '1'))
#   # 
#   Q = '1'
#   data['1', ] = runif(r, min = .95, max = 1)
#   
#   repeat{
#     if(length(Q) == 0) break;
#     
#     x = Q[1] # node
#     Q = Q[-1]
# 
#     ch = children(DataFrameToMatrix(tree), x)
#     
#     if(length(ch) == 0) next;
#     
#     
#     for(s in 1:r) {
#       x.CCF = data[x, s]
#       ch = sample(ch, length(ch))
#       chR = sample(ch, sample(1:length(ch)))
#       
#       # print(chR)
#       
#       # ccf = runif(length(chR), min = max(0, x.CCF - .5), max = x.CCF-0.4)/length(chR)
#       ccf = runif(length(chR), min = 0.02, max = x.CCF/(length(chR) + identify))
#       
#       data[chR, s] = ccf
#     }
#     
#     ch = sample(ch, sample(1:length(ch)))
#     
#     Q = c(Q, ch)
#   }
#   data[is.na(data)] = 0
#   data[data <= 0.02] = 0
#   
#   
#   data = round(data, 2)
# 
#   missing = rownames(data[rowSums(data) < 0.02, ])
#   data = data[rowSums(data) >= 0.02, ]
# 
#   tree = tree[!(tree$from %in% missing), , drop = FALSE]
#   tree = tree[!(tree$to %in% missing), , drop = FALSE]
#   
#   # tree = rbind(data.frame(from = 'GL', to = '1', stringsAsFactors = FALSE), tree)
# 
#   # Apply noise to CCF
#   sample.data = lapply(
#     noise,
#     function(n) {
#       
#       if(n > 0) {
#         v = apply(data, c(1,2), 
#                   function(w){ 
#                     if(w > 0.02) {
#                       r = rnorm(1, mean = 0, sd = n)
#                       if(w + r < 0.02) w = 0.02
#                       else w = w + r
#                     }
#                     w 
#                     } )
#       }
#       else v = data
# 
#       for(c in 1:ncol(v)) if(v[1, c] < max(v[, c])) { v[1, c] = max(v[, c]) + 0.01;  v[, c] = v[, c] / max(v[, c])  }
#       
#       v = round(v, 2)
#       
#       v[v > 1] = 1
#       v[v < 0.02] = 0
#       
#       v
#     }
#   )
#   names(sample.data) = noise
#   # print(sample.data)
#   
# 
#   return(list(tree = tree, CCF = data, sample.data = sample.data))
# }


# It is ambiguous a linear model with CCF that resambles branching, and the other way round.
# However, generating linear confounding models is much easier. We do this here.
rev_generate_true_model = function(n, r, noise = c(0, 0.02, 0.05, 0.1), identify = TRUE)
{

  hops = expand.grid(from = 'GL', to = 1, stringsAsFactors = FALSE)
  
  for(i in 2:n){
    hops = rbind(hops,
                 data.frame(from = i-1, to = i, stringsAsFactors = FALSE)
    )
  }
  # plot(graph_from_adjacency_matrix(DataFrameToMatrix(hops)))
  hops = hops[hops$from != 'GL', , drop = FALSE]
  
  tree = hops
  clones = unique(unlist(tree))
  
  data = matrix(0, nrow = length(clones), ncol = r)
  colnames(data) = paste('R', 1:r, sep = '')
  rownames(data) = clones
  
  data['1', ] = runif(r, min = .95, max = 1)

  for(reg in colnames(data)){
    for(idx in 1:nrow(data))
    {
      if(idx == 1) next; 
      
      if(!identify)
      {
        p = data[idx-1, reg] / (nrow(data) - idx + 1)
        if(0.01 * idx < p - 0.02)
          data[idx, reg] = runif(1, min = 0 + 0.01 * idx, p - 0.02)
        else data[idx, reg] = 0
      }
      else
      {
        p = data[idx-1, reg] / (nrow(data) - idx + 1)
        if(idx == nrow(data))  data[idx, reg] = runif(1, min = 0, max = data[idx-1, reg] - 0.02)
        else data[idx, reg] = runif(1, min = p + 0.02, max = data[idx-1, reg] - 0.02)
      }
    }
  }

  data[is.na(data)] = 0
  data[data <= 0.02] = 0
  
  data = round(data, 2)
  
  # print(tree)  
  # print(data)
  # print(missing)
  
  missing = rownames(data[rowSums(data) < 0.02, ])
  data = data[rowSums(data) >= 0.02, , drop = FALSE]

  # print(tree)  
  tree = tree[tree$from %in% rownames(data), , drop = FALSE]
  tree = tree[tree$to %in% rownames(data), , drop = FALSE]
  # print(tree)  
  
  # tree = rbind(data.frame(from = 'GL', to = '1', stringsAsFactors = FALSE), tree)
  
  # Apply noise to CCF
  sample.data = lapply(
    noise,
    function(n) {
      
      if(n > 0) {
        v = apply(data, c(1,2), 
                  function(w){ 
                    if(w > 0.02) {
                      r = rnorm(1, mean = 0, sd = n)
                      if(w + r < 0.02) w = 0.02
                      else w = w + r
                    }
                    w 
                  } )
      }
      else v = data
      
      for(c in 1:ncol(v)) if(v[1, c] < max(v[, c])) { v[1, c] = max(v[, c]) + 0.01;  v[, c] = v[, c] / max(v[, c])  }
      
      v = round(v, 2)
      
      v[v > 1] = 1
      v[v < 0.02] = 0
      
      v
    }
  )
  names(sample.data) = noise
  # print(sample.data)
  
  
  return(list(tree = tree, CCF = data, sample.data = sample.data))
}


rev_generate_denovo_sequencing = function(n.confounded, n.total, b, r, noise, 
                                          n.truncalMuts = 1,
                                          n.subclonalMuts = 30,
                                          n.childrenPerMut = 3,
                                          n.driverPerCloneMuts = 3
                                          )
{
  driversTree = NULL
  subcl.tree = make_tree(n.subclonalMuts, children = n.childrenPerMut) 
  subcl.tree = set.vertex.attribute(subcl.tree, 'name', value = paste('S', 1:n.subclonalMuts, sep =''))
    
  cl.trunk = make_tree(n.truncalMuts, children = 1)
  cl.trunk = set.vertex.attribute(cl.trunk, 'name', value = paste('C', 1:n.truncalMuts, sep =''))
  cl.trunk = set.vertex.attribute(cl.trunk, 'name', value = paste('S', 1, sep =''), index = n.truncalMuts)
    
  driversTree = cl.trunk + subcl.tree
    
  random_code = paste(LETTERS[sample(1:15, 3, replace = T)], collapse = '')
  driversTree = set.vertex.attribute(driversTree, "name", value = paste(random_code, 1:vcount(driversTree), sep =''))
  
  # Longest path in this tree
  driversTree.p = all_simple_paths(driversTree, 
                          from = V(driversTree)$name[1], 
                          to = leaves(as_adj(driversTree, sparse = FALSE)))
  
  ldriversTree = sapply(driversTree.p, function(w) length(w$name))
  mldriversTree = max(ldriversTree)
  
  # A path in T1 to a path in driversTree
  morph = function(T1){
    
    # print(T1)
    # repeat{
      GT1 = graph_from_adjacency_matrix(DataFrameToMatrix(T1))
      T1.p = all_simple_paths(GT1, 
                              from = V(GT1)$name[1], 
                              to = leaves(as_adj(GT1, sparse = FALSE)))
      
      lT1 = sapply(T1.p, function(w) length(w$name))
      T1.p = T1.p[which(lT1 == max(lT1))]
      
      driversTree.p = all_simple_paths(driversTree, 
                              from = V(driversTree)$name[1], 
                              to = leaves(as_adj(driversTree, sparse = FALSE)))
      
      maxL = max(unlist(lapply(T1.p, function(x) length(x$name))))
      
      T1.p = T1.p[[sample(1:length(T1.p), 1)]]$name
      driversTree.p = driversTree.p[[sample(1:length(driversTree.p), 1)]]$name
      
      mapping = NULL
      idx = idxD = 1
      repeat {
        if(length(T1.p) < idx) break;
        if(length(driversTree.p) < idxD) break;
        
        # if(idx > 1 && runif(1) > .2) {idx = idx + 1; next;}
        
        attach = sample(1:n.driverPerCloneMuts, 1)
        attach = idxD:(idxD + attach-1) 
        attach = attach[attach <= length(driversTree.p)]
        # cat(idxD, '  ', attach, '\n')
  
        mapping = append(mapping, list(driversTree.p[attach]))
        names(mapping)[length(mapping)] = T1.p[idx]
        
        idxD = max(attach) + 1
        # idx = idx + 2
        idx = idx + 1
      }

     # if(length(mapping) > min(2, maxL-1)) break;
    # }
    
    mapping
  }
  
  # Confounded models
  cat(cyan('Generate and morph'), n.confounded, cyan('confounded models: '))
  true_models = mappings = NULL
  for(i in 1:n.confounded){
    repeat{
      suppressWarnings({model = rev_generate_true_model(mldriversTree, r, noise, identify = FALSE)})
      
      model.p = graph_from_adjacency_matrix(DataFrameToMatrix(model$tree))
      model.p = all_simple_paths(model.p, 
                                 from = V(model.p)$name[1], 
                                 to = leaves(as_adj(model.p, sparse = FALSE)))
     
      lmodel.p = max(sapply(model.p, function(w) length(w$name)))

      # print(lmodel.p)
      
      if(lmodel.p < mldriversTree) next;
      
      
      map =  morph(model$tree) 
      # print(map)
      
      

      if(length(map) >= 3) break;
    }
    mappings = append(mappings, list(map))
    true_models = append(true_models, list(model))
    cat('+')
  }
  
  # NON Confounded models
  n.nconfounded = n.total - n.confounded
  cat(cyan('\nGenerate and morph'), n.nconfounded, cyan(' non-confounded models: '))
  for(i in 1:n.nconfounded){
    repeat{
      
      model = rev_generate_true_model(mldriversTree + 2 , r, noise, identify = TRUE)
      
      model.p = graph_from_adjacency_matrix(DataFrameToMatrix(model$tree))
      model.p = all_simple_paths(model.p, 
                                 from = V(model.p)$name[1], 
                                 to = leaves(as_adj(model.p, sparse = FALSE)))
      
      lmodel.p = max(sapply(model.p, function(w) length(w$name)))
      
      if(lmodel.p < mldriversTree) next;
      
      map =  morph(model$tree) 

      if(length(map) >= 3) break;
    }
    mappings = append(mappings, list(map))
    true_models = append(true_models, list(model))
    cat('+')
  }

  names(true_models) = paste('SynthPatient', 1:length(true_models), sep = '-')
  names(mappings) = names(true_models)
  
  # Transform everything in a dataframe  
  generate_dataset = function(ccf, attachment, patient) 
  {
    df = data.frame(matrix(ncol = 7, nrow = 0), stringsAsFactors = FALSE)
    colnames(df) = c('Misc', 'patientID', 'variantID', 'cluster', 'is.driver', 'is.clonal', 'CCF')
    
    clones = rownames(ccf)
    
    ccf = sapply(colnames(ccf), function(c) paste(c, ccf[, c], sep = ':'))
    ccf = apply(ccf, 1, paste, collapse = ';')
    names(ccf) = clones
    
    for(c in clones) {
      entry = data.frame(
        Misc = '--', 
        patientID = patient, 
        variantID = 'Dontcare', 
        cluster = c, 
        is.driver = FALSE, 
        is.clonal = ifelse(c == '1', TRUE, FALSE), 
        CCF = ccf[c], 
        stringsAsFactors = FALSE)
      df = rbind(df, entry)
    }
    
    for(c in 1:length(attachment)) 
    {
      for(v in attachment[[c]]){
        entry = data.frame(
          Misc = '--', 
          patientID = patient, 
          variantID = v, 
          cluster = names(attachment)[c], 
          is.driver = TRUE, 
          is.clonal = ifelse(c == '1', TRUE, FALSE), 
          CCF = ccf[c],
          stringsAsFactors = FALSE)
        df = rbind(df, entry)
      }
    }
    
    df
  }
  
  datasets = NULL
  for(n in noise) {
    cat('\nDataset noise @', n)
    d = lapply(names(true_models), 
           function(p)
             generate_dataset(true_models[[p]]$sample.data[[as.character(n)]], mappings[[p]], p))
    datasets = append(datasets, list(Reduce(rbind, d)))
    cat('OK')
    
  }
  names(datasets) = noise

  return(list(datasets = datasets, true.models = true_models, driversTree = driversTree))
}


rev_batch_denovoTest = function(seed = 12345,
                                n.confounded = 6, n.total = 10, b = 8, r = 3, 
                                noise = c(0, 0.01, 0.05, 0.1), 
                                initial.solution = NA,
                                n.truncalMuts = 1,
                                n.subclonalMuts = 30,
                                n.childrenPerMut = 3,
                                n.driverPerCloneMuts = 1)
{
  if(!is.null(seed)) set.seed(seed)
  
  cat(bgRed('rev_batch_denovoTest\n\n'), red(
    'NUMBER OF CONFOUNDED', n.confounded, '\n',
    'TOTAL', n.total, '\n',
    'REGIONS', r, '\n',
    'MUTS PER CLUSTER', n.driverPerCloneMuts, '\n'
  ))
  
  # Generate cohort
  denovo = rev_generate_denovo_sequencing(n.confounded = n.confounded, n.total = n.total, 
                                          b = b, r = r, noise = noise, 
                                          n.truncalMuts = n.truncalMuts,
                                          n.subclonalMuts = n.subclonalMuts,
                                          n.childrenPerMut = n.childrenPerMut,
                                          n.driverPerCloneMuts = n.driverPerCloneMuts)
  
  print(denovo)

  # Some auxiliary functions
  TRACERx.CCF.parser = function(x)
  {
    tk = strsplit(x, ';')[[1]]
    tk = unlist(strsplit(tk, ':'))
    
    samples = tk[seq(1, length(tk), 2)]
    
    values = tk[seq(2, length(tk), 2)]
    names(values) = samples
    
    return(values)  
  }
  
  ### Recurrence threshold
  dataset = denovo$datasets[[1]]

  data.split = dataset[dataset$is.driver, ]
  data.split = split(data.split, f = data.split$variantID)
  head(data.split)
  
  occurrencesCount = lapply(data.split, function(x) unique(x['patientID']))
  recurrentDrivers = occurrencesCount[unlist(lapply(occurrencesCount, function(x) nrow(x) >= 2))]
  
  variantIDs.recurrentDrivers = names(recurrentDrivers)
  driversToSkip = setdiff(names(occurrencesCount), variantIDs.recurrentDrivers)
  
  cat(red('Below recurrence threshold : \n'))
  print(driversToSkip, quote = FALSE)
  
  cat(green('Above recurrence threshold : \n'))
  print(variantIDs.recurrentDrivers, quote = FALSE)
  
  exit.codes = NULL
  for(ns in noise) {
    
    ns = as.character(ns)
    # print(ns)
    # print(names(denovo$datasets))
    
    dataset = denovo$datasets[[ns]]

    # We change the DRIVER annotation to FALSE for all entries in "driversToSkip"
    dataset[dataset$variantID %in% driversToSkip, 'is.driver'] = FALSE

    cohort = revolver_cohort(
      dataset = dataset,
      CCF.parser = TRACERx.CCF.parser,
      annotation = 'Synthetic cohort',
      options = list(
        ONLY.DRIVER = FALSE, 
        MIN.CLUSTER.SIZE = 0)
    )

  if(revolver_check_cohort(cohort, return.value = TRUE)) stop('Skewed... aborting... ')
  
  for(patient in cohort$patients) {
    cohort = compute_rev_phylogenies(
      cohort, 
      patient, 
      options = list(sspace.cutoff = 2000,
                     n.sampling = 500,
                     store.max = 100,
                     overwrite = FALSE))
    save(cohort, file='temp.RData')
  }

  fit = revolver_fit(cohort, initial.solution = initial.solution, restarts = 15, transitive.orderings = T, verbose = F)
  
  result = list(cohort = cohort, 
                fit = fit, denovo = denovo,
                n.confounded = n.confounded, 
                n.total = n.total, 
                b = b, 
                r = r, 
                noise = ns,
                n.driverPerCloneMuts = n.driverPerCloneMuts)
  
  performance = rev_analyze_experiment_denovo(result)
  result$performance = performance
  
  code = paste(sample(LETTERS, 8), collapse = '')
  exit.codes = c(exit.codes, code)
  
  save(result, file = paste('Experiment', code, '.Rdata', sep = '-'))
  }
  
  return(exit.codes)
}

rev_analyze_experiment_denovo = function(result)
{
  fit = result$fit
  patients = fit$patients
  
  tO = as_adj(result$denovo$driversTree, sparse = FALSE)
  tO = rbind(
    MatrixToDataFrame(tO), 
    data.frame(from = 'GL', to = root(tO), stringsAsFactors = FALSE)
  )

  whichIs = function(p){
    key = sort(DataFrameToEdges(result$denovo$true.models[[p]]$tree))
    
    bool = sapply(fit$phylogenies[[p]], function(w) {
      w = MatrixToDataFrame(w$adj_mat)
      w = w[w$from != 'GL', , drop = F]
      w = sort(DataFrameToEdges(w))
      all(w == key)
    })
    
    if(all(!bool)) return('*')
    bool = as.integer(bool)
    which(bool == 1, arr.ind = TRUE)
  }
  
  performance.fit = data.frame(
    fitID = fit$fit$solutionID,
    trueID = sapply(patients, whichIs))
  
  performance.fit$NumPhylo = sapply(patients, function(w){length(fit$phylogenies[[w]])})
  performance.fit$NumTR = sapply(patients, rev_count_information_transfer_comb, x = fit)
  performance.fit$NumDrv = sapply(patients, function(w) {
    r = fit$phylogenies[[w]][[1]]
    nrow(r$dataset[r$dataset$is.driver, ])
    })
  
  w =  lapply(
    1:nrow(performance.fit),
    function(w) {
      patient = rownames(performance.fit)[w]
      
      fT = fit$fit$phylogenies[[patient]]$transfer$driver
      tOfT = tO[tO$from %in% unlist(fT), ]
      tOfT = tOfT[tOfT$to %in% unlist(fT), ]
      
      TP = length(intersect(DataFrameToEdges(fT), DataFrameToEdges(tOfT)))
      FP = length(setdiff(DataFrameToEdges(fT), DataFrameToEdges(tOfT)))
      FN = length(setdiff(DataFrameToEdges(tOfT), DataFrameToEdges(fT)))
      
      return(data.frame(
                  TLdrvTP = TP/nrow(fT),
                  TLdrvFP = FP/nrow(fT),
                  TLdrvFN = FN/nrow(fT),
                  stringsAsFactors = FALSE
                ))
              })
  
  performance.fit = cbind(performance.fit, Reduce(rbind, w))
  
  w =  lapply(
    1:nrow(performance.fit),
    function(w) {
      patient = rownames(performance.fit)[w]
      
      fT = fit$phylogenies[[patient]][[1]]$transfer$driver
      tOfT = tO[tO$from %in% unlist(fT), ]
      tOfT = tOfT[tOfT$to %in% unlist(fT), ]
      
      TP = length(intersect(DataFrameToEdges(fT), DataFrameToEdges(tOfT)))
      FP = length(setdiff(DataFrameToEdges(fT), DataFrameToEdges(tOfT)))
      FN = length(setdiff(DataFrameToEdges(tOfT), DataFrameToEdges(fT)))
      
      return(data.frame(
        NTLdrvTP = TP/nrow(fT),
        NTLdrvFP = FP/nrow(fT),
        NTLdrvFN = FN/nrow(fT),
        stringsAsFactors = FALSE
      ))
    })
  
  performance.fit = cbind(performance.fit, Reduce(rbind, w))
  
  list(performance.fit,
       overall = list(
         TLdrvTP = mean(performance.fit$TLdrvTP),
         TLdrvFP = mean(performance.fit$TLdrvFP),
         TLdrvFN = mean(performance.fit$TLdrvFN),
         NTLdrvTP = mean(performance.fit$NTLdrvTP),
         NTLdrvFP = mean(performance.fit$NTLdrvFP),
         NTLdrvFN = mean(performance.fit$NTLdrvFN),
         NumPhylo = mean(performance.fit$NumPhylo),
         NumTR = mean(performance.fit$NumTR),
         NumDrv = mean(performance.fit$NumDrv),
         NonTopFit = sum(as.integer(performance.fit$fitID != 1))/nrow(performance.fit)
       )
  )
  
}



########### Clustering functions
rev_batch_denovoTest = function(seed = 12345,
                                n.confounded = 6, n.total = 10, b = 8, r = 3, 
                                noise = c(0, 0.01, 0.05, 0.1), 
                                initial.solution = NA,
                                n.truncalMuts = 1,
                                n.subclonalMuts = 30,
                                n.childrenPerMut = 3,
                                n.driverPerCloneMuts = 1)
{
  if(!is.null(seed)) set.seed(seed)
  
  cat(bgRed('rev_batch_denovoTest\n\n'), red(
    'NUMBER OF CONFOUNDED', n.confounded, '\n',
    'TOTAL', n.total, '\n',
    'REGIONS', r, '\n',
    'MUTS PER CLUSTER', n.driverPerCloneMuts, '\n'
  ))
  
  # Generate cohort
  denovo = rev_generate_denovo_sequencing(n.confounded = n.confounded, n.total = n.total, 
                                          b = b, r = r, noise = noise, 
                                          n.truncalMuts = n.truncalMuts,
                                          n.subclonalMuts = n.subclonalMuts,
                                          n.childrenPerMut = n.childrenPerMut,
                                          n.driverPerCloneMuts = n.driverPerCloneMuts)
  
  print(denovo)
  
  # Some auxiliary functions
  TRACERx.CCF.parser = function(x)
  {
    tk = strsplit(x, ';')[[1]]
    tk = unlist(strsplit(tk, ':'))
    
    samples = tk[seq(1, length(tk), 2)]
    
    values = tk[seq(2, length(tk), 2)]
    names(values) = samples
    
    return(values)  
  }
  
  ### Recurrence threshold
  dataset = denovo$datasets[[1]]
  
  data.split = dataset[dataset$is.driver, ]
  data.split = split(data.split, f = data.split$variantID)
  head(data.split)
  
  occurrencesCount = lapply(data.split, function(x) unique(x['patientID']))
  recurrentDrivers = occurrencesCount[unlist(lapply(occurrencesCount, function(x) nrow(x) >= 2))]
  
  variantIDs.recurrentDrivers = names(recurrentDrivers)
  driversToSkip = setdiff(names(occurrencesCount), variantIDs.recurrentDrivers)
  
  cat(red('Below recurrence threshold : \n'))
  print(driversToSkip, quote = FALSE)
  
  cat(green('Above recurrence threshold : \n'))
  print(variantIDs.recurrentDrivers, quote = FALSE)
  
  exit.codes = NULL
  for(ns in noise) {
    
    ns = as.character(ns)
    # print(ns)
    # print(names(denovo$datasets))
    
    dataset = denovo$datasets[[ns]]
    
    # We change the DRIVER annotation to FALSE for all entries in "driversToSkip"
    dataset[dataset$variantID %in% driversToSkip, 'is.driver'] = FALSE
    
    cohort = revolver_cohort(
      dataset = dataset,
      CCF.parser = TRACERx.CCF.parser,
      annotation = 'Synthetic cohort',
      options = list(
        ONLY.DRIVER = FALSE, 
        MIN.CLUSTER.SIZE = 0)
    )
    
    if(revolver_check_cohort(cohort, return.value = TRUE)) stop('Skewed... aborting... ')
    
    for(patient in cohort$patients) {
      cohort = compute_rev_phylogenies(
        cohort, 
        patient, 
        options = list(sspace.cutoff = 2000,
                       n.sampling = 500,
                       store.max = 100,
                       overwrite = FALSE))
      save(cohort, file='temp.RData')
    }
    
    fit = revolver_fit(cohort, initial.solution = initial.solution, restarts = 15, transitive.orderings = T, verbose = F)
    
    result = list(cohort = cohort, 
                  fit = fit, denovo = denovo,
                  n.confounded = n.confounded, 
                  n.total = n.total, 
                  b = b, 
                  r = r, 
                  noise = ns,
                  n.driverPerCloneMuts = n.driverPerCloneMuts)
    
    performance = rev_analyze_experiment_denovo(result)
    result$performance = performance
    
    code = paste(sample(LETTERS, 8), collapse = '')
    exit.codes = c(exit.codes, code)
    
    save(result, file = paste('Experiment', code, '.Rdata', sep = '-'))
  }
  
  return(exit.codes)
}

