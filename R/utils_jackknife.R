analyse_jackknife = function(x)
{
  results = x$jackknife$results
  resamples = x$jackknife$params$resamples
  leave.out = x$jackknife$params$leave.out
  
  npatients = x$n$patients
  patients = x$patients
  
  #Co-clustering probability
  pio::pioTit('Patients co-clustering frequency')
  
  # A function to compute co-occurrences of samples in the clustering's results
  cumulative_co_occurrences = function(x, M) 
  {
    cluster.labels = Cluster(x) %>% pull(cluster) %>% unique
    
    for(cluster in cluster.labels) 
    {
      membership = Cluster(x) %>%
        filter(cluster == !!cluster) %>%
        pull(patientID)
      
      print(membership)
      print(Cluster(x) %>%
              filter(cluster == !!cluster))
      
      pairs = combn(membership, 2, simplify = F)
      
      for(p in 1:length(pairs)) {
        M[pairs[[p]][1], pairs[[p]][2]] = M[pairs[[p]][1], pairs[[p]][2]] + 1
        M[pairs[[p]][2], pairs[[p]][1]] = M[pairs[[p]][2], pairs[[p]][1]] + 1
      }
    }
    M
  }
  
  co.clustering = matrix(0, nrow = npatients, ncol = npatients)
  rownames(co.clustering) = colnames(co.clustering) = patients
  
  # Cumulative computation with prograss bar
  pb = dplyr::progress_estimated(length(results), 3)
  
  for(p in 1:length(results)) {
    pb$tick()$print()
    
    co.clustering = cumulative_co_occurrences(results[[p]], co.clustering)
  }
  
  co.clustering %>% pioDisp
  
  # Store normalized co-clustering probability
  x$jackknife$co_clustering = co.clustering/resamples
  
  # median jackknife statistics per cluster
  median_stats = c()
  
  pio::pioTit('Median per cluster')
  
  clusters = Cluster(x) %>% pull(cluster) %>% unique
  
  for(cluster in clusters) 
  {
    members = Cluster(x) %>% filter(cluster == !!cluster) %>% pull(patientID)
    median_stats[cluster] = median(unlist(x$jackknife$co_clustering[members, members]))
  }
  
  print(median_stats)
  
  x$jackknife$co_clustering_medians = median_stats
  
  ##################### EDGE ESTIMATORS
  pio::pioTit('Edge frequency across resamples')
  
  # From every resample, probability of one trajectory
  trajectories = lapply(x$jackknife$results,
                        function(w) {
                          # trajectories in this resample
                          # from the IT, top rankings
                          this_patient = lapply(
                            w$patients,
                            ITransfer,
                            x = x,
                            rank = 1,
                            type = 'drivers',
                            data = 'trees'
                          )
                          # Pooled and made unique
                          Reduce(bind_rows, this_patient) %>% 
                            distinct()
                        })
  
  trajectories = Reduce(bind_rows, trajectories) %>%
    group_by(from, to) %>%
    summarise(prob_resamp = n()/resamples) %>%
    ungroup() %>%
    arrange(desc(prob_resamp))
  
  # From every resample, number of patients with the trajectory
  trajectories_counts = lapply(seq_along(x$jackknife$results),
                        function(w) {
                          this_resample = x$jackknife$results[[w]]
                          
                          # trajectories in this resample
                          # from the IT, top rankings
                          this_patient = lapply(
                            this_resample$patients,
                            ITransfer,
                            x = x,
                            rank = 1,
                            type = 'drivers',
                            data = 'trees'
                          )
                          
                          # Pooled and made unique
                          Reduce(bind_rows, this_patient) %>% 
                            group_by(from, to) %>%
                            summarise(n = n()/this_resample$n$patients) %>%
                            arrange(desc(n)) %>%
                            mutate(resample = w) %>%
                            ungroup()
                        })
  
  trajectories_counts = Reduce(bind_rows, trajectories_counts) %>%
    group_by(from, to) %>%
    summarise(num_patients = median(n)) %>%
    arrange(desc(num_patients))  %>%
    ungroup()
    
  # Bond tibble
  trajectories = trajectories %>%
    full_join(trajectories_counts, by = c('from', 'to'))
  
  trajectories %>% pioDisp
  
  x$jackknife$trajectories = trajectories
  
  x
}
