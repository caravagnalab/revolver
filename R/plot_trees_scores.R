plot_trees_scores = function(x,
         patients = x$patients,
         transfer_palette = distinct_palette_few
         cex = 1)
{
  # plot of a single patient
  single_plot = function(p)
  {
    n = length(x$phylogenies[[p]])
    
    scores = lapply(1:n,
                function(w) {
                  data.frame(
                    patient = p,
                    rank = w,
                    score = Phylo(x, p, w)$score,
                    IT = paste(sort(DataFrameToEdges(
                      ITransfer(x, p, w, type = 'drivers')
                    )),
                    collapse = ':')
                  )
                })
    
    scores = Reduce(rbind, scores) %>% 
      as_tibble()
    
    factor(scores$IT)
    
    st = Stats_trees(x, p)
    nc = st$combInfTransf
    
    ggplot(scores,
           aes(y = score, x = rank, fill = IT)) +
      geom_bar(stat = 'identity') +
      guides(color = FALSE, fill = guide_legend(title = "Information Transfer", label = FALSE)) +
      theme_minimal(base_size = 10 * cex) +
      scale_fill_manual(values = transfer_palette(nc)) +
      theme(
        legend.position = 'bottom',
        legend.key.size = unit(3, 'mm')
      ) +
      labs(
        title = p,
        x = "Tree rank",
        y = "Tree score",
        subtitle = paste0(nc, ' combination(s) of transfer.')
      )
  }
  
  lapply(patients, single_plot)
}
