#' Plot a barplot of the scored trees for a patient.
#' 
#' @description Returns a barplot where each tree is
#' a bar, with height proportional to the tree score and
#' colour to its particular combination of transfer. If
#' the function is called on a vector of patient ids, the
#' result is a list of plot.
#'
#' @param x A REVOLVER cohort object
#' @param patients The plots 
#' @param transfer_palette A function that can be used to sample
#' an arbitrary number of colours to identify the transfers.
#' @param cex Cex of the plot.
#' @param ... Extra parameters, unused.
#'
#' @return A \code{ggplot} plot.
#' @export
#'
#' @examples
#' data(TRACERx_cohort)
#' 
#' plot_trees_scores(TRACERx_cohort, 'CRUK0002')
plot_trees_scores = function(x,
                             patient,
                             transfer_palette = distinct_palette_few,
                             cex = 1,
                             ...)
{
  # # plot of a single patient
  # single_plot = function(p)
  # {
  #   if(!has_patient_trees(x,p))
  #     stop("Patient ", p, " does not have trees, aborting.")
  #     
  #   n = length(x$phylogenies[[p]])
  #   
  #   scores = lapply(1:n,
  #                   function(w) {
  #                     data.frame(
  #                       patient = p,
  #                       rank = w,
  #                       score = Phylo(x, p, w, data = 'trees')$score,
  #                       IT = paste(sort(DataFrameToEdges(
  #                         ITransfer(x, p, w, type = 'drivers')
  #                       )),
  #                       collapse = ':')
  #                     )
  #                   })
  #   
  #   scores = Reduce(rbind, scores) %>%
  #     as_tibble()
  #   
  #   factor(scores$IT)
  #   
  #   st = Stats_trees(x, p)
  #   nc = st$combInfTransf
  #   
  #   ggplot(scores,
  #          aes(y = score, x = rank, fill = IT)) +
  #     geom_bar(stat = 'identity') +
  #     guides(color = FALSE,
  #            fill = guide_legend(title = "Information Transfer", label = FALSE)) +
  #     theme_minimal(base_size = 10 * cex) +
  #     scale_fill_manual(values = transfer_palette(nc)) +
  #     theme(legend.position = 'bottom',
  #           legend.key.size = unit(3, 'mm')) +
  #     labs(
  #       title = p,
  #       x = "Tree rank",
  #       y = "Tree score",
  #       subtitle = paste0(nc, ' combination(s) of transfer.')
  #     )
  # }
  # 
  # lapply(patients, single_plot)
  
  p = patient
  
    if(!has_patient_trees(x,p))
      stop("Patient ", p, " does not have trees, aborting.")
    
    n = length(x$phylogenies[[p]])
    
    scores = lapply(1:n,
                    function(w) {
                      data.frame(
                        patient = p,
                        rank = w,
                        score = Phylo(x, p, w, data = 'trees')$score,
                        IT = paste(sort(ctree:::DataFrameToEdges(
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
      guides(color = FALSE,
             fill = guide_legend(title = "Information Transfer", label = FALSE)) +
      theme_minimal(base_size = 10 * cex) +
      scale_fill_manual(values = transfer_palette(nc)) +
      theme(legend.position = 'bottom',
            legend.key.size = unit(3, 'mm')) +
      labs(
        title = p,
        x = "Tree rank",
        y = "Tree score",
        subtitle = paste0(nc, ' combination(s) of transfer.')
      )
}


