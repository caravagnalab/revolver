#' Plot the scores for a patient's trees.
#'
#' @description Plots a barplot where each tree is
#' a bar, with height proportional to the tree score and
#' colour to its particular combination of transfer.
#'
#' @param x A REVOLVER object.
#' @param patient The id of a patient.
#' @param transfer_palette A function that can be used to sample
#' an arbitrary number of colours to identify the transfers.
#'
#' @family Plotting functions
#'
#' @return A \code{ggplot} plot.
#' @export
#'
#' @examples
#' # Data released in the 'evoverse.datasets'
#' data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')
#'
#' plot_patient_trees_scores(TRACERx_NEJM_2017_REVOLVER, 'CRUK0001')
#' plot_patient_trees_scores(TRACERx_NEJM_2017_REVOLVER, 'CRUK0002')
plot_patient_trees_scores = function(x,
                             patient,
                             transfer_palette = distinct_palette_few,
                             ...)
{
  p = patient

    if(!has_patient_trees(x,p))
      stop("Patient ", p, " does not have trees, aborting.")

    n = length(x$phylogenies[[p]])


    scores = lapply(1:n,
                    function(w) {
                      # Score
                      score = Phylo(x, p, w, data = 'trees')$score

                      # Inf Transf - coded as string
                      IT =  ITransfer(x, p, w, type = 'drivers')
                      IT = paste(sort(ctree:::DataFrameToEdges(IT)), collapse = ':')

                      data.frame(
                        patient = p,
                        rank = w,
                        score = score,
                        IT = IT
                      )
                    })

    scores = Reduce(rbind, scores) %>%
      as_tibble()

    st = Stats_trees(x, p)
    nc = st$combInfTransf

    ggplot(scores,
           aes(y = score, x = rank, fill = IT)) +
      geom_bar(stat = 'identity') +
      guides(color = FALSE,
             fill = guide_legend(title = "Information Transfer", label = FALSE)) +
      scale_fill_manual(values = transfer_palette(nc)) +
      my_ggplot_theme() +
      labs(
        title = p,
        x = "Tree rank",
        y = "Tree score",
        subtitle = paste0(nc, ' combination(s) of transfer.')
      )
}


