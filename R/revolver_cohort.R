#' Construct a REVOLVER cohort object (S3 class \code{"rev_cohort"}).
#'
#' @param dataset A dataframe in the specified format (see Online manual).
#' @param CCF.parser A function to parse the format for the encoding of CCF
#' or binary values for each sequenced region. A possible function is available
#' inside REVOLVER; since it is not exported but is available with
#' \code{revolver:::CCF.parser} (the default of this parameter).
#' @param options A list of 2 parameters that should be a boolean value for
#' \code{ONLY.DRIVER} (use only driver SNVs), and \code{MIN.CLUSTER.SIZE}, the minimum cluster size.
#' @param annotation String for annotation of this cohort. This will be prompted
#'                   in every print for this object.
#'
#' @return An object of class \code{"rev_cohort"}
#'
#' @aliases revolver_cohort
#'
#' @examples
#' data(CRC)
#' cohort = revolver_cohort(CRC, options = list(ONLY.DRIVER = FALSE, MIN.CLUSTER.SIZE = 0))
#'
#' @import crayon
#'
#' @export
revolver_cohort = function(dataset,
                           CCF.parser = revolver:::CCF.parser,
                           options = list(ONLY.DRIVER = FALSE, MIN.CLUSTER.SIZE = 10),
                           annotation = 'My REVOLVER dataset')
{
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Some header print and the output object are created, before checking for input consistency
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  
  pio::pioHdr(
    paste('REVOLVER Cohort constructor'),
    c(
      `Use only alterations annotated as driver` = options$ONLY.DRIVER,
      `Filter: minimum number of alterations in a cluster` = options$MIN.CLUSTER.SIZE
    ),
    prefix = '\t'
  )
  
  # The output object will be this, of class rev_cohort
  obj <-
    structure(
      list(
        patients = NULL,
        variantIDs = NULL,
        variantIDs.driver = NULL,
        numVariants = NULL,
        annotation = annotation,
        dataset = NULL,
        CCF = NULL,
        n = NULL,
        CCF.parser = CCF.parser
      ),
      class = "rev_cohort",
      call = match.call()
    )
  
  # Check input and stop on error
  dataset = check_input(dataset, CCF.parser)
  dataset$id = paste0('__mut_id_', 1:nrow(dataset))
  
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Input data is subset according to the options parameter
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if (options$ONLY.DRIVER)
    dataset = dataset %>%
    filter(is.driver)
  
  # by cluster size
  grp = dataset %>%
    group_by(patientID, cluster) %>%
    summarize(cluster_size = n()) %>%
    ungroup()
  
  dataset = left_join(dataset, grp, by = c('patientID', 'cluster'))
  
  if (options$MIN.CLUSTER.SIZE > 0)
  {
    pio::pioTit('Checking size of each patient\'s cluster -- removing those below cutoff.')
    
    cat('Rows before filtering:', nrow(dataset), '\n\n')
    
    cat('Removed\n\n')
    print(dataset %>%
            filter(cluster_size < options$MIN.CLUSTER.SIZE))
    
    dataset = dataset %>%
      filter(cluster_size >= options$MIN.CLUSTER.SIZE)
    
    cat('\nRows after filtering:', nrow(dataset), '\n')
  }
  
  if (nrow(dataset) == 0)
    stop("Your dataset is empty, aborting.")
  
  pio::pioTit('REVOLVER input data')
  pio::pioDisp(dataset)
  
  obj$patients = unique(dataset$patientID)
  obj$dataset = dataset
  obj$variantIDs = unique(dataset$variantID)
  obj$variantIDs.driver = unique(dataset %>% filter(is.driver) %>% pull(variantID))
  obj$numVariants = nrow(dataset)
  
  obj$n = list(
    patients = length(obj$patients),
    variants = nrow(obj$dataset),
    drivers = length(obj$variantIDs.driver)
  )
  
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Data of each patient needs to be unwrapped and stored
  # This takes some times
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  pio::pioTit('Extracting dataset for each patient (this may take some time)')
  
  obj$dataset = obj$dataset %>%
    group_by(patientID) %>%
    nest() %>%
    select(data) %>%
    unlist(recursive = F)
  names(obj$dataset) = obj$patients
  
  # For each patient extract its explicit region CCF data
  foo = function(w)
  {
    values = revolver:::CCF.parser(w)
    tv = as_tibble(as.numeric(values))
    tv$sample = names(values)
    
    tv %>%
      spread(sample, value)
  }
  
  obj$dataset = lapply(obj$patients,
                       function(pat)
                       {
                         d = obj$dataset[[pat]] %>%
                           mutate(patientID = pat) %>%
                           select(id,
                                  Misc,
                                  patientID,
                                  variantID,
                                  is.driver,
                                  is.clonal,
                                  cluster,
                                  cluster_size,
                                  CCF)
                         
                         pio::pioStr(d$patientID[1],
                                     paste0(nrow(d), ' entries'),
                                     suffix =  '\n',
                                     prefix = '\n')
                         
                         CCF_values = d %>%
                           group_by(id) %>%
                           do(foo(.$CCF)) %>%
                           ungroup()
                         
                         d %>%
                           full_join(CCF_values, by = 'id')
                       })
  names(obj$dataset) = obj$patients
  
  pio::pioTit('Extracting clones\' table for each patient')
  
  # Then we also create a table with the clones for each patient
  obj$CCF = lapply(names(obj$dataset),
                   function(pat)
                   {
                     d = obj$dataset[[pat]]
                     
                     pio::pioStr(pat, paste0(nrow(d), ' entries'), suffix =  '')
                     
                     d = d %>%
                       select(-Misc, -variantID, -CCF, -id, -patientID, -cluster_size)
                     
                     CCFclone_headers = d %>%
                       group_by(cluster) %>%
                       summarise(
                         nMuts = n(),
                         is.driver = any(is.driver),
                         is.clonal = all(is.clonal)
                       ) %>%
                       arrange(desc(nMuts))
                     
                     CCFclone_values = d %>%
                       select(-is.driver, -is.clonal) %>%
                       reshape2::melt(id = 'cluster') %>%
                       group_by(cluster, variable) %>%
                       summarise(CCF = median(as.numeric(value))) %>%
                       spread(variable, CCF) %>%
                       ungroup()
                     
                     pio::pioStr(' --> ', paste0(nrow(CCFclone_headers), ' clones'), suffix =  '\n')
                     
                     CCFclone_headers %>% full_join(CCFclone_values, by = 'cluster')
                   })
  names(obj$CCF) = names(obj$dataset)
  
  
  
  return(obj)
}


# Check input format - fails on error
check_input = function(dataset, CCF.parser)
{
  # Dataset with certain columns required
  required.cols = c('Misc',
                    'patientID',
                    'variantID',
                    'cluster',
                    'is.driver',
                    'is.clonal',
                    'CCF')
  
  
  # check for types
  types = sapply(dataset, class)
  
  if (!is.data.frame(dataset) |
      !all(required.cols %in% colnames(dataset)) |
      types['Misc'] != 'character' |
      types['patientID'] != 'character' |
      types['variantID'] != 'character' |
      types['cluster'] != 'character' |
      types['CCF'] != 'character' |
      types['is.driver'] != 'logical' |
      types['is.clonal'] != 'logical')
  {
    print(str(dataset))
    
    stop(
      'Input `dataset` should be a dataframe in with the following columns:\n
      \t- Misc: any annotation (string),
      \t- patientID: a unique ID to identify a patient (string)
      \t- variantID: a unique ID to identify a variant across patients (string)
      \t- cluster: to putative clone in this patient\'s tumour where the variant is found (string)
      \t- is.driver: if the mutation is a driver to be correlated (logical)
      \t- is.clonal: if the mutation is a clonal (truncal) in this patient (logical)
      \t- CCF: Cancer Cell Fraction or binary annotation for this mutations across the biopsies of this patient\'s tumour (string)

      For more information see the manual: https://github.com/caravagn/revolver'
    )
  }
  
  if (any(dataset[, required.cols[-1]] == ''))
    stop("Entries empty string ('') in the dataset are allowed only in the `Misc` field.")
  
  if (!is.function(CCF.parser))
  {
    print(revolver:::CCF.parser)
    stop(
      'You need to provide a function to parse CCFs, see "CCF.parser" available in this package.'
    )
  }
  
  dataset[, required.cols] %>% as_tibble()
}

#' Print a \code{"rev_cohort"} object
#'
#' @param x obj of class \code{"rev_cohort"}
#' @param digits number of output digits
#'
#' @return nothing
#' @export print.rev_cohort
#' @import crayon
#'
#' @examples
#' data(CRC.cohort)
#' CRC.cohort
print.rev_cohort <-
  function(x, ...)
  {
    stopifnot(inherits(x, "rev_cohort"))
    
    pio::pioHdr('REVOLVER - Repeated Evolution in Cancer',
                toPrint = NULL,
                suffix = '\n')
    
    pio::pioStr("Dataset :",
                x$annotation,
                suffix = '\n',
                prefix = '\n')
    
    pio::pioStr(
      "Cohort  :",
      paste0(
        x$n$patients,
        ' patients, ',
        x$n$variants,
        ' variants and ',
        x$n$drivers,
        ' driver events.'
      ),
      suffix = '\n\n'
    )
    
    # pio::pioTit('Available computations')
    
    
    pio::pioStr('Trees per patient    :',
                ifelse(is.null(x$phylogenies), red('NO'), green('YES')),
                suffix = '\n')
    
    pio::pioStr('Fit via TL           :',
                ifelse(is.null(x$fit), red('NO'), green('YES')), 
                suffix = '\n')
    
    pio::pioStr('REVOLVER clustering  :',
                ifelse(is.null(x$cluster), red('NO'), green('YES')),
                suffix = '\n')
    
    pio::pioStr('Jackknife statistics :',
        ifelse(is.null(x$jackknife), red('NO'), green('YES')),
        suffix = '\n')
    

    
    # pio::pioTit("Cohort summary statistics")
    # print(Stats(x))
    #
    # pio::pioTit("Drivers summary statistics")
    # print(Stats_drivers(x))
    
    
    pio::pioStr(
      "",
      "\nFor summary statistics see `?Stats_*(x)` with * = {cohort, drivers, trees, fits, clusters, ...}",
      suffix = '\n'
    )
    
    
    # if (!is.null(x$phylogenies))
    # {
    #   pio::pioTit('Summary for data and models', prefix = '\t')
    #
    #   # cat(cyan('\n\tTL Model Fit  :'), ifelse(is.null(x$fit), red('NO\n'), green('YES ')), '\n')
    #
    #   longest.name = max(nchar(names(x$phylogenies)))
    #
    #   patf = function(w, fit, clusters) {
    #     patient = names(x$phylogenies)[w]
    #
    #     s = paste(
    #       yellow(sprintf(
    #         paste('%', longest.name, 's', sep = ''), patient
    #       )),
    #       ': ',
    #       sprintf(
    #         'k = %3s | t = %3s | n = %2s | r = %2s | m = %4s | d = %2s',
    #         length(x$phylogenies[[w]]),
    #         rev_count_information_transfer_comb(x, names(x$phylogenies)[w]),
    #         x$phylogenies[[w]][[1]]$numNodes,
    #         x$phylogenies[[w]][[1]]$numRegions,
    #         nrow(x$phylogenies[[w]][[1]]$dataset),
    #         nrow(x$phylogenies[[w]][[1]]$dataset[x$phylogenies[[w]][[1]]$dataset$is.driver,])
    #       ),
    #       '\t'
    #     )
    #
    #     cat('\t', s)
    #
    #     if (fit)
    #     {
    #       stas = stats.rev_phylo(x$fit$phylogenies[[w]])
    #
    #       rank = x$fit$solutionID[w]
    #       score = format(x$fit$phylogenies[[w]]$score,
    #                      scientific = TRUE,
    #                      digits = 2)
    #       gofit = format(stas$gofit, scientific = TRUE, digits = 2)
    #
    #       cat(
    #         bgBlue('[ Fit ]', sprintf(' # %3s', yellow(rank)), ' '),
    #         '| g =',
    #         red(sprintf('%8s', gofit)),
    #         '| f =',
    #         green(sprintf('%8s', score)),
    #         '  '
    #       )
    #     }
    #
    #     if (clusters)
    #     {
    #       cluster = x$cluster$clusters[patient]
    #
    #       cat(bgMagenta('[ Cluster', sprintf('%3s', yellow(cluster)), ']'))
    #
    #     }
    #
    #     cat('\n')
    #   }
    #
    #   cat('\n')
    #   sapply(
    #     1:length(x$phylogenies),
    #     patf,
    #     fit = !is.null(x$fit),
    #     cluster = !is.null(x$cluster)
    #   )
    #
    #
    #   cat(
    #     '\n\tLegend \n\t\t k : phylogenies  \n\t\t t : combinations of information transfer \n\t\t n : groups (nodes of the tree) \n\t\t r : regions (inputs per patient) \n\t\t m : number of alterations \n\t\t d : number of driver alterations\n'
    #   )
    #
    #   if (!is.null(x$fit))
    #     cat(
    #       '\t\t # : number of the solution selection (out of k)  \n\t\t g : goodness-of-fit \n\t\t f : score of the model\n'
    #     )
    #
    #   if (!is.null(x$cluster))
    #   {
    #     pio::pioTit('Summary for clusters', prefix = '\t')
    #
    #     ids = names(x$cluster$labels.colors)
    #
    #     medians = rep(NA, length(ids))
    #
    #     if (!is.null(x$jackknife)) {
    #       medians = x$jackknife$cluster.medians
    #     }
    #
    #     size = as.vector(table(x$cluster$clusters))
    #
    #     df = data.frame(Cluster = ids,
    #                     n = size,
    #                     Jackknife = medians)
    #     apply(df, 1, function(w) {
    #       cat('\n\t',
    #           paste(
    #             yellow(sprintf('%4s', w['Cluster'])),
    #             ': ',
    #             sprintf('n = %4s | jackknife = %4s',
    #                     w['n'],
    #                     w['Jackknife'])
    #           ))
    #
    #     })
    #
    #     cat('\n\n\tLegend \n\t\t         n : number of patients in the cluster')
    #
    #     if (!is.null(x$jackknife))
    #       cat('\n\t\t jackknife : median co-clustering probability estiamted with jackknife\n')
    #
    #     # assignments = x$cluster$clusters
    #     # edges = x$jackknife$edges
    #     # pio::pioTit('Summary for edges', prefix = '\t')
    #     # colnames(edges)[3] = 'Jackknife'
    #     # pio::pioDisp(edges)
    #   }
    # }
    
    cat('\n')
    revolver_check_cohort(x)
    
  }



#' Cohort summary scatterplot.
#'
#' @details Returns a scatterplot of the tumour mutational burden. at the
#' clonal and subclonal level. Each dot is sized by the number of
#' clones with drivers, and coloured by the numnber of drivers.
#'
#' @param x A REVOLVER cohort.
#' @param cex Cex of the plot.
#'
#' @return A \code{ggplot2} object of the plot.
#' @export
#'
#' @examples
#' data(Breast.fit)
plot.rev_cohort = function(x, cex = 1, ...)
{
  ggplot(
    Stats(x),
    aes(
      color = numDriverMutations,
      size = numClonesWithDriver,
      x = numTruncalMutations,
      y = numSubclonalMutations
    )
  ) +
    stat_density_2d(size = .2 * cex, color = 'gray') +
    geom_point(alpha = .8) +
    geom_rug(size = .3 * cex) +
    scale_color_distiller(palette = 'Spectral') +
    theme_minimal(base_size = 10 * cex) +
    theme(legend.position = 'bottom',
          legend.key.size = unit(3 * cex, 'mm')) +
    scale_y_log10() +
    scale_x_log10() +
    coord_cartesian(clip = 'off') +
    labs(
      title = x$annotation,
      subtitle = paste("Cohort summary"),
      x = "Truncal mutations",
      y = "Subclonal mutations"
    ) +
    guides(
      color = guide_colorbar("Drivers", barwidth  = unit(3 * cex, 'cm')),
      size = guide_legend("Clones with drivers", barwidth  = unit(3 * cex, 'cm'))
    )
  
}