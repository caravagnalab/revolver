#' Construct a REVOLVER cohort.
#'
#' @param dataset A dataframe in the specified format (see the package vignettes).
#' @param CCF_parser A function to parse the format for the encoding of CCF
#' or binary values for each sequenced region. A possible function is available
#' inside REVOLVER; \code{\link{CCF_parser}} (the default of this parameter).
#' @param ONLY.DRIVER If true, uses only annotated driver events.
#' @param MIN.CLUSTER.SIZE Discard clusters that have less than this number of entries.
#' @param annotation Brief cohort description.
#'
#' @return An object of the S3 class \code{"rev_cohort"} that represents a \code{REVOLVER} cohort.
#'
#' @import tidyverse
#' @import tidygraph
#' @import pio
#' @import easypar
#' @import RColorBrewer
#' @import crayon
#' @import cluster
#' @import dendextend
#' @import dynamicTreeCut
#' @import igraph
#' @import ggpubr
#' @import ggrepel
#' @import clisymbols
#' @import evoverse.datasets
#'
#' @export
#' @family Cohort creation
#'
#' @examples
#' # Example cohort creation with the TRACERx data
#' data(TRACERx_data)
#'
#' # To speed up the process we use only 2 patients
#' TRACERx_data = TRACERx_data %>%
#'    filter(patientID %in% c('CRUK0001', 'CRUK0002'))
#'
#' cohort = revolver_cohort(TRACERx_data, annotation = 'A toy REVOLVER dataset')
#'
#' # The S3 print/ plot for this cohort
#' print(cohort)
#' plot(cohort)
revolver_cohort = function(dataset,
                           CCF_parser = revolver::CCF_parser,
                           ONLY.DRIVER = FALSE,
                           MIN.CLUSTER.SIZE = 10,
                           annotation = 'My REVOLVER dataset')
{
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Some header print and the output object are created, before checking for input consistency
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  options = list(ONLY.DRIVER = ONLY.DRIVER, MIN.CLUSTER.SIZE = MIN.CLUSTER.SIZE)

  pio::pioHdr(
    paste('REVOLVER ~ Cohort constructor'),
    toPrint = c(
      ` Use only drivers` = ifelse(options$ONLY.DRIVER == 0, TRUE, FALSE),
      `Reject clusters with size below` = options$MIN.CLUSTER.SIZE
    ),
    prefix = paste0(clisymbols::symbol$radio_on, ' ')
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
        CCF_parser = CCF_parser
      ),
      class = "rev_cohort",
      call = match.call()
    )

  # Check input and stop on error
  dataset = revolver:::check_input(dataset, CCF_parser)
  dataset$id = paste0('__mut_id_', 1:nrow(dataset))

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Input data is subset according to the options parameter
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if (options$ONLY.DRIVER)
    dataset = dataset %>% filter(is.driver)

  # by cluster size
  grp = dataset %>%
    group_by(patientID, cluster) %>%
    summarize(cluster_size = n(), .groups = 'keep') %>%
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
    ungroup() %>%
    select(data) %>%
    unlist(recursive = F)
  names(obj$dataset) = obj$patients

  # For each patient extract its explicit region CCF data
  foo = function(w)
  {
    values = CCF_parser(w)
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

  if (!is.function(CCF_parser))
  {
    stop(
      'You need to provide a function to parse CCFs, see "CCF_parser" available in this package.'
    )
  }

  dataset[, required.cols] %>% as_tibble()
}



