#' Plot the fit penalty
#' 
#' @description 
#' 
#' Plot the fit penalty as a barplot, for each one of a set of desired
#' driver events, where the bar represents the counts of each trajectory
#' in the data. This function allows also to filter out entries that have
#' been seen below a predetermined cutoff, and tests for significance in the
#' association A -> B via a one-sided Fisher 2x2 test adjusted for the number of
#' comparison (marginal count of B-ended trajectories). The tests are carried
#' out by function \code{revolver:::enrichment_test_incoming_edge}, which can
#' be used to obtain a tidy representation of the tests' results.
#'
#' @param x A REVOLVER cohort with fit models.
#' @param drivers The list of drivers to use; by default all of them. If the
#' entry is a subset of the actual list of all drivers, all the entries in the
#' penalty data structure \code{x$fit$penalty} will be used if they involve at
#' least one event from \code{drivers}.
#' @param min.occurrences The penalty data structure will be filtered for
#' \code{count} values above this threshold.
#' @param alpha_level The significance level for the enrichment Fisher test.
#' @param drivers_palette A function that can return, for an input number,
#' a number of colours.
#' @param cex The cex of the plot.
#' @param ... Other unused parameters.
#'
#' @return A ggplot object for this plot.
#' 
#' @export
#'
#' @examples
#' TODO
plot_penalty = function(x,
                        drivers = x$variantIDs.driver,
                        min.occurrences = 0,
                        alpha_level = 0.05,
                        drivers_palette = colorRampPalette(RColorBrewer::brewer.pal(n = 9, "Set1")),
                        cex = 1,
                        ...
                        )
{
  # TODO - check input x for fits

  # Penalty for the required nodes
  E = x$fit$penalty %>%
    filter(to %in% drivers, count >= min.occurrences)

  # Actual drivers
  drivers = unique(c(E$from, E$to))

  N = length(drivers)

  drivers_colors = c("black", drivers_palette(N))
  names(drivers_colors) = c("GL", drivers)

  # Enrichment test for association
  tests = enrichment_test_incoming_edge(
    E, 
    alpha_level)

  # Labels for significant hits
  tests_labels = tests %>%
    filter(psign) %>%
    group_by(to) %>%
    summarise(
      label = paste0(' ', from, ' (p = ', format(p.value, scientific = TRUE, digits = 2), ')', collapse = ', ')
    ) %>%
    left_join(
      # Positioning of the label at the end of each cumulative bar
      E %>% 
        group_by(to) %>%
        summarise(count = sum(count)),
      by = 'to'
    )
  
  # Plot
  ggplot(E,
         aes(y = count, x = to)) +
    geom_bar(aes(fill = from), stat = 'identity') +
    coord_flip() +
    theme_minimal(base_size = 10 * cex) +
    scale_fill_manual(values = drivers_colors) +
    theme(
      legend.position = 'bottom',
      legend.key.size = unit(3, 'mm')
    ) +
    labs(
      title = paste0("Counts for trajectories detected at least ", min.occurrences, " times"),
      subtitle =
        paste0("Enrichment via adjusted one-sided Fisher test at alpha level ", alpha_level),
      y = "Counts",
      x = "Drivers"
    ) +
    guides(fill = guide_legend(title = 'Ancestor', nrow = 2)) +
    geom_text(
      data = tests_labels,
      aes(label = label), 
      color = 'black', 
      size = 2 * cex, 
      hjust = 0
      )
}



# This fuction performs an enrichment test (Fisher one-sided)
# for each edge A -> B  testing if the amount of observations
# with A upstream B is significant. The tests are adjusted
# via Bonferroni for the number of distinct tests carried out
# for each possible pair A/B, with B fixed (i.e. the number of
# possible events evolutionary upstream of B).
# 
# All tests are summarised in a tibble via broom
enrichment_test_incoming_edge = function(E, alpha_level = 0.05)
{
  pio::pioTit("Enrichment test for incoming edges")
  
  # Single test
  single_test = function(from, to)
  {
    # 2 x 2 contingency table
    CTb =
      matrix(
        c(
          E %>% filter(from == !!from, to == !!to) %>% summarise(count = sum(count)) %>% pull(count),
          E %>% filter(from != !!from, to == !!to) %>% summarise(count = sum(count)) %>% pull(count),
          E %>% filter(from == !!from, to != !!to) %>% summarise(count = sum(count)) %>% pull(count),
          E %>% filter(from != !!from, to != !!to) %>% summarise(count = sum(count)) %>% pull(count)
          ),
        nrow = 2,
        dimnames = list(From = c(from, paste0("Not_", from)),
                        To = c(to, paste0("Not_", to))
        )
      )
  
    # Broom + FT
    ft = fisher.test(CTb, alternative = 'greater')
    ft = broom::tidy(ft)
  
    ft$from = from
    ft$to = to
    
    ft$POS_POS = CTb[1,1]
    ft$POS_NEG = CTb[1,2]
    ft$NEG_POS = CTb[2,1]
    ft$NEG_NEG = CTb[2,2]
  
    ft
  }
  
  # All tests for all entries of E
  tests = apply(
    data.frame(E, stringsAsFactors = FALSE), 
    1, 
    function(w) single_test(w['from'], w['to'])
    )
  
  tests = Reduce(bind_rows, tests)
  
  # Number of tests by target gene is used for MHT
  N_by_gene = tests %>% 
    group_by(to) %>%
    summarise(N = n())
  
  N_by_gene_v = N_by_gene %>% pull(N)
  names(N_by_gene_v) = N_by_gene %>% pull(to)
  
  # MHT
  tests = tests %>%
    mutate(
      alpha_level = !!alpha_level,
      N = N_by_gene_v[to],
      psign = p.value < alpha_level/N
    ) %>%
    arrange(p.value)
}

# plot_enrichment_test = function(x, alpha_level = 0.05)
# {
#   tests = enrichment_test_incoming_edge(
#     E = x$fit$penalty  %>%
#       filter(count > 3), 
#     alpha_level)
#   
#   ggplot(tests %>%
#            group_by(to) %>%
#            filter(any(psign)),
#          aes(x = POS_POS, y = p.value, color = to, fill = to, shape =  psign)
#          ) +
#     geom_point(size = 1.5) +
#     # scale_color_gradient(low = 'steelblue', high = 'darkred', breaks = 0.05) +
#     # scale_color_manual(values = c(`TRUE` = 'orange', `FALSE` = 'black')) +
#     scale_y_log10() +
#     # scale_x_log10() +
#     coord_cartesian(clip = 'off') +
#     geom_label_repel(
#       data = tests %>% filter(psign),
#       aes(label = paste0(from, ' -> ', to), fill = NA),
#       size = 2,
#       force = 2,
#       nudge_x = 0.3,
#       nudge_y = 0.3
#     ) +
#     theme_minimal(base_size = 10 * cex) +
#     theme(
#       legend.position = 'bottom',
#       legend.key.size = unit(3, 'mm')
#     ) +
#     # geom_hline(yintercept = alpha/nrow(tests), size = .5, color = 'red', linetype = 'dashed') +
#     labs(
#       title = paste0("P-values from the association test"),
#       y = "Node different from A as parent (counts)",
#       x = "A as parent (counts)",
#       caption = paste0("One-sided Fisher test")
#     ) 
#     # facet_wrap(~to)
# 
# }
