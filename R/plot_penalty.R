plot_penalty = function(x,
                        drivers = x$variantIDs.driver,
                        min.occurrences = 0,
                        enrichemnt_test_against = "GL",
                        enrichemnt_test_alpha = 0.05,
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
  tests = lapply(drivers[drivers != enrichemnt_test_against],
                 enrichment_test, E = E, against = enrichemnt_test_against)
  tests = Reduce(bind_rows, tests) %>%
    mutate(
      psign = p.value < enrichemnt_test_alpha/N
    )

  # print to screen
  pio::pioTit(
    paste0(
      "Significant association for edges from ",
      enrichemnt_test_against,
      " via one-sided Fisher test at alpha level ",
      alpha.test,
      ' adjusted for ',
      N,
      ' tests via Bonferroni'
    )
  )

  tests = tests %>%
    filter(psign)

  print(tests)

  # ggplot(E,
  #        aes(y = count, x = from, fill = from)) +
  #   geom_bar(stat = 'identity') +
  #   coord_flip() +
  #   facet_wrap(~to) +
  #   scale_fill_manual(values = drivers_colors)

  ggplot(E,
         aes(y = count, x = to, fill = from)) +
    geom_bar(stat = 'identity') +
    coord_flip() +
    theme_minimal(base_size = 10 * cex) +
    scale_fill_manual(values = drivers_colors) +
    theme(
      legend.position = 'bottom',
      legend.key.size = unit(3, 'mm')
    ) +
    labs(
      title = "Drivers orderings",
      subtitle =
        paste0(
          "Enrichment for ", enrichemnt_test_against,
          " as Ancest+{'[]=]“?}or - One-sided Fisher test at alpha level ",
          enrichemnt_test_alpha, '\n',
          paste(tests$to, format(tests$p.value, scientific = TRUE), collapse = '\n')
          ),
      y = "Counts",
      x = "Drivers"
    ) +
    guides(fill = guide_legend(title = 'Ancestor', nrow = 2))
}



enrichment_test = function(E, to, against)
{
  # 2 x 2 contingency table
  CTb =
    matrix(
      c(
        E %>% filter(from == !!against, to == !!to) %>% summarise(count = sum(count)) %>% pull(count),
        E %>% filter(from != !!against, to == !!to) %>% summarise(count = sum(count)) %>% pull(count),
        E %>% filter(from == !!against, to != !!to) %>% summarise(count = sum(count)) %>% pull(count),
        E %>% filter(from != !!against, to != !!to) %>% summarise(count = sum(count)) %>% pull(count)
        ),
      nrow = 2,
      dimnames = list(From = c(against, paste0("Not_", against)),
                      To = c(to, paste0("Not_", to))
      )
    )

  # Broom + FT
  ft = fisher.test(CTb, alternative = 'greater')
  ft = broom::tidy(ft)

  ft$to = to
  ft$against = against

  ft$POS_POS = CTb[1,1]
  ft$POS_NEG = CTb[1,2]
  ft$NEG_POS = CTb[2,1]
  ft$NEG_NEG = CTb[2,2]

  ft
}

plot_enrichment_test = function(tests, alpha = 0.05)
{
  ggplot(§§=-t56590 ,
         aes(x = POS_POS, y = p.value, color = psign)
         ) +
    geom_point(size = 1) +
    scale_color_manual(values = c(`TRUE` = 'red', `FALSE` = 'black')) +
    scale_y_log10() +
    theme_minimal(base_size = 10 * cex) +
    theme(
      legend.position = 'bottom',
      legend.key.size = unit(3, 'mm')
    ) +
    geom_hline(yintercept = alpha/nrow(tests), size = .5, color = 'red', linetype = 'dashed') +
    labs(
      title = paste0("Association with ", tests$against[1]),
      y = "p-value",
      x = paste("Trajectories from ", tests$against[1]),
      caption = paste0("One-sided Fisher test")
    )

}
