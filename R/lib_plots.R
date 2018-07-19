#' @importFrom graphics layout par
plot_fitStats = function(x, palette = 'Accent', file = NA)
{
  IDs = sort(x$fit$solutionID)
  NumSols = sapply(names(IDs), function(w){ length(fit$phylogenies[[w]]) })
  InfTransf = sapply(names(IDs), rev_count_information_transfer_comb, x = x)

  ord = sort(IDs + NumSols + InfTransf)
  IDs = IDs[names(ord)]
  NumSols = NumSols[names(ord)]
  InfTransf = InfTransf[names(ord)]

  layout(matrix(c(1,2,3,1,2,3,4,5,6), byrow = T, ncol = 3))

  barplot(
    IDs,
    col = sapply(IDs, function(w){
      if(w > 1) 'darkred'
      else 'darkgray'
      }),
    border = NA,
    pch = 19,
    xlab = NA,
    ylab = NA,
    main = ('Solution ID'),
    horiz = T,
    las = 2,
    cex.names = .5
  )
  legend('bottomright', legend = c('top-rank', 'non top-rank'), col = c('darkgray', 'darkred'), pch = 19, bty = 'n')

  barplot(
    # log(NumSols),
    NumSols,
    col = sapply(log(NumSols), function(w){
      if(w > 0) 'steelblue'
      else 'darkgray'
    }),
    border = NA,
    pch = 19,
    xlab = NA,
    ylab = NA,
    main = '# Solutions',
    horiz = T,
    las = 2,
    log = 'x',
    cex.names = .5
  )
  legend('bottomright', legend = c('> 1', '= 1'), col = c('steelblue', 'darkgray'), pch = 19, bty = 'n')


  # abline(v = 1, col = 'red', lty = 2)

  barplot(
    # log(InfTransf),
    InfTransf,
    col = sapply(log(InfTransf), function(w){
      if(w > 0) 'forestgreen'
      else 'darkgray'
    }),
    border = NA,
    pch = 19,
    xlab = NA,
    ylab = NA,
    main = ('# Inf. Tr.'),
    horiz = T,
    las = 2,
    log = 'x',
    cex.names = .5
  )
  legend('bottomright', legend = c('= 1', '> 1'), col = c('darkgray', 'forestgreen'), pch = 19, bty = 'n')


  hist(IDs, breaks = 22, col = 'black', border = NA, xlab = NA, main = '')
  hist(NumSols, breaks = 22, col = 'black', border = NA, xlab = NA, main = '')
  hist(InfTransf, breaks = 22, col = 'black', border = NA, xlab = NA, main = '')


  par(mfrow = c(1,1))

  if(!is.na(file)) dev.copy2pdf(file = file)

}
