resample.test.greater <- function(x, y, nBoots=10000){
    n.x = length(x)
    n.y = length(y)
    x.bmean = rep(NA, nBoots)
    y.bmean = rep(NA, nBoots)
    t = rep(NA, nBoots)
    for (i in 1:nBoots){
        x.bmean[i] = mean(sample(x, n.x, replace=TRUE))
        y.bmean[i] = mean(sample(y, n.y, replace=TRUE))
        t[i] = x.bmean[i] - y.bmean[i]
    }
    p = (sum(t < 0)+1)/(length(t)+1)
    return(p)
}

resample.test.pooled <- function(x, y, nBoots=10000){
    mean.diff = mean(x) - mean(y)
    z = c(x,y)
    n.x = length(x)
    n.y = length(y)
    labs = c(rep('x',n.x), rep('y', n.y))
    t = rep(NA, nBoots)
    for (i in 1:nBoots){
        labs.i = sample(labs)
        x.mean = mean(z[which(labs.i=='x')])
        y.mean = mean(z[which(labs.i=='y')])
        t[i] = x.mean - y.mean
    }
    p = (sum(t > mean.diff)+1)/(length(t)+1)
    return(p)
}


generate.boot.nonparametric <- function(variants, cluster.col.name='cluster',
                          vaf.col.names=NULL, vaf.in.percent=TRUE,
                          num.boots=1000,
                          zero.sample=NULL){
    cat('Generating boostrap samples...\n')
    boot.means = NULL
    if (is.null(vaf.col.names)){
        vaf.col.names = setdiff(colnames(variants), cluster.col.name)
    }

    # if no cluster or no sample provided, return NULL
    clusters = unique(variants[[cluster.col.name]])
    num.clusters = length(clusters)
    num.samples = length(vaf.col.names)
    if (num.samples == 0 || num.clusters == 0){return(NULL)}
    if (vaf.in.percent){
        variants[,vaf.col.names] = variants[,vaf.col.names]/100.00
    }
    # make separate data frame for each cluster
    v = list()
    for (cl in clusters){
        v1 = variants[variants[[cluster.col.name]]==cl, vaf.col.names]
        v[[as.character(cl)]] = v1
    }

    # debug
    # vv <<- v


    # generate bootstrap samples for each cluster, each sample
    num.variants.per.cluster = table(variants[[cluster.col.name]])
    #print(num.variants.per.cluster)

    boot.means = list()
    clusters = as.character(clusters)
    zeros = c()
    for (vaf.col.name in vaf.col.names){
        #cat('Booting sample: ', vaf.col.name, '\n')
        sample.boot.means = matrix(NA, nrow=num.boots, ncol=num.clusters)
        colnames(sample.boot.means) = clusters
        rownames(sample.boot.means) = seq(1, num.boots)
        for (cl in clusters){
            # learn zero samples from data, if a cluster has median VAF = 0, consider
            # it as a sample generated from true VAF = 0
            vafs = v[[cl]][[vaf.col.name]]
            if (median(vafs)==0){zeros = c(zeros, vafs)}

            boot.size = num.variants.per.cluster[cl]
            #cat('Booting cluster: ', cl, 'boot.size=', boot.size, '\n')
            for (b in 1:num.boots){
                s.mean = mean(sample(v[[cl]][[vaf.col.name]], boot.size,
                                     replace=TRUE))
                sample.boot.means[b, cl] = s.mean
            }
        }
        boot.means[[vaf.col.name]] = sample.boot.means
    }

    #generate bootstrap means for zero sample
    if (is.null(zero.sample)){
        zero.sample = zeros
    }
    if (length(zero.sample) > 0){
        zero.sample.boot.means = rep(NA, num.boots)
        zero.sample.boot.size = length(zero.sample)
        for (b in 1:num.boots){
            s.mean = mean(sample(zero.sample, zero.sample.boot.size,
                                 replace=TRUE))
            zero.sample.boot.means[b] = s.mean
        }
        boot.means$zero.means = zero.sample.boot.means
    }
    return(boot.means)
}

boot.vaf.ci <- function(boot, sample, cluster, alpha=0.05){
    vaf.means = boot[[sample]][,as.character(cluster)]
    mean.vaf = mean(vaf.means)
    upper.vaf = quantile(vaf.means, 1-alpha/2)
    lower.vaf = quantile(vaf.means, alpha/2)
    upper.vaf.fmt = sprintf('%0.2f%%', upper.vaf)
    lower.vaf.fmt = sprintf('%0.2f', lower.vaf)
    vaf.ci.str = paste0(lower.vaf.fmt, ' - ', upper.vaf.fmt)
    return(list(vaf.mean=mean.vaf,
                vaf.lower=lower.vaf, vaf.upper=upper.vaf,
                vaf.ci.str=vaf.ci.str))
}


subclonal.test <- function(vaf.col.name, parent.cluster, sub.clusters=NULL,
                           boot=NULL, cdf=NULL, min.cluster.vaf=0, alpha=0.05,
                           alternative='greater'){
    # debug
    # cat('subclonal.test: sample=', vaf.col.name, 'parent.cluster=', parent.cluster,
    #    'sub.clusters=', paste(sub.clusters, collapse=','),'\n')

    # if min.cluster.vaf provided,
    #zero.vaf = ifelse(is.null(sub.clusters) & !is.null(min.cluster.vaf),
    #    min.cluster.vaf, 0)
    if (is.null(boot) || length(boot) == 0){
        #cat('NO BOOT\n')
        # test using absolute values of VAF
        # debug
        # cat('No bootstrap! absolute value comparison.\n')
        min.cluster.vaf = ifelse(is.null(min.cluster.vaf), 0, min.cluster.vaf)
        zero.vaf = min.cluster.vaf
        if (is.null(sub.clusters)){
            mean.free.vaf = cdf$vaf[cdf$lab==parent.cluster]
            p = ifelse(mean.free.vaf > zero.vaf, 1, 0)
        }else{
            mean.free.vaf = (cdf$vaf[cdf$lab==parent.cluster] -
                sum(cdf$vaf[cdf$lab %in% sub.clusters]))
            p = ifelse(mean.free.vaf >= 0, 1, 0)
        }
        #debug
        #if (parent.cluster == 'c1' && (all(c('c1a', 'c1aC1') %in% sub.clusters))){
        #  cat('mean.free.vaf=', mean.free.vaf, 'p=', p, '\n')
        #stop()
        #}
        #zz <<- mean.free.vaf
        free.vaf.mean=mean.free.vaf
        upper.free.vaf = NA
        lower.free.vaf = NA
        upper.free.vaf.fmt = NA
        lower.free.vaf.fmt = NA
        confident.level = NA
        confident.level.non.negative = NA
        free.vaf.ci.str = NA

    }else{
        num.boots = nrow(boot[[1]])
        if (num.boots == 0){return(NULL)}
        if (is.null(sub.clusters)){
            #free.vaf = boot[[vaf.col.name]][,parent.cluster]
            # -boot$zero.means
            free.vaf = boot[[vaf.col.name]][,parent.cluster]
            #cat('debug: AAA\n')
        }else{

            free.vaf = apply(boot[[vaf.col.name]], 1,
                             function(row) (row[parent.cluster] -
                                                sum(row[sub.clusters])))
            #cat('debug: BBB\n')
            #bbb <<- boot
            #cat(vaf.col.name, parent.cluster, '\n')
            #print(sub.clusters)
        }
        zz <<- free.vaf
        zero.vaf = 0
        # p = probability that clone has non-negative ccf
        # also equal p-value of test to reject Ho: ccf < 0
        # TODO: change to >= in free.vaf > zero.vaf???
        p = sum(free.vaf > zero.vaf)/length(free.vaf)
        mean.free.vaf = mean(free.vaf)
        upper.free.vaf = quantile(free.vaf, 1-alpha/2)
        lower.free.vaf = quantile(free.vaf, alpha/2)
        upper.free.vaf.fmt = sprintf('%0.2f%%', upper.free.vaf)
        lower.free.vaf.fmt = sprintf('%0.2f', lower.free.vaf)
        confident.level = 1 - alpha
        confident.level.non.negative = ifelse(lower.free.vaf >= 0,
            confident.level, sum(free.vaf >= 0 &
            free.vaf <= upper.free.vaf)/length(free.vaf)
        )
        free.vaf.ci.str = ifelse(lower.free.vaf > 0,
                                 paste0(lower.free.vaf.fmt, ' - ',
                                        upper.free.vaf.fmt),
                                 paste0('0 - ', upper.free.vaf.fmt))
    }
    #debug
    # cat('p-value =', p, '\n')
    #cat('CI =', lower.free.vaf, '-', upper.free.vaf, 'free.mean=', mean.free.vaf, '\n')

    # default was less in the past, now just need to recalc p as 1-p
    # so new p is prob that ccf is negative
    # or the p.value to reject Ho: ccf < 0
    if (alternative == 'greater'){p = 1 - p}

    return(list(free.vaf.ci=free.vaf.ci.str,
                free.vaf.mean=mean.free.vaf,
                free.vaf.lower=lower.free.vaf,
                free.vaf.upper=upper.free.vaf,
                free.vaf.confident.level=confident.level,
                free.vaf.confident.level.non.negative=confident.level.non.negative,
                p.value=p))
}


fisher.p <- function(pvals, max.p=1){
    return(pchisq( -2*sum(log(pvals)), df=length(pvals), lower.tail=FALSE))
}

 combine.p <- function(pvals, method='fisher', max.p=1, w=NULL){
    if (method == 'z'){
        library(metap)
    }
    if (is.null(w)){
        w = rep(1, length(pvals))
    }
    keep = !is.na(pvals) & pvals <= max.p
    pvals = pvals[keep]
    w = w[keep]
    n = length(pvals)
    if (n < 1){
        return(NA)
    } else if (n == 1){
        return(pvals[1])
    }else if (method == 'fisher'){
        return(fisher.p(pvals, max.p))
    }else if (method == 'z'){
        return(sumz(pvals, weights=w)$p)
    }
}

