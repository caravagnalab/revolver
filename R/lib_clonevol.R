###########################################################################################################
# I downloaded the ClonEvol repo and modified slighlty the main running function to skip the consensus
# step at the end of the algorithm. We want to return just the trees per region. Afterwards, I removed all
# functions that are not used WITH THE CURRENT execution setting for ClonEvol
#
# GIULIO CARAVAGNA
###########################################################################################################

###############################################################################
# Clonevol: Inferring and visualizing clonal evolution in multi-sample cancer
# sequencing
# Created by: Ha Dang <haxdangATgmailDOTcom>
# Date: Dec. 25, 2014
#
# Purposes:
#  Infer and visualize clonal evolution in multi cancer samples
#  using somatic mutation clusters and their variant allele frequencies
#
###############################################################################

make.clonal.data.frame <- function (vafs, labels, add.normal=FALSE,
                                    #normal.clone.color='#f0f0f0', # very light gray
                                    normal.clone.color='#e5f5f9', # very light blue
                                    founding.label=NULL, colors=NULL){
    v = data.frame(lab=as.character(labels), vaf=vafs, stringsAsFactors=FALSE)
    if (is.null(colors)){
        #colors = c('#a6cee3', '#b2df8a', '#cab2d6', '#fdbf6f', '#fb9a99',
        #           '#d9d9d9','#999999', '#33a02c', '#ff7f00', '#1f78b4',
        #           '#fca27e', '#ffffb3', '#fccde5', '#fb8072', '#b3de69',
        #           'f0ecd7', rep('#e5f5f9',1))
        # colors = get.clonevol.colors(nrow(v))
        colors = NULL
        # if normal clone added, set it color
        if (v$lab[1] == '0'){colors = c(normal.clone.color, colors)}
    }
    clone.colors = colors[seq(1,nrow(v))]
    v$color = clone.colors
    v = v[order(v$vaf, decreasing=TRUE),]
    if (!is.null(founding.label)){
        #make founding clone first in data frame
        v1 = v[v$lab == founding.label,]
        v2 = v[v$lab != founding.label,]
        v = rbind(v1, v2)
    }
    # add dummy normal cluster to cover
    if (add.normal){
        v = rbind(data.frame(vaf=0.5, lab='0', color=colors[length(colors)],
                             stringsAsFactors=FALSE), v)
    }

    v$parent = NA
    v$ancestors = '-'
    v$occupied = 0
    v$free = v$vaf
    v$free.mean = NA
    v$free.lower = NA
    v$free.upper = NA
    v$free.confident.level = NA
    v$free.confident.level.non.negative = NA
    v$p.value = NA
    v$num.subclones = 0
    v$excluded = FALSE
    #rownames(v) = seq(1,nrow(v))
    rownames(v) = v$lab
    #print(str(v))
    #print(v)
    return(v)
}

is.ancestor <- function(v, a, b){
    #cat('Checking if', a, '---ancestor-->', b, '\n')
    if (is.na(b) || b == '-1'){
        return(FALSE)
    }else{
        par = v$parent[v$lab == b]
        if(is.na(par)){
            return(FALSE)
        }else if (par == a){
            return(TRUE)
        }else{
            return(is.ancestor(v, a, par))
        }
    }
}


estimate.ccf <- function(vx, sample, i, boot, min.cluster.vaf,
    alpha, t=NULL, sub.clusters=NULL){
    if (is.null(t)){
        t = subclonal.test(sample,
           as.character(vx[i,]$lab),
           sub.clusters=sub.clusters, boot=boot,
           cdf=vx,
           min.cluster.vaf=min.cluster.vaf,
           alpha=alpha)
    }
    vx$free.mean[i] = t$free.vaf.mean
    vx$free.lower[i] = t$free.vaf.lower
    vx$free.upper[i] = t$free.vaf.upper
    vx$p.value[i] = t$p.value
    vx$free.confident.level[i] =
        t$free.vaf.confident.level
    vx$free.confident.level.non.negative[i] =
        t$free.vaf.confident.level.non.negative
    return(vx)

}


enumerate.clones <- function(v, sample=NULL, variants=NULL,
                             founding.cluster = NULL,
                             ignore.clusters=NULL,
                             subclonal.test.method='bootstrap',
                             boot=NULL,
                             p.value.cutoff=0.05,
                             alpha=0.05,
                             min.cluster.vaf=0){
    cat(sample, ': Enumerating clonal architectures...\n')
    vv = list() # to hold list of output clonal models
    #cat('*********: p : ', p.value.cutoff, '\n')
    findParent <- function(v, i){
        #print(i)
        if (i > nrow(v)){
            #debug
            #print(v)

            v$is.zero = ifelse(v$free.lower >= 0, FALSE, TRUE)
            # determine subclone
            clone.stat = determine.subclone(v, v$lab[!is.na(v$parent)
                                            & v$parent == '-1'])

            v$is.subclone = clone.stat$is.sub[v$lab]
            v$is.founder = clone.stat$is.founder[v$lab]
            rownames(v) = v$lab
            vv <<- c(vv, list(v))
        }else{
            #print(head(v))
            vaf = v[i,]$vaf
            if (!is.na(v[i,]$parent) && v[i,]$parent == '-1'){# root
                vx = v
                # estimate CCF for root if it does not have subclones nested yet,
                # just in case there is no other clone to be nested
                if (vx$num.subclones[i] == 0){
                    vx = estimate.ccf(vx, sample, i, boot, min.cluster.vaf, alpha, t=NULL)
                }
                findParent(vx, i+1)
            }else if (v[i,]$excluded){
                vx = v
                vx$parent[i] = NA
                findParent(vx, i+1)
            }else{
                #for (j in 1:(i-1)){
                for (j in 1:nrow(v)){
                    parent.cluster = as.character(v[j,]$lab)
                    current.cluster = as.character(v[i,]$lab)
                    is.ancestor = is.ancestor(v, current.cluster,
                                              parent.cluster)
                    #print(v)
                    #cat('i=', i, 'j=', j, '\n')
                    if (i != j && !v$excluded[j] && !is.ancestor){
                        # assign cluster in row j as parent of cluster in row i
                        sub.clusters = as.character(c(v$lab[!is.na(v$parent) &
                                                    v$parent == parent.cluster],
                                                    current.cluster))
                        #print(str(v))
                        # debug
                        # cat('Testing...', sample, '-', j, parent.cluster,
                          # 'sub clusters:', sub.clusters, '\n')
                        t = subclonal.test(sample, parent.cluster, sub.clusters,
                                           boot=boot,
                                           cdf=v,
                                           min.cluster.vaf=min.cluster.vaf,
                                           alpha=alpha)
                        # hdng: test direction changed to greater, so p = 1 - p, sign flipped to <
                        # if a clonal nesting do not violate sum rule, this is unresolvable
                        # so it will be recorded as a temporary solution, later, other samples
                        # come in, we may find one or a few models resolvable, then applying
                        # cross rule (matching between samples) will solve the model
                        if(t$p.value < 1 - p.value.cutoff){
                            vx = v
                            # debug
                            #cat(i, '<-', j, 'vaf=', vaf, '\n')
                            #print(head(v))
                            #print(head(vx))
                            vx$p.value[j] = t$p.value
                            vx$free.mean[j] = t$free.vaf.mean
                            vx$free.lower[j] = t$free.vaf.lower
                            vx$free.upper[j] = t$free.vaf.upper
                            vx$free.confident.level[j] =
                                t$free.vaf.confident.level
                            vx$free.confident.level.non.negative[j] =
                                t$free.vaf.confident.level.non.negative
                            vx$free[j] = vx$free[j] - vaf
                            vx$occupied[j] = vx$occupied[j] + vaf
                            vx$num.subclones[j] = length(sub.clusters)
                            #vx$parent[i] = vx[j,]$lab
                            vx$parent[i] = parent.cluster
                            vx$ancestors[i] = paste0(vx$ancestors[j],
                                paste0('#',parent.cluster,'#'))

                            # calculate confidence interval for vaf estimate of
                            # the subclone if it does not contain other
                            # subclones (will be overwrite later
                            # if subclones are added to this subclone)
                            #if (is.na(vx$free.lower[i])){
                            if (vx$num.subclones[i] == 0){
                                #t = subclonal.test(sample,
                                #       as.character(vx[i,]$lab),
                                #       sub.clusters=NULL, boot=boot,
                                #       cdf=vx,
                                #       min.cluster.vaf=min.cluster.vaf,
                                #       alpha=alpha)
                                #vx$free.mean[i] = t$free.vaf.mean
                                #vx$free.lower[i] = t$free.vaf.lower
                                #vx$free.upper[i] = t$free.vaf.upper
                                #vx$p.value[i] = t$p.value
                                #vx$free.confident.level[i] =
                                #    t$free.vaf.confident.level
                                #vx$free.confident.level.non.negative[i] =
                                #    t$free.vaf.confident.level.non.negative
                                vx = estimate.ccf(vx, sample, i, boot,
                                            min.cluster.vaf, alpha, t=NULL)
                            }
                            findParent(vx, i+1)
                        }
                    }
                }
            }
        }
    }

    # exclude some cluster with VAF not significantly diff. from zero
    # print(v)
    cat('Determining if cluster VAF is significantly positive...\n')
    if (is.null(min.cluster.vaf)){
        cat('No min.cluster.vaf provided. Using bootstrap\n')
    }else{
        cat('Exluding clusters whose VAF < min.cluster.vaf=',
                min.cluster.vaf, '\n', sep='')
    }
    for (i in 1:nrow(v)){
        cl = as.character(v[i,]$lab)
        if (is.null(min.cluster.vaf)){
            # test if VAF of this cluster cl is > 0
            t = subclonal.test(sample, parent.cluster=cl, sub.clusters=NULL,
                           boot=boot, min.cluster.vaf=min.cluster.vaf,
                           alpha=alpha)
            # if not, exclude from analysis
            v[i,]$excluded = ifelse(t$p.value > p.value.cutoff, TRUE, FALSE)
        }else{
            # if the median/mean (estimated earlier) VAF < e, do
            # not consider this cluster in this sample
            v[i,]$excluded = ifelse(v[i,]$vaf < min.cluster.vaf, TRUE, FALSE)
        }

    }

    cat('Non-positive VAF clusters:',
        paste(v$lab[v$excluded], collapse=','), '\n')

    # also exlude clusters in the ignore.clusters list
    if (!is.null(ignore.clusters)){
        ignore.idx = v$lab %in% as.character(ignore.clusters)
        v$excluded[ignore.idx] = TRUE
        cat('User ignored clusters: ',
            paste(v$lab[ignore.idx], collapse=','), '\n')
    }
    #print(v)

    # if normal sample (0) is included, the normal sample
    # will be root (polyclonal model), otherwise find the
    # founding clone and place it first
    if (v[1,]$lab == 0 || v[1,]$lab == '0'){
        v[1,]$parent = -1
        findParent(v, 2)
    }else{
        #print(founding.cluster)
        if (is.null(founding.cluster)){
            max.vaf = max(v$vaf)
            roots = rownames(v)[v$vaf == max.vaf]
        }else{
            roots = rownames(v)[v$lab == founding.cluster]
        }
        # debug
        #cat('roots:', paste(roots, collapse=','), '\n')
        for (r in roots){
            #print(roots)
            vr = v
            vr[r,]$parent = -1
            #print(vr)
            findParent(vr,1)
        }
    }

    return(vv)
}




# REQUIRED
determine.subclone <- function(v, r){
    rownames(v) = v$lab
    next.clones = c(r)
    is.sub = rep(NA, nrow(v))
    names(is.sub) = v$lab
    is.founder = rep(NA, nrow(v))
    names(is.founder) = v$lab

    v$is.zero = ifelse(v$free.lower >= 0, FALSE, TRUE)

    # if no confidence interval estimated (no bootstrap model)
    if (all(is.na(v$free.lower))){
        v$is.zero = ifelse(v$free.mean > 0, TRUE, FALSE)
    }

    while (length(next.clones) > 0){
        cl = next.clones[1]
        children = v$lab[!is.na(v$parent) & v$parent == cl]
        next.clones = c(next.clones[-1], children)
        par = v[cl, 'parent'];
        if (!is.na(par) && (par == '-1' || par=='0')){
            # founding clone in monoclonal model, or clones
            # coming out of normal clone is not subclone
            is.sub[cl] = FALSE
            is.founder[cl] = TRUE
        }
        #if (v[cl, 'free.lower'] <= 0 && v[cl, 'num.subclones'] == 1){
        if (v[cl, 'is.zero'] && v[cl, 'num.subclones'] == 1){
            is.sub[children] = is.sub[cl]
            is.founder[children] = TRUE
        }else{
            is.sub[children] = TRUE
            if(v[cl, 'is.zero']){
                is.founder[children] = TRUE
            }else{
                is.founder[children] = FALSE
            }
        }
    }
    is.founder = is.founder & !v$is.zero
    return(list(is.sub=is.sub, is.founder=is.founder, is.zero=v$is.zero))
}


# REQUIRED
infer.clonal.models <- function(c=NULL, variants=NULL,
                                cluster.col.name='cluster',
                                founding.cluster=NULL,
                                ignore.clusters=NULL,
                                vaf.col.names=NULL,
                                ccf.col.names=NULL,
                                vaf.in.percent=TRUE,
                                depth.col.names=NULL,
                                weighted=FALSE,
                                sample.names=NULL,
                                sample.groups=NULL,
                                model='monoclonal',
                                cancer.initiation.model=NULL,
                                subclonal.test='bootstrap',
                                cluster.center='median',
                                subclonal.test.model='non-parametric',
                                seeding.aware.tree.pruning=FALSE,
                                merge.similar.samples=FALSE,
                                clone.colors=NULL,
                                random.seed=NULL,
                                boot=NULL,
                                num.boots=1000,
                                p.value.cutoff=NULL,
                                sum.p.cutoff=0.01,
                                cross.p.cutoff=NULL,
                                alpha=NULL,
                                min.cluster.vaf=0.01,
                                score.model.by='probability',
                                verbose=TRUE){
    # backward compatible with old p.value.cutoff
    if (!is.null(p.value.cutoff)){sum.p.cutoff = p.value.cutoff}
    if (is.null(alpha)){alpha = sum.p.cutoff}
    #if (is.null(cross.p.cutoff)){cross.p.cutoff = sum.p.cutoff}

    if (is.null(vaf.col.names) && is.null(ccf.col.names)){
        # check format of input, find vaf column names
        if(!is.null(c)){
            cnames = colnames(c)
        }else if(!is.null(variants)){
            cnames = colnames(variants)
        }else{
            stop('ERROR: Need at least parameter c or variants\n')
        }
        if (!(cluster.col.name %in% cnames && length(cnames) >= 2)){
            stop('ERROR: No cluster column and/or no sample\n')
        }
        vaf.col.names = setdiff(cnames, cluster.col.name)
    }

    if (!is.null(vaf.col.names) && !is.null(ccf.col.names)){
        cat('WARN: Both VAF and CCF provided. Will analyze using CCF.\n')
    }

    if (!is.null(ccf.col.names)){
        # calculate VAF as half of CCF, and use ccf.col.names
        # for vaf.col.names
        cat('Calculate VAF as CCF/2\n')
        variants[,ccf.col.names] = variants[,ccf.col.names]/2
        vaf.col.names = ccf.col.names
    }

    # convert cluster column to character (allow both number and string as
    # cluster or clone IDs)
    founding.cluster = as.character(founding.cluster)
    if (!is.null(c)) {
        c[[cluster.col.name]] = as.character(c[[cluster.col.name]])
    }
    if (!is.null(variants)){
        variants[[cluster.col.name]] = as.character(variants[[cluster.col.name]])
    }
    if (is.null(c) && is.null(variants)){
        stop('ERROR: No variant clustering result provided via c or variants parameters.\n')
    }

    # backward compatible between params: cancer.initiation.model and  model
    if (!is.null(cancer.initiation.model)){model = cancer.initiation.model}


    if (is.null(sample.names)){
        sample.names = vaf.col.names
    }

    nSamples = length(sample.names)
    n.vaf.cols = length(vaf.col.names)

    if (nSamples != n.vaf.cols || nSamples == 0){
        stop('ERROR: sample.names and vaf.col.names have different length
         or both have zero length!\n')
    }

    if (nSamples >= 1 && verbose){
        for (i in 1:nSamples){
            cat('Sample ', i, ': ', sample.names[i], ' <-- ',
                vaf.col.names[i], '\n', sep='')
        }
    }

    # if polyclonal model, add normal as founding clone
    add.normal = NA
    if (model == 'monoclonal'){
        add.normal = FALSE
    }else if (model == 'polyclonal'){
        add.normal = TRUE
        founding.cluster = '0'
        # add a faked normal clone with VAF = norm(mean=50, std=10)
        # TODO: follow the distribution chosen by user
        if (add.normal){
            tmp = variants[rep(1,100),]
            tmp[[cluster.col.name]] = founding.cluster
            vaf50 = matrix(rnorm(100, 50, 10), ncol=1)[,rep(1, length(vaf.col.names))]
            tmp[, vaf.col.names] = vaf50
            variants = rbind(tmp, variants)
        }

    }
    if (is.na(add.normal)){
        stop(paste0('ERROR: Model ', model, ' not supported!\n'))
    }
    if(verbose){cat('Using ', model, ' model\n', sep='')}

    # prepare cluster data and infer clonal models for individual samples
    if (is.null(c)){
        c = estimate.clone.vaf(variants, cluster.col.name,
                               vaf.col.names, vaf.in.percent=vaf.in.percent,
                               method=cluster.center)
    }
    vv = list()
    for (i in 1:nSamples){
        s = vaf.col.names[i]
        sample.name = sample.names[i]
        v = make.clonal.data.frame(c[[s]], c[[cluster.col.name]],
            colors=clone.colors)
        if (subclonal.test == 'none'){
            #models = enumerate.clones.absolute(v)
            models = enumerate.clones(v, sample=s,
                                      founding.cluster=founding.cluster,
                                      min.cluster.vaf=min.cluster.vaf,
                                      ignore.clusters=ignore.clusters)
        }
        # else if (subclonal.test == 'bootstrap'){
        #     # # if (is.null(boot)){
        #     #     #boot = generate.boot(variants, vaf.col.names=vaf.col.names,
        #     #     #                     vaf.in.percent=vaf.in.percent,
        #     #     #                     num.boots=num.boots)
        #     #
        #     #     boot = generate.boot(variants, vaf.col.names=vaf.col.names,
        #     #                          depth.col.names=depth.col.names,
        #     #                          vaf.in.percent=vaf.in.percent,
        #     #                          num.boots=num.boots,
        #     #                          bootstrap.model=subclonal.test.model,
        #     #                          cluster.center.method=cluster.center,
        #     #                          weighted=weighted,
        #     #                          random.seed=random.seed)
        #     #     #bbb <<- boot
        #     # }
        #
        #     models = enumerate.clones(v, sample=s, variants, boot=boot,
        #                               founding.cluster=founding.cluster,
        #                               ignore.clusters=ignore.clusters,
        #                               min.cluster.vaf=min.cluster.vaf,
        #                               p.value.cutoff=sum.p.cutoff,
        #                               alpha=alpha)
        # }

        if(verbose){cat(s, ':', length(models),
                        'clonal architecture model(s) found\n\n')}
        if (length(models) == 0){
            print(v)
            message(paste('ERROR: No clonal models for sample:', s,
                       '\nCheck data or remove this sample, then re-run.
                       \nAlso check if founding.cluster was set correctly!'))
            return(NULL)
        }else{
            vv[[sample.name]] = models
        }
    }

    # infer clonal evolution models,  an accepted model must satisfy
    # the clone-subclonal relationship
    matched = NULL
    scores = NULL
    probs = NULL
    if (nSamples == 1 && length(vv[[1]]) > 0){
        num.models = length(vv[[1]])
        matched = data.frame(x=1:num.models)
        colnames(matched) = c(sample.names[1])
        scores = data.frame(x=rep(0,num.models))
        probs = data.frame(x=rep(0,num.models))
        colnames(scores) = c(sample.names[1])
        colnames(probs) = c(sample.names[1])
        merged.trees = list()
        merged.traces = NULL
        for (i in 1:num.models){
            scores[i,1] = get.model.score(vv[[1]][[i]])
            probs[i,1] = get.model.non.negative.ccf.prob(vv[[1]][[i]])
            merged.trees = c(merged.trees, list(vv[[1]][[i]]))
        }
        #scores$model.score = scores[, 1]
        probs$model.prob = probs[,1]
    }
    # if (nSamples >= 2){
    #
    #
    #     z = find.matched.models(vv, sample.names, sample.groups,
    #         merge.similar.samples=merge.similar.samples)
    #
    #     matched = z$matched.models
    #     scores = z$scores
    #     merged.trees = z$merged.trees
    #     merged.traces = z$merged.traces
    #     vv = z$models
    #     if (!is.null(matched)){
    #         rownames(matched) = seq(1,nrow(matched))
    #         colnames(matched) = sample.names
    #         matched = as.data.frame(matched)
    #         rownames(scores) = seq(1,nrow(matched))
    #         colnames(scores) = sample.names
    #         scores = as.data.frame(scores)
    #         # TODO: this prod is supposed to be prod of prob, but with 2nd scoring
    #         # strategy using max.p and pval combination, this is prod of max.p of
    #         # individual samples
    #         # scores$model.score = apply(scores, 1, prod)
    #
    #         probs = matrix(nrow=nrow(matched), ncol=(ncol(matched)))
    #         colnames(probs) = colnames(matched)
    #         for (model.idx in 1:nrow(matched)){
    #             for(samp in colnames(matched)){
    #                 probs[model.idx,samp] = get.model.non.negative.ccf.prob(
    #                     vv[[samp]][[matched[model.idx,samp]]])
    #             }
    #         }
    #         probs = as.data.frame.matrix(probs)
    #         probs$model.prob = apply(probs[, sample.names], 1, prod)
    #     }
    # }
    if (!is.null(matched)) {
        # sort models by score
        #idx = 1:nrow(scores)
        #if (score.model.by == 'metap'){
        #    # cross sample p-value
        #    print(scores); print(str(scores))
        #    idx = order(scores$max.clone.ccf.combined.p, deceasing=T)
        #}else if (score.model.by == 'probability'){
            # non-neg ccf probability
            idx = order(probs$model.prob, decreasing=TRUE)
        #}
        matched = matched[idx, ,drop=FALSE]
        scores = scores[idx, , drop=FALSE]
        probs = probs[idx, , drop=FALSE]
        merged.trees = merged.trees[idx]
        merged.traces = merged.traces[idx]
    }
    # num.matched.models = ifelse(is.null(matched), 0, nrow(matched))
    # if (verbose){ cat(paste0('Found ', num.matched.models,
    #                          ' consensus model(s)\n'))}
    # # trim and remove redundant merged.trees
    # cat('Pruning consensus clonal evolution trees....\n')
    # trimmed.merged.trees = trim.clone.trees(merged.trees,
    #     seeding.aware.tree.pruning=seeding.aware.tree.pruning)$unique.trees
    # cat('Number of unique pruned consensus trees:', length(trimmed.merged.trees), '\n')
    merged.trees = merged.traces = trimmed.merged.trees = num.matched.models = NULL
    results = list(models=vv, matched=list(index=matched,
        merged.trees=merged.trees, merged.traces=merged.traces,
        scores=scores, probs=probs, trimmed.merged.trees=trimmed.merged.trees),
        num.matched.models=num.matched.models)
    # cat('Scoring models...\n')
    # results = cross.rule.score(results, rank=(score.model.by=='metap'))
    # if (!is.null(cross.p.cutoff)){
    #     num.sig.models = sum(results$matched$scores$max.clone.ccf.combined.p <= cross.p.cutoff)
    #     cat(num.sig.models, 'model(s) with p-value <=', cross.p.cutoff, '\n')
    #     if(num.sig.models == 0){
    #         message('\n***WARN: Inra-tumor heterogeneity could result in a clone (eg. founding)
    #          that is not present (no cells) in any samples, although  detectable via
    #          clonal marker variants due to that its subclones are distinct across
    #          samples. Therefore, a model with a higher p-value for the CCF of such
    #          a clone can still be biologically consistent, interpretable, and
    #          interesting! Manual investigation of those higher p-value models
    #          is recommended\n\n')
    #     }
    # }

    # record data and params used
    results$variants = variants
    results$params = list(cluster.col.name=cluster.col.name,
                          sum.p.cutoff=sum.p.cutoff,
                          cross.p.cutoff=cross.p.cutoff,
                          alpha=alpha,
                          min.cluster.vaf=min.cluster.vaf,
                          vaf.col.names=vaf.col.names,
                          vaf.in.percent=vaf.in.percent,
                          sample.groups=sample.groups,
                          num.boots=num.boots,
                          bootstrap.model=subclonal.test.model,
                          cluster.center.method=cluster.center,
                          merge.similar.samples=merge.similar.samples,
                          score.model.by=score.model.by,
                          random.seed=random.seed
                          )
    results$boot=boot

    return(results)
}



# REQUIRED
estimate.clone.vaf <- function(v, cluster.col.name='cluster',
                               vaf.col.names=NULL,
                               vaf.in.percent=TRUE,
                               method='median',
                               ref.count.col.names=NULL,
                               var.count.col.names=NULL,
                               depth.col.names=NULL){
    #clusters = sort(unique(v[[cluster.col.name]]))
    clusters = unique(v[[cluster.col.name]])
    clone.vafs = NULL

    if (is.null(vaf.col.names)){
        vaf.col.names = setdiff(colnames(v), cluster.col.name)
    }

    for (cl in clusters){
        #cat('cluster: ', cl, '\n')
        #print(str(v))
        #print(str(clusters))
        #print(str(cl))
        #print(length(vaf.col.names))
        is.one.sample = length(vaf.col.names) == 1
        if (method == 'median'){
            if (is.one.sample){
                median.vafs = median(v[v[[cluster.col.name]]==cl,vaf.col.names])
                names(median.vafs) = vaf.col.names
            }else{
                median.vafs = apply(v[v[[cluster.col.name]]==cl,vaf.col.names],
                                    2, median)
            }
        }else if (method == 'mean'){
            if (is.one.sample){
                median.vafs = mean(v[v[[cluster.col.name]]==cl,vaf.col.names])
                names(median.vafs) = vaf.col.names
            }else{
                median.vafs = apply(v[v[[cluster.col.name]]==cl,vaf.col.names],
                                    2, mean)
            }
        }
        #print(str(median.vafs))
        median.vafs = as.data.frame(t(median.vafs))
        #print(str(median.vafs))
        #median.vafs[[cluster.col.name]] = cl
        if (is.null(clone.vafs)){
            clone.vafs = median.vafs
        }else{
            clone.vafs = rbind(clone.vafs, median.vafs)
        }
    }
    clone.vafs = cbind(clusters, clone.vafs)
    colnames(clone.vafs)[1] = cluster.col.name
    #clone.vafs = clone.vafs[order(clone.vafs[[cluster.col.name]]),]
    if (vaf.in.percent){
        clone.vafs[,vaf.col.names] = clone.vafs[,vaf.col.names]/100.00
    }
    return(clone.vafs)
}


#' # REQUIRED BUT I FORCED ITS REMOVAL AS IT WAS ALWAYS EXPORTED IN OUR PACKAGE
#' get.clonevol.colors = function(num.colors, strong.color=FALSE) {
#'     colors = c('#cccccc',
#'                #'#b3b3b3',
#'                #'#999999',
#'                '#a6cee3', '#b2df8a', '#cab2d6','#ff99ff', '#fdbf6f', '#fb9a99',
#'                #'#bf812d',
#'                '#bbbb77',
#'                '#cf8d30',
#'                '#41ae76', '#ff7f00',
#'                '#3de4c5',
#'                #'#aaaa55',
#'                '#ff1aff',
#'                #'#c51b7d',
#'                '#9933ff', '#3690c0','#8c510a', '#666633', '#54278f',
#'                '#e6e600',
#'                '#e31a1c', '#00cc00', '#0000cc', '#252525',
#'                #'#ef6548',
#'                '#fccde5',  '#d9d9d9',
#'                #'#33a02c', '#3f007d', '#1f78b4',
#'                '#f0ecd7', '#ffffb3',
#'                #'#ffcc00',
#'                #'#fca27e', '#fb8072',
#'                '#ffff00', rep('#e5f5f9',10000))
#'     if(strong.color){
#'         colors = c('#e41a1c', '#377eb8', '#4daf4a', '#984ea3',
#'             '#ff7f00', 'black', 'darkgray', rep('lightgray',10000))
#'         colors[1:3] = c('red', 'blue', 'green')
#'     }
#'     if (num.colors > length(colors)){
#'         stop('ERROR: Not enough colors!\n')
#'     }else{
#'         return(colors[1:num.colors])
#'     }
#' }



# REQUIRED
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


