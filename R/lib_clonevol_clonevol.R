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
        colors=get.clonevol.colors(nrow(v))
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



match.sample.clones <- function(v1, v2){
    compatible = TRUE
    for (i in 1:nrow(v2)){
        vi = v2[i,]
        parent2 = vi$parent
        if (is.na(parent2)){next}
        parent1 = v1[v1$lab == vi$lab,]$parent
        #debug
        #cat(vi$lab, ' par1: ', parent1, 'par2: ', parent2, '\n')
        if (!is.na(parent1) && parent1 != parent2){
            compatible = FALSE
            break
        }
    }
    return(compatible)
}

generate.fill.points <- function(x, y, num.points=50){
    n = length(x)
    k = floor(n/2)
    #z = c(1,2,k,k+1,k+2,k+3,n)
    gen.points <- function(x1, x2, y11, y12, y21, y22, step=NULL){
        xmin = min(x1,x2)
        xmax = max(x1,x2)
        ymid = (y11 + y12)/2
        rx = c()
        ry = c()
        if (is.null(step)){step = (xmax-xmin)/10}
        xx = seq(xmin, xmax, step)
        xx = xx[-length(xx)]
        for (xi in xx){
            yi = y11 + (y21 - y11)/(xmax - xmin)*(xi-xmin)
            yr = runif(num.points, ymid-(yi-ymid), yi)
            xr = runif(num.points, xi, xi+step)
            rx = c(rx, xr)
            ry = c(ry, yr)
        }
        return(list(x=rx, y=ry))
    }

    t1 = gen.points(x[1], x[2], y[1], y[1], y[2], y[n])
    t2 = gen.points(x[2], x[k], y[2], y[k], y[k+2], y[k+3])
    return(list(x=c(t1$x, t2$x), y=c(t1$y, t2$y)))

}


draw.clone <- function(x, y, wid=1, len=1, col='gray',
                       clone.shape='bell',
                       bell.curve.step = 0.25,
                       label=NA, cell.frac=NA,
                       #cell.frac.position='top.out',
                       cell.frac.position='right.mid',
                       cell.frac.top.out.space = 0.75,
                       cell.frac.side.arrow.width=1.5,
                       cell.frac.angle=NULL,
                       cell.frac.side.arrow=TRUE,
                       cell.frac.side.arrow.col='black',
                       variant.names=NULL,
                       variant.color='blue',
                       variant.angle=NULL,
                       text.size=1,
                       border.color='black',
                       border.width=1,
                       wscale=1
                       ){
    beta = min(wid/5, (wid+len)/20)
    gamma = wid/2
    #cat('len=',len, '\n')

    if (clone.shape == 'polygon'){
        xx = c(x, x+beta, x+len, x+len, x+beta)
        yy = c(y, y+gamma, y+gamma, y-gamma, y-gamma)
        polygon(xx, yy, border=border.color, col=col, lwd=border.width)
    }else if(clone.shape == 'bell'){
        beta = min(wid/5, (wid+len)/10, len/3)
        beta = beta*wscale
        xx0= c(x, x+beta, x+len, x+len, x+beta)
        yy0 = c(y, y+gamma, y+gamma, y-gamma, y-gamma)
        #polygon(xx, yy, border='black', col=col, lwd=0.2)

        gamma.shift = min(bell.curve.step, 0.5*gamma)

        # this is to prevent a coeff from being NaN when curve is generated below
        #if (beta <= 0.25){beta = 0.3}
        zeta = min(0.25, max(beta-0.1,0))
        #print(len)

        # shorter time, curve earlier
        zeta = zeta*len/7
        zeta = zeta*wscale

        x0=x+zeta; y0=0; x1=x+beta; y1 = gamma - gamma.shift
        n = 3; n = 1 + len/3
        a = ((y0^n-y1^n)/(x0-x1))^(1/n)
        b = y0^n/a^n - x0
        c = y
        #cat('a=', a, 'b=', b, 'c=', c, 'gamma=', gamma, 'len=', len, 'x0=',
        #    x0, 'x1=', x1, 'y0=', y0, 'y1=', y1,'\n')
        #curve(a*(x+b)^(1/n)+c, n=501, add=TRUE, col=col, xlim=c(x0,x1))
        #curve(-a*(x+b)^(1/n)+c, n=501, add=TRUE, col=col, xlim=c(x0,x1))

        beta0 = beta/5
        if (x0+beta0 > x1){beta0 = (x1-x0)/10}
        gamma0 = gamma/10

        xx = seq(x0+beta0,x1,(x1-x0)/100)
        yy = a*(xx+b)^(1/n)+c
        yy = c(y, yy, y+gamma, y-gamma, -a*(rev(xx)+b)^(1/n)+c)
        xx = c(x, xx, x+len, x+len, rev(xx))
        polygon(xx, yy, border=border.color, col=col, lwd=border.width)

        # generate some points to depict cells
        # buggy, does not work yet, and looks ugly
        if(FALSE){
            cells = generate.fill.points(xx, yy)
            parNew = par('new')
            par(new=TRUE)
            co = par('usr')
            plot(cells$x, cells$y, pch=20, axes=FALSE, col='black', cex=0.1,
                xlim=co[1:2], ylim=co[3:4])
            par(new=parNew)
        }

        #xxx <<- xx; yyy <<- yy; ccc <<- cells
        #pdf('tmp.pdf');plot(xxx,yyy); co =par('usr'); par(new=T); plot(ccc$x, ccc$y, col='red', xlim=co[1:2], ylim=co[3:4], axes=F); dev.off(); dev.off(); dev.off()
        #stop()

    }else if (clone.shape == 'triangle'){
        #TODO: this does not work well yet. Implement!
        xx = c(x, x+len, x+len)
        yy = c(y/10, y+gamma, y-gamma)
        y = y/10
        polygon(xx, yy, border='black', col=col, lwd=0.2)
    }else if (clone.shape == 'parabol'){
        # TODO: Resovle overlapping (ie. subclone parabol expand outside of
        # parent clone parabol)
        x0=x; y0=0; x1=x+len; y1=gamma
        n = 3; n = 1 + len/3
        a = ((y0^n-y1^n)/(x0-x1))^(1/n)
        b = y0^n/a^n - x0
        c = y
        #cat('a=', a, 'b=', b, 'c=', c, 'gamma=', gamma, 'len=', len, 'x0=',
        #    x0, 'x1=', x1, 'y0=', y0, 'y1=', y1,'\n')
        curve(a*(x+b)^(1/n)+c, n=501, add=TRUE, col=col, xlim=c(x0,x1))
        curve(-a*(x+b)^(1/n)+c, n=501, add=TRUE, col=col, xlim=c(x0,x1))
        xx = seq(x0,x1,(x1-x0)/100)
        yy = a*(xx+b)^(1/n)+c
        yy = c(yy, -a*(rev(xx)+b)^(1/n)+c)
        xx = c(xx, rev(xx))
        #print(xx)
        #print(yy)
        polygon(xx, yy, col=col)
    }


    if (!is.na(label)){
        text(x+0.2*text.size*len/3.5, y, label, cex=text.size, adj=c(0,0.5))
    }
    if (!is.na(cell.frac)){
        cell.frac.x = 0
        cell.frac.y = 0
        angle = 0
        adj = c(0,0)
        if (cell.frac.position == 'top.left'){
            cell.frac.x = max(x+beta, x + 0.4)
            cell.frac.y = y+gamma#-0.3*text.size
            adj = c(0, 1)
        }else if (cell.frac.position == 'top.right'){
            cell.frac.x = x+len
            cell.frac.y = y+gamma
            adj = c(1, 1)
        }else if (cell.frac.position == 'top.mid'){
            cell.frac.x = x+beta+(len-beta)/2
            cell.frac.y = y+gamma
            adj = c(0.5, 1)
        }else if (cell.frac.position == 'right.mid'){
            cell.frac.x = x+len
            cell.frac.y = y
            adj = c(0, 0.5)
            angle = 45
        }else if (cell.frac.position == 'right.top'){
            cell.frac.x = x+len
            cell.frac.y = y+gamma
            adj = c(0, 1)
            angle = 45
        }else if (cell.frac.position == 'side'){
            angle = atan(gamma/beta)*(180/pi)# - 5
            cell.frac.x = x+beta/3+0.3*text.size
            cell.frac.y = y+gamma/2-0.3*text.size
            adj = c(0.5, 0.5)
        }else if (cell.frac.position == 'top.out'){
            cell.frac.x = x+len
            cell.frac.y = y.out
            adj = c(1, 0.5)
            # increase y.out so next time, text will be plotted a little higher
            # to prevent overwritten! Also, x.out.shift is distance from arrow
            # to polygon
            y.out <<- y.out + cell.frac.top.out.space
            x.out.shift <<- x.out.shift + 0.1
        }

        if (cell.frac.position == 'top.right' && clone.shape == 'bell'){
            angle = atan(gamma.shift/len*w2h.scale)*(180/pi)
        }
        #debug
        #cat('x=', cell.frac.x, 'y=', cell.frac.y, '\n')
        if (is.null(cell.frac.angle)){
            cell.frac.angle = angle
        }
        text(cell.frac.x, cell.frac.y, cell.frac,
             cex=text.size*0.7, srt=cell.frac.angle, adj=adj)
        if(cell.frac.side.arrow && cell.frac.position=='top.out'){
            # draw arrow
            x0 = cell.frac.x
            y0 = cell.frac.y
            x1 = cell.frac.x+x.out.shift
            y1 = y0
            x2 = x2 = x1
            y2 = y+gamma
            x3 = x0
            y3 = y2
            segments(x0, y0, x1, y1, col=cell.frac.side.arrow.col,
                     lwd=cell.frac.side.arrow.width)
            segments(x1, y1, x2, y2, col=cell.frac.side.arrow.col,
                     lwd=cell.frac.side.arrow.width)
            arrows(x0=x2, y0=y2, x1=x3, y1=y3,col=cell.frac.side.arrow.col,
                   length=0.025, lwd=cell.frac.side.arrow.width)
        }

        if (!is.null(variant.names)){
            if (is.null(variant.angle)){
                variant.angle = atan(gamma/beta)*(180/pi)# - 5
            }
            variant.x = x+beta/3+0.3*text.size
            variant.y = y+gamma/2-1.5*text.size
            variant.adj = c(0.5, 1)
            text(variant.x, variant.y, paste(variant.names, collapse='\n'),
                 cex=text.size*0.54, srt=variant.angle, adj=variant.adj,
                 col=variant.color)
        }
    }
}


 rescale.vaf <- function(v, down.scale=0.99){
    #v = vx
    #print(v)
    #cat('Scaling called.\n')
    rescale <- function(i){
        #print(i)
        parent = v[i,]$lab
        parent.vaf = v[i, ]$vaf
        subclones.idx = which(v$parent == parent)
        sum.sub.vaf = sum(v[subclones.idx,]$vaf)
        scale = ifelse(sum.sub.vaf > 0, parent.vaf/sum.sub.vaf, 1)
        #debug
        #cat('parent.vaf=', parent.vaf, ';sum.sub.vaf=', sum.sub.vaf,
        #    ';scale=', scale, '\n')
        for (idx in subclones.idx){
            if(scale < 1){
                #cat('Scaling...\n')
                #print(v)
                v[idx,]$vaf <<- down.scale*scale*v[idx,]$vaf
                #print(v)
            }
        }
        for (idx in subclones.idx){
            rescale(idx)
        }

    }
    root.idx = which(v$parent == -1)
    rescale(root.idx)
    return(v)
}


set.position <- function(v){
    v$y.shift = 0
    max.vaf = max(v$vaf)
    scale = 0.5/max.vaf
    #debug
    for (i in 1:nrow(v)){
        vi = v[i,]
        subs = v[!is.na(v$parent) & v$parent == vi$lab,]
        if (nrow(subs) == 0){next}
        vafs = subs$vaf
        margin = (vi$vaf - sum(vafs))/length(vafs)*scale
        sp = 0
        if (margin > 0){
            margin = margin*0.75
            sp = margin*0.25
        }
        spaces = rep(sp, length(vafs))
        if (length(spaces) >= 2){
            for (j in 2:length(spaces)){
                spaces[j] = sum(vafs[1:j-1]+margin)
            }
        }else{
            # re-centering if only 1 subclone inside another
            spaces = (vi$vaf-vafs)/2
        }
        #debug
        #print(subs)
        v[!is.na(v$parent) & v$parent == vi$lab,]$y.shift = spaces
    }
    #print(v)
    return(v)
}


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


get.cell.frac.ci <- function(vi, include.p.value=TRUE, sep=' - '){
    cell.frac = NULL
    is.zero = NULL
    is.subclone = NULL
    if('free.lower' %in% colnames(vi)){
        cell.frac.lower = ifelse(vi$free.lower == 0, '0',
                             gsub('\\.[0]+$|0+$', '',
                                  sprintf('%0.1f', 200*vi$free.lower)))
        cell.frac.upper = ifelse(vi$free.upper >= 0.5, '100%',
                             gsub('\\.[0]+$|0+$', '',
                                  sprintf('%0.1f%%', 200*vi$free.upper)))
        cell.frac = paste0(cell.frac.lower, sep , cell.frac.upper)
        if(include.p.value){
            cell.frac = paste0(cell.frac, '(',sprintf('%0.2f',
                            vi$free.confident.level),
                            #',p=', 1-vi$p.value,
                           ')')
            cell.frac = paste0(cell.frac, '/p=', sprintf('%0.3f', vi$p.value))
        }
        rownames(vi) = vi$lab

        # if only one clone as root, all cell.frac is NA
        # this is a dirty fix for output display
        # TODO: assign cell.frac for clone with zero subclone when
        # enumarating the models in enumerate.clones function.
        if(all(is.na(cell.frac.lower))){
            cell.frac[1] = '100-100%'
        }

    #if ('free.lower' %in% colnames(vi)){
        is.zero = ifelse(vi$free.lower >= 0, FALSE, TRUE)
        rownames(vi) = vi$lab
        names(is.zero) = vi$lab
        is.subclone = determine.subclone(vi,
            vi$lab[!is.na(vi$parent) & vi$parent == '-1'])$is.sub
    #}
    }

    #debug
    #print(vi$free.confident.level.non.negative)
    #if (vi$free.confident.level.non.negative == 0.687){
    #    print(vi$free.confident.level.non.negative)
    #    print(cell.frac)
    #    print(vi)
    #    vii <<- vi
    #}

    return(list(cell.frac.ci=cell.frac, is.zero.cell.frac=is.zero, is.subclone=is.subclone))
}


draw.sample.clones <- function(v, x=1, y=0, wid=30, len=9,
                               clone.shape='bell',
                               bell.curve.step=0.25,
                               bell.border.width=1,
                               clone.time.step.scale=1,
                               label=NULL, text.size=1,
                               cell.frac.ci=FALSE,
                               disable.cell.frac=FALSE,
                               zero.cell.frac.clone.color=NULL,
                               zero.cell.frac.clone.border.color=NULL,
                               nonzero.cell.frac.clone.border.color=NULL,
                               nonzero.cell.frac.clone.border.width=NULL,
                               zero.cell.frac.clone.border.width=NULL,
                               top.title=NULL,
                               clone.height=TRUE,
                               cell.frac.top.out.space=0.75,
                               cell.frac.side.arrow.width=1.5,
                               variants.to.highlight=NULL,
                               variant.color='blue',
                               variant.angle=NULL,
                               show.time.axis=TRUE,
                               color.node.by.sample.group=FALSE,
                               color.border.by.sample.group=TRUE,
                               show.clone.label=TRUE,
                               wscale=1){
    v = v[!v$excluded,]
    if (adjust.clone.height){
        #cat('Will call rescale.vaf on', label, '\n')
        #print(v)
        v = rescale.vaf(v)

    }
    # scale VAF so that set.position works properly, and drawing works properly
    max.vaf = max(v$vaf)
    scale = 0.5/max.vaf
    v$vaf = v$vaf*scale
    max.vaf = max(v$vaf)
    high.vaf = max.vaf - 0.02
    low.vaf = 0.2
    y.out <<- wid*max.vaf/2+0.5
    x.out.shift <<- 0.1

    wscale = wscale*clone.time.step.scale

    #print(v)


    draw.sample.clone <- function(i){
        vi = v[i,]
        #debug
        #cat('drawing', vi$lab, '\n')
        if (vi$vaf > 0){
            #if (vi$parent == 0){# root
            #if (is.na(vi$parent)){
            if (!is.na(vi$parent) && vi$parent == -1){
                xi = x
                yi = y
                leni = len
            }else{
                # for bell curve, needs to shift x further to make sure
                # bell of subclone falls completely in its parent bell
                x.shift = 1 * ifelse(clone.shape=='bell', 1.2, 1)

                if (vi$y.shift + vi$vaf >= high.vaf && vi$vaf < low.vaf){
                    x.shift = 2*x.shift
                }
                if (clone.shape=='triangle'){
                    x.shift = x.shift + 1
                }
                par = v[v$lab == vi$parent,]

                if (vi$vaf < 0.05 && par$num.subclones > 1){x.shift = x.shift*2}
                x.shift = x.shift*clone.time.step.scale*len/7*wscale
                xi = par$x + x.shift

                yi = par$y - wid*par$vaf/2 + wid*vi$vaf/2 + vi$y.shift*wid
                leni = par$len - x.shift
            }
            #cell.frac.position = ifelse(vi$free.lower < 0.05 & vi$vaf > 0.25, 'side', 'top.right')
            #cell.frac.position = ifelse(vi$free.lower < 0.05, 'top.out', 'top.right')
            cell.frac.position = ifelse(vi$free < 0.05, 'top.out', 'top.right')
            #cell.frac.position = ifelse(vi$free < 0.05, 'top.out', 'right.mid')
            #cell.frac.position = ifelse(vi$free < 0.05, 'top.out', 'top.out')
            #cell.frac.position = ifelse(vi$num.subclones > 0 , 'right.top', 'right.mid')
            #cell.frac.position = 'top.mid'
            cell.frac = paste0(gsub('\\.[0]+$|0+$', '',
                                    sprintf('%0.2f', vi$free.mean*2*100)), '%')
            if(cell.frac.ci && !disable.cell.frac){
                cell.frac = get.cell.frac.ci(vi, include.p.value=TRUE)$cell.frac.ci
            }else if (disable.cell.frac){
                cell.frac = NA
            }
            variant.names = variants.to.highlight$variant.name[
                variants.to.highlight$cluster == vi$lab]
            if (length(variant.names) == 0) {
                variant.names = NULL
            }
            clone.color = vi$color
            border.color='black'
            if (color.border.by.sample.group){
                border.color = vi$sample.group.color
            }else if (color.node.by.sample.group){
                clone.color = vi$sample.group.color
            }
            if (!is.null(zero.cell.frac.clone.border.color) & vi$is.zero){
                border.color = zero.cell.frac.clone.border.color
                if (border.color == 'fill'){border.color = clone.color}
            }
            if (!is.null(zero.cell.frac.clone.color) & vi$is.zero){
                clone.color = zero.cell.frac.clone.color
            }

            if (!is.null(nonzero.cell.frac.clone.border.color) & !vi$is.zero){
                border.color = nonzero.cell.frac.clone.border.color
                if (border.color == 'fill'){border.color = clone.color}
            }
            if (!is.null(nonzero.cell.frac.clone.border.width) & !vi$is.zero){
                bell.border.width = nonzero.cell.frac.clone.border.width
            }
            if (!is.null(zero.cell.frac.clone.border.width) & vi$is.zero){
                bell.border.width = zero.cell.frac.clone.border.width
            }

            clone.label = ""
            if (show.clone.label){clone.label = vi$lab}
            draw.clone(xi, yi, wid=wid*vi$vaf, len=leni, col=clone.color,
                       clone.shape=clone.shape,
                       bell.curve.step=bell.curve.step,
                       border.width=bell.border.width,
                       #label=vi$lab,
                       label=clone.label,
                       cell.frac=cell.frac,
                       cell.frac.position=cell.frac.position,
                       cell.frac.side.arrow.col=clone.color,
                       text.size=text.size,
                       cell.frac.top.out.space=cell.frac.top.out.space,
                       cell.frac.side.arrow.width=cell.frac.side.arrow.width,
                       variant.names=variant.names,
                       variant.color=variant.color,
                       variant.angle=variant.angle,
                       border.color=border.color,
                       wscale=wscale)
            v[i,]$x <<- xi
            v[i,]$y <<- yi
            v[i,]$len <<- leni
            for (j in 1:nrow(v)){
                #cat('---', v[j,]$parent,'\n')
                if (!is.na(v[j,]$parent) && v[j,]$parent != -1 &&
                        v[j,]$parent == vi$lab){
                    draw.sample.clone(j)
                }
            }
        }
        # draw time axis
        if (show.time.axis && i==1){
            axis.y = -9
            arrows(x0=x,y0=axis.y,x1=10,y1=axis.y, length=0.05, lwd=0.5)
            text(x=10, y=axis.y-0.75, label='time', cex=1, adj=1)
            segments(x0=x,y0=axis.y-0.2,x1=x, y1=axis.y+0.2)
            text(x=x,y=axis.y-0.75,label='Cancer initiated', cex=1, adj=0)
            segments(x0=x+len,y0=axis.y-0.2,x1=x+len, y1=axis.y+0.2)
            text(x=x+len, y=axis.y-0.75, label='Sample taken', cex=1, adj=1)
        }
    }
    plot(c(0, 10),c(-10,10), type = "n", xlab='', ylab='', xaxt='n',
         yaxt='n', axes=FALSE)
    if (!is.null(label)){
        text(x-1*wscale, y, label=label, srt=90, cex=text.size, adj=c(0.5,1))
        #text(x, y, label=label, srt=90, cex=text.size, adj=c(0.5,1))
    }
    if (!is.null(top.title)){
        text(x, y+10, label=top.title, cex=(text.size), adj=c(0,0.5))
    }

    # move root to the first row and plot
    root = v[!is.na(v$parent) & v$parent == -1,]
    v = v[is.na(v$parent) | v$parent != -1,]
    v = rbind(root, v)
    v = set.position(v)
    v$x = 0
    v$y = 0
    v$len = 0

    #debug
    #print(v)

    draw.sample.clone(1)
}


make.graph <- function(v, cell.frac.ci=TRUE, node.annotation='clone', node.colors=NULL){
    library(igraph)
    #v = v[!is.na(v$parent),]
    #v = v[!is.na(v$parent) | v$vaf != 0,]
    v = v[!is.na(v$parent),]
    #rownames(v) = seq(1,nrow(v))
    rownames(v) = v$lab
    g = matrix(0, nrow=nrow(v), ncol=nrow(v))
    rownames(g) = rownames(v)
    colnames(g) = rownames(v)
    if (nrow(v) == 0){# single sample, no pruned tree
        return(NULL)
    }
    for (i in 1:nrow(v)){
        par.lab = v[i,]$parent
        if (!is.na(par.lab) && par.lab != -1){
            #if(par.lab != 0){
            par.idx = rownames(v[v$lab == par.lab,])
            #debug
            #cat(par.lab, '--', par.idx, '\n')
            g[par.idx, i] = 1
            #}
        }
    }
    #print(g)
    g <- graph.adjacency(g)
    cell.frac = gsub('\\.[0]+$|0+$', '', sprintf('%0.2f%%', v$free.mean*2*100))
    if(cell.frac.ci){
        cell.frac = get.cell.frac.ci(v, include.p.value=FALSE,
                                     sep=' -\n')$cell.frac.ci
    }
    labels = v$lab
    colors = v$color
    if (!is.null(node.colors)){
        colors = node.colors[labels]
    }
    if (!is.null(cell.frac) && !all(is.na(cell.frac))){
        labels = paste0(labels,'\n', cell.frac)
    }

    # add sample name
    # trick to strip off clone having zero cell.frac and not a founding clone of a sample
    # those samples prefixed by 'o*'
    remove.founding.zero.cell.frac = FALSE
    if (node.annotation == 'sample.with.cell.frac.ci.founding.and.subclone'){
        node.annotation = 'sample.with.cell.frac.ci'
        remove.founding.zero.cell.frac = TRUE
    }
    if (node.annotation != 'clone' && node.annotation %in% colnames(v)){
        # this code is to add sample to its terminal clones only, obsolete
        #leaves = !grepl(',', v$sample)
        #leaves = v$is.term
        #if (any(leaves)){
        #    #labels[leaves] = paste0('\n', labels[leaves], '\n', v$leaf.of.sample[leaves])
        #}
        has.sample = !is.na(v[[node.annotation]])
        samples.annot = v[[node.annotation]][has.sample]
        if (!cell.frac.ci){
            # strip off cell.frac.ci
            #tmp = unlist(strsplit(',', samples.annot))
            #tmp = gsub('\\s*:\\s*[^:].+', ',', tmp)
            samples.annot = gsub('\\s*:\\s*[^:]+(,|$)', ',', v[[node.annotation]][has.sample])
            samples.annot = gsub(',$', '', samples.annot)
        }
        if (remove.founding.zero.cell.frac){
            samples.annot = gsub('o\\*[^,]+(,|$)', '', samples.annot)
        }
        labels[has.sample] = paste0(labels[has.sample], '\n', samples.annot)
    }
    V(g)$name = labels
    V(g)$color = colors
    return(list(graph=g, v=v))
}




draw.sample.clones.all <- function(x, outPrefix, object.to.plot='polygon',
                                   ignore.clusters=NULL){
    pdf(paste0(outPrefix, '.pdf'), width=6, height=6)
    for(i in 1:length(x)){
        xi = x[[i]]
        #xi = scale.cell.frac(xi, ignore.clusters=ignore.clusters)
        if (object.to.plot == 'polygon'){
            draw.sample.clones(xi, cell.frac.ci=TRUE)
        }else{
            plot.tree(xi, node.shape='circle', node.size=35, cell.frac.ci=TRUE)
        }
    }
    dev.off()
    cat(outPrefix, '\n')
}


plot.tree <- function(v, node.shape='circle', display='tree',
                      node.size=50,
                      node.colors=NULL,
                      color.node.by.sample.group=FALSE,
                      color.border.by.sample.group=TRUE,
                      show.legend=TRUE,
                      tree.node.text.size=1,
                      cell.frac.ci=TRUE,
                      node.prefix.to.add=NULL,
                      title='',
                      #show.sample=FALSE,
                      node.annotation='clone',
                      node.label.split.character=NULL,
                      node.num.samples.per.line=NULL,
                      out.prefix=NULL,
                      graphml.out=FALSE,
                      out.format='graphml'){
    library(igraph)
    grps = NULL
    grp.colors = 'black'
    if (color.border.by.sample.group){
        color.node.by.sample.group = FALSE #disable coloring node by group if blanket is used
        #grps = list()
        #for (i in 1:nrow(v)){
        #    grps = c(grps, list(i))
        #}
        grp.colors = v$sample.group.color
        # get stronger color for borders
        #uniq.colors = unique(grp.colors)
        #border.colors = get.clonevol.colors(length(uniq.colors), T)
        #names(border.colors) = uniq.colors
        #grp.colors = border.colors[grp.colors]
        #v$sample.group.border.color = grp.colors
    }else if (color.node.by.sample.group){
        node.colors = v$sample.group.color
        names(node.colors) = v$lab
    }
    #x = make.graph(v, cell.frac.ci=cell.frac.ci, include.sample.in.label=show.sample, node.colors)
    x = make.graph(v, cell.frac.ci=cell.frac.ci, node.annotation=node.annotation, node.colors)
    #print(v)
    g = x$graph
    v = x$v
    root.idx = which(!is.na(v$parent) & v$parent == '-1')
    #cell.frac = gsub('\\.[0]+$|0+$', '', sprintf('%0.2f%%', v$free*2*100))


    #V(g)$color = v$color
    #display = 'graph'
    if(display == 'tree'){
        layout = layout.reingold.tilford(g, root=root.idx)
    }else{
        layout = NULL
    }

    vertex.labels = V(g)$name
    #vlabs <<- vertex.labels
    if (!is.null(node.label.split.character)){
        num.splits = sapply(vertex.labels, function(l)
            nchar(gsub(paste0('[^', node.label.split.character, ']'), '', l)))

        # only keep the node.label.split.char in interval of node.num.samples.per.line
        # such that a block of node.num.samples.per.line samples will be grouped and
        # kept in one line
        if (!is.null(node.num.samples.per.line)){
            for (i in 1:length(vertex.labels)){
                vl = unlist(strsplit(vertex.labels[i], node.label.split.character))
                sel = seq(min(node.num.samples.per.line, length(vl)),
                    length(vl),node.num.samples.per.line)
                vl[sel] = paste0(vl[sel], node.label.split.character)
                vl[-sel] = paste0(vl[-sel], ';')
                vertex.labels[i] = paste(vl, collapse='')
            }
            num.splits = length(sel) + 1
        }
        extra.lf = sapply(num.splits, function(n) paste(rep('\n', n), collapse=''))
        vertex.labels = paste0(extra.lf, gsub(node.label.split.character, '\n',
            vertex.labels))
    }

    plot(g, edge.color='black', layout=layout, main=title,
         edge.arrow.size=0.75, edge.arrow.width=0.75,
         vertex.shape=node.shape, vertex.size=node.size,
         vertex.label.cex=tree.node.text.size,
         #vertex.label.color=sample(c('black', 'blue', 'darkred'), length(vertex.labels), replace=T),
         vertex.label=vertex.labels,
         #mark.groups = grps,
         #mark.col = 'white',
         #mark.border = grp.colors,
         vertex.frame.color=grp.colors)
         #, vertex.color=v$color, #vertex.label=labels)
    if ((color.node.by.sample.group || color.border.by.sample.group) & show.legend &
            'sample.group' %in% colnames(v)){
        vi = unique(v[!v$excluded & !is.na(v$parent),
            c('sample.group', 'sample.group.color')])
        vi = vi[order(vi$sample.group),]
        if (color.border.by.sample.group){
            legend('topright', legend=vi$sample.group, pt.cex=3, cex=1.5,
                 pch=1, col=vi$sample.group.color)
        }else{
            legend('topright', legend=vi$sample.group, pt.cex=3, cex=1.5,
                 pch=16, col=vi$sample.group.color)
        }
        legend('bottomleft', legend=c('*  sample founding clone',
                                    '\u00B0  zero cellular fraction',
                                   '\u00B0* ancestor of sample founding clone'
                                  ),
                                  pch=c('', '', ''))
        # events on each clone legend
        if ('events' %in% colnames(v)){
            ve = v[v$events != '',]
            ve = ve[order(as.integer(ve$lab)),]
            # only print 5 events per-line
            ve$events = insert.lf(ve$events, 5, ',')
            legend('topleft', legend=paste0(sprintf('%2s', ve$lab), ': ', ve$events),
                    pt.cex=2, cex=1, pch=19, col=ve$color)
        }
    }

    # remove newline char because Cytoscape does not support multi-line label
    V(g)$name = gsub('\n', ' ', V(g)$name, fixed=TRUE)
    if (!is.null(node.prefix.to.add)){
        V(g)$name = paste0(node.prefix.to.add, V(g)$name)
    }

    if (!is.null(out.prefix)){
        out.file = paste0(out.prefix, '.', out.format)
        #cat('Writing tree to ', out.file, '\n')
        if (graphml.out){
            write.graph(g, file=out.file, format=out.format)
        }
    }

    return(g)

}


write.tree <- function(v, out.file, out.format='tabular'){
    v = v[, c('parent', 'lab', '')]
}

get.model.score <- function(v){
    #return(prod(v$p.value[!is.na(v$p.value)]))
    return(max(v$p.value[!is.na(v$p.value)]))
}

get.model.max.p.value <- function(v){
    #return(prod(v$p.value[!is.na(v$p.value)]))
    return(max(v$p.value[!is.na(v$p.value)]))
}


get.model.non.negative.ccf.prob <- function(v){
    return(prod(1 - v$p.value[!is.na(v$p.value)]))
}



get.subclones.across.samples <- function(x, matched.model.index){
    samples = names(x$models)
    tree = x$matched$merged.trees[[matched.model.index]]
    labs = tree$lab[!tree$excluded]
    # look for subclones in each sample for each clone
    subs = NULL
    for (s in samples){
        m = x$models[[s]][[x$matched$index[matched.model.index, s]]]
        for (cl in labs){
            sc = m$lab[!is.na(m$parent) & m$parent == cl & !m$excluded]
            p = m$p.value[m$lab == cl]
            if (length(sc) > 0){
                sc = paste(sort(sc), collapse=',')
                r = data.frame(lab=cl, sample=s, subclones=sc, p=p,
                    stringsAsFactors=FALSE)
                if(is.null(subs)){subs = r}else{subs = rbind(subs,r)}
            }
        }
    }

    subs = subs[order(subs$lab),]
    return(subs)
}


cross.rule.score <- function(x, meta.p.method='fisher', exhaustive.mode=FALSE,
                             rank=TRUE, boot=NULL){
    if (!is.null(x$matched) && x$num.matched.models > 0 && ncol(x$matched$index) > 1){
        samples = names(x$models)
        num.models = nrow(x$matched$index)
        x$matched$scores$max.clone.ccf.combined.p = NA
        x$matched$clone.ccf.pvalues = list()
        # foreach matched model, recalc score by combining p values across
        # samples for each clone
        for (i in 1:num.models){
            trees = NULL
            t = x$match$merged.trees[[i]]
            p = NULL
            for (s in samples){
                mi = x$models[[s]][[x$match$index[i,s]]]
                mi = mi[!mi$excluded & !is.na(mi$parent), c('lab', 'p.value')]
                colnames(mi) = c('lab', s)
                if (is.null(p)){
                    p = mi
                }else{
                    p = merge(p, mi, all=TRUE)
                }
            }
            #ppp <<- p
            if (ncol(p) == 2){#single sample
                p$cmb.p = apply(p[,c(2,2)], 1, combine.p, method=meta.p.method)
            }else{
                p$cmb.p = apply(p[,-1], 1, combine.p, method=meta.p.method)
            }
            # model score = max (combined p of each clone)
            x$matched$scores$max.clone.ccf.combined.p[i] = max(p$cmb.p)
            x$matched$merged.trees[[i]]$clone.ccf.combined.p = p$cmb.p
            # save the whole pvalue matrix
            x$matched$clone.ccf.pvalues[[i]] = p
        }
        # order matched models by new score
        idx = seq(1,nrow(x$matched$scores))
        if (rank){
            idx = order(x$matched$scores$max.clone.ccf.combined.p)
        }
        x$matched$index = x$matched$index[idx,, drop=FALSE]
        x$matched$scores = x$matched$scores[idx,, drop=FALSE]
        x$matched$probs = x$matched$probs[idx,, drop=FALSE]
        # order merged trees
        tmp = list()
        for (i in idx){
            tmp = c(tmp, list(x$matched$merged.trees[[i]]))
        }
        x$matched$merged.trees = tmp
        # order merged traces
        tmp = list()
        for (i in idx){
            tmp = c(tmp, list(x$matched$merged.traces[[i]]))
        }
        x$matched$merged.traces = tmp
        # order pvalues
        tmp = list()
        for (i in idx){
            tmp = c(tmp, list(x$matched$clone.ccf.pvalues[[i]]))
        }
        x$matched$clone.ccf.pvalues = tmp

        # remove previous obsolete model scores (which was very small probability)
        #x$matched$scores$model.prob = x$matched$scores$model.score
        #x$matched$scores$model.score = NULL

    }
    return(x)
}


merge.clone.trees <- function(trees, samples=NULL, sample.groups=NULL,
                              merge.similar.samples=FALSE){
    # to keep track of what samples merged with what samples
    merged.trace = NULL

    if (merge.similar.samples){
        # remove sample having tree similar to tree of another sample
        z = trim.clone.trees(trees, samples,
                             remove.sample.specific.clones=FALSE)
        merged.trace = z$merged.trace
        trees = z$unique.trees
        samples = names(trees)
        sample.groups = sample.groups[samples]
    }
    n = length(trees)
    merged = NULL
    if (is.null(samples)){samples = seq(1,n)}
    #leaves = c()
    lf = NULL
    ccf.ci = NULL
    ccf.ci.nonzero = NULL
    subclones = NULL #subclonal status
    cgrp = NULL #grouping clones based on sample groups
    #let's group all samples in one group if sample groups not provided
    if (is.null(sample.groups)){
        sample.groups = rep('group1', length(samples))
        names(sample.groups) = samples
    }

    #TODO: there is a bug in infer.clonal.models that did not give
    # consistent ancestors value across samples, let's discard this
    # column now for merging, but later need to fix this.
    #v = trees[[i]][, c('lab', 'color', 'parent', 'ancestors', 'excluded')]
    key.cols = c('lab', 'color', 'parent', 'excluded')

    for (i in 1:n){
        v = trees[[i]]
        s = samples[i]

        # get cell.frac
        v = v[!v$excluded & !is.na(v$parent),]
        if (nrow(v) == 0){stop('ERROR: Something wrong. No clones left after filter. They might have been excluded.\n')}
        # TODO: scale.cell.frac here works independent of plot.clonal.models
        # which has a param to ask for scaling too. Make them work together nicely.
        # also, clonal tree data frame v is now have two more column indicating
        # if a clone is.subclone or is.zero cell frac, so getting these info via
        # get.cell.frac.ci is redundant and potentially create inconsistency if code changes
        # TODO: utilize is.subclone and is.zero columns in v
        cia = get.cell.frac.ci(scale.cell.frac(v), sep='-')
        #ci = data.frame(lab=v$lab, sample.with.cell.frac.ci=paste0(ifelse(cia$is.subclone,
        #    '', '*'), s, ' : ', cia$cell.frac.ci), stringsAsFactors=F)
        ci = data.frame(lab=v$lab, sample.with.cell.frac.ci=paste0(ifelse(v$is.founder,
            '*', ''), s, ' : ', cia$cell.frac.ci), stringsAsFactors=FALSE)
        #ci.nonzero = ci[!is.na(cia$is.zero.cell.frac) & !cia$is.zero.cell.frac,]
        ci.nonzero = ci[!cia$is.zero.cell.frac,]
        ci$sample.with.cell.frac.ci[cia$is.zero.cell.frac] = paste0('\u00B0',
            ci$sample.with.cell.frac.ci[cia$is.zero.cell.frac])
        if (is.null(ccf.ci)){ccf.ci = ci}else{ccf.ci = rbind(ccf.ci, ci)}
        if (is.null(ccf.ci.nonzero)){ccf.ci.nonzero = ci.nonzero}else{
            ccf.ci.nonzero = rbind(ccf.ci.nonzero, ci.nonzero)}
        v$sample = s
        #v$sample[!cia$is.subclone] = paste0('*', v$sample[!cia$is.subclone])
        v$sample[v$is.founder] = paste0('*', v$sample[v$is.founder])
        v$sample[cia$is.zero.cell.frac] = paste0('\u00B0', v$sample[cia$is.zero.cell.frac])
        this.leaves = v$lab[!is.na(v$parent) & !(v$lab %in% v$parent)]
        this.lf = data.frame(lab=this.leaves, leaf.of.sample=s,
                             stringsAsFactors=FALSE)
        if (is.null(lf)){lf = this.lf}else{lf = rbind(lf, this.lf)}
        #leaves = c(leaves, this.leaves)

        #clone grouping only non.zero cell.frac clones
        vz = v[!cia$is.zero.cell.frac,]
        if (!is.null(sample.groups)){#this is uneccesary if given default grouping above
            cg = data.frame(lab=vz$lab, sample.group=sample.groups[s],
                stringsAsFactors=FALSE, row.names=NULL)
            if (is.null(cgrp)){cgrp = cg}else{cgrp = rbind(cgrp, cg)}
        }

        # keep only key.cols and sample cols for merging
        v = v[, c(key.cols, 'sample')]
        if (is.null(merged)){merged = v}else{merged = rbind(merged, v)}
    }
    merged = merged[!is.na(merged$parent),]

    #merged = unique(merged)
    merged = aggregate(sample ~ ., merged, paste, collapse=',')
    #leaves = unique(leaves)
    lf = aggregate(leaf.of.sample ~ ., lf, paste, collapse=',')
    lf$is.term = TRUE
    merged = merge(merged, lf, all.x=TRUE)

    ccf.ci = aggregate(sample.with.cell.frac.ci ~ ., ccf.ci, paste, collapse=',')
    merged = merge(merged, ccf.ci, all.x=TRUE)

    ccf.ci.nonzero = aggregate(sample.with.cell.frac.ci ~ ., ccf.ci.nonzero,
                                paste, collapse=',')
    colnames(ccf.ci.nonzero) = c('lab', 'sample.with.nonzero.cell.frac.ci')
    merged = merge(merged, ccf.ci.nonzero, all.x=TRUE)


    if (!is.null(cgrp)){
        #print(cgrp)
        cgrp = unique(cgrp)
        cgrp = cgrp[order(cgrp$sample.group),]
        cgrp = aggregate(sample.group ~ ., cgrp, paste, collapse=',')
        #print(cgrp)
        sample.grps = unique(cgrp$sample.group)
        sample.group.colors = get.clonevol.colors(length(sample.grps),
                                                  strong.color=TRUE)
        names(sample.group.colors) = sample.grps
        cgrp$sample.group.color = sample.group.colors[cgrp$sample.group]
        merged = merge(merged, cgrp, all.x=TRUE)
    }


    #leaves = unique(lf$lab)
    #merged$is.term = FALSE
    #merged$is.term[merged$lab %in% leaves] = TRUE
    merged$is.term[is.na(merged$is.term)] = FALSE
    merged$num.samples = sapply(merged$sample, function (l)
        length(unlist(strsplit(l, ','))))
    merged$leaf.of.sample.count = sapply(merged$leaf.of.sample, function (l)
        length(unlist(strsplit(l, ','))))
    merged$num.samples[is.na(merged$num.samples)] = 0
    merged$leaf.of.sample.count[is.na(merged$leaf.of.sample)] = 0
    rownames(merged) = merged$lab
    return (list(merged.tree=merged, merged.trace=merged.trace))
}



compare.clone.trees <- function(v1, v2, compare.seeding.models=TRUE){
    res = FALSE
    if (nrow(v1) == nrow(v2) &&
        all(v1$parent == v2$parent) &&
        all(v1$lab == v2$lab)){
        if (compare.seeding.models){
            #cat('\n**** Compare seeding models.\n')
            res = all(v1$sample == v2$sample)
        }else{
             res = TRUE
        }

    }
    return(res)
}


trim.clone.trees <- function(merged.trees, remove.sample.specific.clones=TRUE, samples=NULL,
                                seeding.aware.tree.pruning=TRUE){
    cat('Seeding aware pruning is: ', ifelse(seeding.aware.tree.pruning, 'on\n', 'off\n'))
    n = length(merged.trees)
    if (n == 0){
        return(list(unique.trees=NULL, merged.trace=NULL))
    }
    # trim off sample specific clones (if requested) and excluded nodes, sort by label
    for (i in 1:n){
        v = merged.trees[[i]]
        v = v[!v$excluded,]
        if (remove.sample.specific.clones){
            v = v[v$num.samples > 1,]
        }
        v = v[order(v$lab),]
        merged.trees[[i]] = v
    }

    #ttt2 <<- merged.trees

    # compare all pair of merged.trees and eliminate trees
    # that are already presented
    idx = seq(1,n)
    i = 1;
    merged.trace = c(1,1)
    while (i < n){
        j = i + 1
        if (i > 1){merged.trace = rbind(merged.trace, c(idx[i],idx[i]))}
        while (j <= n){
            #cat(samples[idx[i]], samples[idx[j]], '\t')
            if(compare.clone.trees(merged.trees[[i]], merged.trees[[j]],
                                    compare.seeding.models=seeding.aware.tree.pruning)){
                #cat('Drop tree', j, '\n')
                merged.trees[[j]] = NULL
                n = n - 1
                merged.trace = rbind(merged.trace, c(idx[i], idx[j]))
                idx = idx[-j]
                #cat('equal\n')
            }else{
                j = j + 1
                #cat('diff\n')
            }
        }
        i = i + 1
    }
    #print(idx)
    #print(samples)

    # today debug
    # mtr <<- merged.trace
    if (is.null(dim(merged.trace))){# one row, 1 tree to merge,
        #need to make matrix for followed code to work
        merged.trace = matrix(merged.trace, nrow=1)
    }

    colnames(merged.trace) = c('sample', 'similar.sample')
    # last sample that were not merged with any other
    merged.trace = as.data.frame.matrix(merged.trace)
    rownames(merged.trace) = NULL
    if (!(idx[n] %in% c(merged.trace$sample, merged.trace$similar.sample))){
        merged.trace = rbind(merged.trace, c(idx[n], idx[n]))
        #cat(idx[n], '\n')
        #stop('HMMM')
    }

    if (!is.null(samples)){
        names(merged.trees) = samples[idx]
        merged.trace$sample = samples[merged.trace$sample]
        merged.trace$similar.sample = samples[merged.trace$similar.sample]
    }


    return(list(unique.trees=merged.trees, merged.trace=merged.trace))
}


compare.clone.trees.removing.leaves <- function(v1, v2, ignore.seeding=FALSE){
    res = 2
    v1 = v1[!v1$excluded,]
    v2 = v2[!v1$excluded,]
    v1 = v1[order(v1$lab),]
    v2 = v2[order(v2$lab),]

    if (nrow(v1) == nrow(v2)){
        if (all(v1$parent == v2$parent)){
            res = 0
        }
    }
    if (res !=0){
        # remove leaves
        v1 = v1[!is.na(v1$parent) & (v1$lab %in% v1$parent),]
        v2 = v2[!is.na(v2$parent) & (v2$lab %in% v2$parent),]
        if (all(v1$parent == v2$parent)){
            if (ignore.seeding){
                res = 1
            }else{
                # check if the samples with non-zero cell frac at each
                # node are the same, if not, trees are different, so
                # this preseves the seeding models
                cat('\n***********Checking seeding models...\n')
                if (all(v1$samples.with.nonzero.cell.frac ==
                        v2$samples.with.nonzero.cell.frac)){
                    res = 1
                }
            }
        }
    }


    return(res)
}






find.matched.models <- function(vv, samples, sample.groups=NULL,
                                merge.similar.samples=FALSE){
    cat('Finding consensus models across samples...\n')
    nSamples = length(samples)
    matched = NULL
    scores = NULL
    # for historical reason, variables were named prim, met, etc., but
    # it does not mean samples are prim, met.
    find.next.match <- function(prim.idx, prim.model.idx,
                                met.idx, met.model.idx,
                                matched.models, model.scores){
        if (met.idx > nSamples){
            if (all(matched.models > 0)){
                matched <<- rbind(matched, matched.models)
                scores <<- rbind(scores, model.scores)
            }
        }else{
            for (j in 1:length(matched.models)){
                match.with.all.models = TRUE
                if (matched.models[j] > 0){
                    if(!match.sample.clones(vv[[j]][[matched.models[j]]],
                                            vv[[met.idx]][[met.model.idx]])){
                        match.with.all.models = FALSE
                        break
                    }
                }
            }
            #debug
            #cat(paste(matched.models[matched.models>0], collapse='-'),
            #          met.model.idx, match.with.all.models, '\n')
            if (match.with.all.models){
                matched.models[[met.idx]] = met.model.idx
                model.scores[[met.idx]] =
                    get.model.score(vv[[met.idx]][[met.model.idx]])
                find.next.match(prim.idx, prim.model.idx, met.idx+1, 1,
                                matched.models, model.scores)
            }#else{
            if (length(vv[[met.idx]]) > met.model.idx){
                matched.models[[met.idx]] = 0
                model.scores[[met.idx]] = 0
                find.next.match(prim.idx, prim.model.idx,
                                met.idx, met.model.idx+1,
                                matched.models, model.scores)
            }
            #find.next.match(prim.idx)
            #}

        }
    }
    for (prim.model in 1:length(vv[[1]])){
        # cat('prim.model =', prim.model, '\n')
        matched.models = c(prim.model,rep(0,nSamples-1))
        model.scores = c(get.model.score(vv[[1]][[prim.model]]),
                         rep(0,nSamples-1))
        find.next.match(1, prim.model, 2, 1, matched.models, model.scores)
    }
    num.models.found = ifelse(is.null(matched), 0, nrow(matched))
    cat('Found ', num.models.found, 'consensus model(s)\n')

	# merge clonal trees
    merged.trees = list()
    merged.traces = list()
    if (num.models.found > 0){
        cat('Generating consensus clonal evolution trees across samples...\n')
        for (i in 1:num.models.found){
            m = list()
            for (j in 1:nSamples){
                m = c(m, list(vv[[j]][[matched[i, j]]]))
            }

            zz = merge.clone.trees(m, samples=samples, sample.groups, merge.similar.samples=merge.similar.samples)
            mt = zz$merged.tree
            trace = zz$merged.trace
            # after merged, assign sample.group and color to individual tree
            #print(mt)
            for (j in 1:nSamples){
                vv[[j]][[matched[i, j]]] = merge(vv[[j]][[matched[i, j]]],
                    mt[, c('lab', 'sample.group', 'sample.group.color')],
                    all.x=TRUE)
            }

            merged.trees = c(merged.trees, list(mt))
            merged.traces = c(merged.traces, list(trace))
        }
    }

    return(list(models=vv, matched.models=matched, merged.trees=merged.trees,
        merged.traces=merged.traces, scores=scores))
}



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
        }else if (subclonal.test == 'bootstrap'){
            if (is.null(boot)){
                #boot = generate.boot(variants, vaf.col.names=vaf.col.names,
                #                     vaf.in.percent=vaf.in.percent,
                #                     num.boots=num.boots)

                boot = generate.boot(variants, vaf.col.names=vaf.col.names,
                                     depth.col.names=depth.col.names,
                                     vaf.in.percent=vaf.in.percent,
                                     num.boots=num.boots,
                                     bootstrap.model=subclonal.test.model,
                                     cluster.center.method=cluster.center,
                                     weighted=weighted,
                                     random.seed=random.seed)
                #bbb <<- boot
            }

            models = enumerate.clones(v, sample=s, variants, boot=boot,
                                      founding.cluster=founding.cluster,
                                      ignore.clusters=ignore.clusters,
                                      min.cluster.vaf=min.cluster.vaf,
                                      p.value.cutoff=sum.p.cutoff,
                                      alpha=alpha)
        }

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
    if (!is.null(matched)){
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


scale.cell.frac <- function(m, ignore.clusters=NULL){
    #max.vaf = max(m$vaf[!m$excluded & !(m$lab %in% as.character(ignore.clusters))])
    max.vaf = m$vaf[!is.na(m$parent) & m$parent=='-1']
    scale = 0.5/max.vaf
    m$vaf = m$vaf*scale
    m$free = m$free*scale
    m$free.mean = m$free.mean*scale
    m$free.lower = m$free.lower*scale
    m$free.upper = m$free.upper*scale
    m$occupied = m$occupied*scale
    return(m)
}


scale.sample.position <- function(xstarts, xstops, plot.total.length=7,
                                    evenly.distribute=TRUE){
    #print(xstarts)
    #print(xstops)
    if (any(xstarts >= xstops)){
        stop('\nERROR: stop positions must be greater than start positions\n')
    }
    # make sure the left-most sample starts an zero
    minx = min(xstarts)
    xstarts = xstarts - minx
    xstops = xstops - minx
    # evenly distribute
    if (evenly.distribute){
        #print(xstarts); print(xstops)
        k = length(xstarts)
        x = data.frame(x=c(xstarts, xstops),index=seq(1,2*k))
        x = x[order(x$x),]
        n = length(unique(x$x))
        x$even.x = 0
        even.x = 0:(n-1)
        xi1 = -1000
        j = 0
        for (i in 1:nrow(x)){
           if (x$x[i] != xi1){
              j = j + 1
              xi1 = x$x[i]
           }
           x$even.x[i] = even.x[j]
        }
        x = x[order(x$index),]
        xstarts = x$even.x[1:k]
        xstops = x$even.x[(k+1):(2*k)]
        #print(xstarts); print(xstops)
    }
    # scale position
    requested.total.length = max(xstops) - min(xstarts)
    xscale = plot.total.length/requested.total.length
    xstarts = xscale*xstarts
    xstops = xscale*xstops
    xlens = xstops - xstarts
    return(list(xstarts=xstarts, xstops=xstops, xlens=xlens))
}


plot.clonal.models <- function(y, out.dir,
                               models=NULL,
                               matched=NULL,
                               samples=NULL, #samples and vaf.col.names now are the same
                               variants=NULL,
                               variants.with.mapped.events=NULL,
                               models.to.plot=NULL,
                               clone.shape='bell',
                               bell.curve.step=0.25,
                               bell.border.width=1,
                               disable.sample.label=FALSE,
                               clone.time.step.scale=1,
                               zero.cell.frac.clone.color=NULL,
                               zero.cell.frac.clone.border.color=NULL,
                               zero.cell.frac.clone.border.width=0.25,
                               nonzero.cell.frac.clone.border.color=NULL,
                               nonzero.cell.frac.clone.border.width=0.25,

                               box.plot=FALSE,
                               cn.col.names=c(),
                               loh.col.names=c(),
                               mapped.events=NULL,
                               box.plot.text.size=1.5,
                               fancy.boxplot=FALSE,
                               fancy.boxplot.highlight=NULL,
                               event.col.name = 'event',
                               fancy.variant.boxplot.vaf.suffix='.VAF',
                               fancy.variant.boxplot.vaf.limits=70,
                               fancy.variant.boxplot.show.cluster.axis.label=FALSE,
                               fancy.variant.boxplot.sample.title.size=8,
                               fancy.variant.boxplot.panel.border.linetypes='solid',
                               fancy.variant.boxplot.panel.border.sizes=0.5,
                               fancy.variant.boxplot.panel.border.colors='black',
                               fancy.variant.boxplot.base_size=8,
                               fancy.variant.boxplot.axis.ticks.length=1,
                               fancy.variant.boxplot.axis.text.angle=0,
                               fancy.variant.boxplot.plot.margin=0.1,
                               fancy.variant.boxplot.horizontal=FALSE,
                               fancy.variant.boxplot.box=FALSE,
                               fancy.variant.boxplot.box.line.type='solid',
                               fancy.variant.boxplot.box.line.size=0.5,
                               fancy.variant.boxplot.box.outlier.shape=1,
                               fancy.variant.boxplot.box.alpha=0.5,
                               fancy.variant.boxplot.violin=FALSE,
                               fancy.variant.boxplot.violin.line.type='dotted',
                               fancy.variant.boxplot.violin.line.size=0.5,
                               fancy.variant.boxplot.violin.fill.color='grey80',
                               fancy.variant.boxplot.violin.alpha=0.5,
                               fancy.variant.boxplot.jitter=TRUE,
                               fancy.variant.boxplot.jitter.width=0.5,
                               fancy.variant.boxplot.jitter.color=NULL,
                               fancy.variant.boxplot.jitter.alpha=0.5,
                               fancy.variant.boxplot.jitter.size=1,
                               fancy.variant.boxplot.jitter.shape=1,
                               fancy.variant.boxplot.jitter.center.method='median',
                               fancy.variant.boxplot.jitter.center.color='black',
                               fancy.variant.boxplot.jitter.center.size=1,
                               fancy.variant.boxplot.jitter.center.linetype='solid',
                               fancy.variant.boxplot.jitter.center.display.value='none',# 'mean', 'median'
                               fancy.variant.boxplot.jitter.center.display.value.text.size=5,
                               fancy.variant.boxplot.highlight=NULL,
                               fancy.variant.boxplot.highlight.color='darkgray',
                               fancy.variant.boxplot.highlight.fill.color='red',
                               fancy.variant.boxplot.highlight.shape=21,
                               fancy.variant.boxplot.highlight.size=1,
                               fancy.variant.boxplot.highlight.color.col.name=NULL,
                               fancy.variant.boxplot.highlight.fill.col.names=NULL,
                               fancy.variant.boxplot.highlight.fill.min=1,
                               fancy.variant.boxplot.highlight.fill.mid=2,
                               fancy.variant.boxplot.highlight.fill.max=3,
                               fancy.variant.boxplot.highlight.fill.min.color='green',
                               fancy.variant.boxplot.highlight.fill.mid.color='black',
                               fancy.variant.boxplot.highlight.fill.max.color='red',
                               fancy.variant.boxplot.highlight.size.names=NULL,
                               fancy.variant.boxplot.max.highlight.size.value=500,
                               fancy.variant.boxplot.size.breaks=c(0,50,100,200,500),
                               fancy.variant.boxplot.highlight.size.legend.title='depth',
                               fancy.variant.boxplot.highlight.note.col.name=NULL,
                               fancy.variant.boxplot.highlight.note.color='blue',
                               fancy.variant.boxplot.highlight.note.size=1,
                               fancy.variant.boxplot.order.by.total.vaf=FALSE,
                               fancy.variant.boxplot.ccf=FALSE,
                               fancy.variant.boxplot.founding.cluster='1',
                               fancy.variant.boxplot.show.cluster.label=TRUE,


                               cluster.col.name = 'cluster',
                               ignore.clusters=NULL,# this param is now deprecated
                               scale.monoclonal.cell.frac=TRUE,
                               adjust.clone.height=TRUE,
                               individual.sample.tree.plot=FALSE,
                               merged.tree.plot=TRUE,
                               merged.tree.clone.as.branch=TRUE,
                               merged.tree.distance.from.bottom=0.01,#in
                               mtcab.tree.rotation.angle=180,
                               mtcab.tree.text.angle=NULL,
                               mtcab.tree.label=NULL,
                               mtcab.branch.angle=15, #mtcab=merged.tree.clone.as.branch
                               mtcab.branch.width=1,
                               mtcab.branch.text.size=0.3,
                               mtcab.node.size=3,
                               mtcab.node.label.size=0.75,
                               mtcab.node.text.size=0.5,
                               mtcab.event.sep.char=',',
                               mtcab.show.event=TRUE,
                               merged.tree.node.annotation='sample.with.nonzero.cell.frac.ci',
                               merged.tree.cell.frac.ci=FALSE,
                               trimmed.merged.tree.plot=TRUE,
                               trimmed.merged.tree.plot.width=NULL,
                               trimmed.merged.tree.plot.height=NULL,
                               tree.node.label.split.character=',',
                               tree.node.num.samples.per.line=NULL,
                               color.node.by.sample.group=FALSE,
                               color.border.by.sample.group=TRUE,
                               variants.to.highlight=NULL,
                               bell.event=FALSE,
                               bell.event.label.color='blue',
                               bell.event.label.angle=NULL,
                               width=NULL, height=NULL, text.size=1,
                               panel.widths=NULL,
                               panel.heights=NULL,
                               tree.node.shape='circle',
                               tree.node.size = 50,
                               tree.node.text.size=1,
                               merged.tree.node.size.scale=0.5,
                               merged.tree.node.text.size.scale=1,
                               out.format='png', resolution=300,
                               overwrite.output=FALSE,
                               max.num.models.to.plot=10,
                               cell.frac.ci=TRUE,
                               cell.frac.top.out.space=0.75,
                               cell.frac.side.arrow.width=1.5,
                               disable.cell.frac=FALSE,
                               show.score=TRUE,
                               show.matched.index=FALSE,
                               show.time.axis=TRUE,
                               xstarts=NULL,
                               xstops=NULL,

                               cell.plot = FALSE,
                               num.cells = 100,
                               cell.layout = 'cloud',
                               cell.border.size=0.1,
                               cell.border.color='black',
                               cell.size = 2,
                               cell.show.frame=FALSE,
                               clone.grouping='random',

                               out.prefix='model')
{
    if (!file.exists(out.dir)){
        dir.create(out.dir)
    }else{
        if (!overwrite.output){
            stop(paste('ERROR: Output directory (', out.dir,
                       ') exists. Quit!\n'))
        }
    }

    library(grid)
    library(gridExtra)

    # backward compatible
    x = y
    models = x$models
    matched = x$matched

    if (is.null(samples)){
        samples = names(models)
    }
    nSamples = length(samples)
    w = ifelse(is.null(width), 7, width)
    h = ifelse(is.null(height), 3*nSamples, height)
    w2h.scale <<- h/w/nSamples*ifelse(box.plot, 2, 1.5)
    # prepare sample positions
    if (is.null(xstarts) || is.null(xstops)){
        xstarts = rep(0, nSamples)
        xstops = rep(1, nSamples)
        names(xstarts) = names(xstops) = samples
    }else{
        xstarts = xstarts[samples]
        xstops = xstops[samples]
    }
    plot.total.length = 7
    if (disable.cell.frac){ plot.total.length = 9 }
    sample.pos = scale.sample.position(xstarts, xstops,
            plot.total.length=plot.total.length)
    #print(sample.pos)

    # debug
    #print(sample.pos)

    # backward compatible
    if (is.null(variants)){variants=x$variants}

    if(box.plot && is.null(variants)){
        box.plot = FALSE
        message('box.plot = TRUE, but variants = NULL. No box plot!')
    }
    if (!is.null(matched$index)){
        scores = matched$scores
        probs = matched$probs
        merged.trees = matched$merged.trees
        merged.traces = matched$merged.traces
        trimmed.trees = matched$trimmed.merged.trees
        # for historical reason, use 'matched' variable to indicate index of matches here
        matched = matched$index
        num.models = nrow(matched)
        if (num.models > max.num.models.to.plot &&
                !is.null(max.num.models.to.plot)){
            message(paste0(num.models,
               ' models requested to plot. Only plot the first ',
               max.num.models.to.plot,
               ' models. \nChange "max.num.models.to.plot" to plot more.\n'))
            matched = head(matched, n=max.num.models.to.plot)
        }

        # get driver events in tree if now driver event is provided for bell plot
        if (is.null(variants.to.highlight) && bell.event){
            mt = x$matched$merged.trees[[1]]
            if ('events' %in% colnames(mt)){
                mt = mt[!is.na(mt$events) & mt$events != '',]
                if (nrow(mt) > 0){
                    variants.to.highlight = mt[, c('lab', 'events')]
                    colnames(variants.to.highlight) = c('cluster', 'variant.name')
                }
            }
        }


        if (out.format == 'pdf'){
            pdf(paste0(out.dir, '/', out.prefix, '.pdf'), width=w, height=h,
                useDingbat=FALSE, title='')
        }

        for (i in 1:nrow(matched)){
            if (!is.null(models.to.plot)){
                if (!(i %in% models.to.plot)){next}
            }
            combined.graph = NULL
            this.out.prefix = paste0(out.dir, '/', out.prefix, '-', i)
            if (out.format == 'png'){
                png(paste0(this.out.prefix, '.png'), width=w,
                    height=h, res=resolution, units='in')
            }else if (out.format == 'pdf.multi.files'){
                pdf(paste0(this.out.prefix, '.pdf'), width=w, height=h,
                    useDingbat=FALSE, title='')
            }else if (out.format != 'pdf'){
                stop(paste0('ERROR: output format (', out.format,
                            ') not supported.\n'))
            }

            #num.plot.cols = ifelse(box.plot, 3, 2)
            #num.plot.cols = 2 + box.plot + merged.tree.plot
            num.plot.cols = 1 + box.plot + individual.sample.tree.plot + cell.plot
            par(mfrow=c(nSamples,num.plot.cols), mar=c(0,0,0,0))
            mat = t(matrix(seq(1, nSamples*num.plot.cols), ncol=nSamples))
            if (merged.tree.plot){mat = cbind(mat, rep(nSamples*num.plot.cols+1,nrow(mat)))}
            if (merged.tree.clone.as.branch){mat = cbind(mat, rep(nSamples*num.plot.cols+merged.tree.plot+1,nrow(mat)))}
            #print(mat)
            if (is.null(panel.widths)){
                ww = rep(1, num.plot.cols)
                if (merged.tree.plot){ww = c(ww , 1.5)}
                if (merged.tree.clone.as.branch){ww = c(ww , 0.5)}
                #ww[length(ww)] = 1
                if (box.plot){
                    ww[1] = 1
                }
            }else{
                if (length(panel.widths) != (num.plot.cols+merged.tree.plot+merged.tree.clone.as.branch)){
                    stop(paste0('ERROR: panel.widths length does not equal # of plots ',
                                num.plot.cols+merged.tree.plot+merged.tree.clone.as.branch, '\n'))
                }else{
                    ww = panel.widths
                }
            }

            hh = rep(1, nSamples)
            layout(mat, ww, hh)



            var.box.plots = NULL
            if (fancy.boxplot){
                if (is.null(variants.with.mapped.events)){
                    variants.with.mapped.events = variants
                }
                if (is.null(fancy.variant.boxplot.jitter.color)){
                    fancy.variant.boxplot.jitter.color = x$matched$merged.trees[[1]]$color
                    names(fancy.variant.boxplot.jitter.color) = x$matched$merged.trees[[1]]$lab
                    fancy.variant.boxplot.jitter.color = fancy.variant.boxplot.jitter.color[unique(variants.with.mapped.events$cluster)]
                    #fancy.variant.boxplot.jitter.color = get.clonevol.colors(length(unique(variants.with.mapped.events$cluster)))
                    #print(fancy.variant.boxplot.jitter.color)
                }
                var.box.plots = variant.box.plot(variants.with.mapped.events,
                    cluster.col.name=cluster.col.name,
                    vaf.col.names=samples,
                    show.cluster.size=FALSE,
                    variant.class.col.name=NULL,
                    vaf.suffix=fancy.variant.boxplot.vaf.suffix,
                    vaf.limits=fancy.variant.boxplot.vaf.limits,
                    show.cluster.axis.label=fancy.variant.boxplot.show.cluster.axis.label,
                    sample.title.size=fancy.variant.boxplot.sample.title.size,
                    panel.border.linetypes=fancy.variant.boxplot.panel.border.linetypes,
                    panel.border.sizes=fancy.variant.boxplot.panel.border.sizes,
                    panel.border.colors=fancy.variant.boxplot.panel.border.colors,
                    base_size=fancy.variant.boxplot.base_size,
                    axis.ticks.length=fancy.variant.boxplot.axis.ticks.length,
                    axis.text.angle=fancy.variant.boxplot.axis.text.angle,
                    plot.margin=fancy.variant.boxplot.plot.margin,
                    horizontal=fancy.variant.boxplot.horizontal,
                    box=fancy.variant.boxplot.box,
                    box.line.type=fancy.variant.boxplot.box.line.type,
                    box.line.size=fancy.variant.boxplot.box.line.size,
                    box.outlier.shape=fancy.variant.boxplot.box.outlier.shape,
                    box.alpha=fancy.variant.boxplot.box.alpha,
                    violin=fancy.variant.boxplot.violin,
                    violin.line.type=fancy.variant.boxplot.violin.line.type,
                    violin.line.size=fancy.variant.boxplot.violin.line.size,
                    violin.fill.color=fancy.variant.boxplot.violin.fill.color,
                    violin.alpha=fancy.variant.boxplot.violin.alpha,
                    jitter=fancy.variant.boxplot.jitter,
                    jitter.width=fancy.variant.boxplot.jitter.width,
                    jitter.color=fancy.variant.boxplot.jitter.color,
                    jitter.alpha=fancy.variant.boxplot.jitter.alpha,
                    jitter.size=fancy.variant.boxplot.jitter.size,
                    jitter.shape=fancy.variant.boxplot.jitter.shape,
                    jitter.center.method=fancy.variant.boxplot.jitter.center.method,
                    jitter.center.color=fancy.variant.boxplot.jitter.center.color,
                    jitter.center.size=fancy.variant.boxplot.jitter.center.size,
                    jitter.center.linetype=fancy.variant.boxplot.jitter.center.linetype,
                    jitter.center.display.value=fancy.variant.boxplot.jitter.center.display.value,
                    jitter.center.display.value.text.size=fancy.variant.boxplot.jitter.center.display.value.text.size,
                    highlight=fancy.variant.boxplot.highlight,
                    highlight.color=fancy.variant.boxplot.highlight.color,
                    highlight.fill.color=fancy.variant.boxplot.highlight.fill.color,
                    highlight.shape=fancy.variant.boxplot.highlight.shape,
                    highlight.size=fancy.variant.boxplot.highlight.size,
                    highlight.color.col.name=fancy.variant.boxplot.highlight.color.col.name,
                    highlight.fill.col.names=fancy.variant.boxplot.highlight.fill.col.names,
                    highlight.fill.min=fancy.variant.boxplot.highlight.fill.min,
                    highlight.fill.mid=fancy.variant.boxplot.highlight.fill.mid,
                    highlight.fill.max=fancy.variant.boxplot.highlight.fill.max,
                    highlight.fill.min.color=fancy.variant.boxplot.highlight.fill.min.color,
                    highlight.fill.mid.color=fancy.variant.boxplot.highlight.fill.mid.color,
                    highlight.fill.max.color=fancy.variant.boxplot.highlight.fill.max.color,
                    highlight.size.names=fancy.variant.boxplot.highlight.size.names,
                    max.highlight.size.value=fancy.variant.boxplot.max.highlight.size.value,
                    size.breaks=fancy.variant.boxplot.size.breaks,
                    highlight.size.legend.title=fancy.variant.boxplot.highlight.size.legend.title,
                    highlight.note.col.name=fancy.variant.boxplot.highlight.note.col.name,
                    highlight.note.color=fancy.variant.boxplot.highlight.note.color,
                    highlight.note.size=fancy.variant.boxplot.highlight.note.size,
                    display.plot=FALSE,
                    order.by.total.vaf=fancy.variant.boxplot.order.by.total.vaf,
                    ccf=fancy.variant.boxplot.ccf,
                    founding.cluster=fancy.variant.boxplot.founding.cluster,
                    show.cluster.label=fancy.variant.boxplot.show.cluster.label,
                    cluster.axis.name=''
                )
            }

            for (k in 1:length(samples)){
                s = samples[k]
                s.match.idx = matched[[s]][i]
                m = models[[s]][[matched[[s]][i]]]
                merged.tree = merged.trees[[i]]
                if (scale.monoclonal.cell.frac){
                    #TODO: auto identify ignore.clusters
                    m = scale.cell.frac(m, ignore.clusters=NULL)
                }
                lab = s
                # turn this on to keep track of what model matched
                if (show.matched.index){
                    lab = paste0(s, ' (', s.match.idx, ')')
                }
                if (show.score){
                    if (x$params$score.model.by == 'metap'){
                        lab = paste0(s, '\n(max.p=',
                                 sprintf('%0.3f', scores[[s]][i]), ')')
                    }else if(x$params$score.model.by == 'probability'){
                        lab = paste0(s, '\n(prob=',
                                 sprintf('%0.3f', probs[[s]][i]), ')')
                    }
                }
                top.title = NULL
                if (k == 1 && show.score){
                    if (x$params$score.model.by == 'metap'){
                        top.title = paste0('Max (clone cross-sample p) = ',
                            scores$max.clone.ccf.combined.p[i])
                    }else if(x$params$score.model.by == 'probability'){
                        top.title = paste0('Prob = ',
                                           probs$model.prob[i])
                    }
                }
                if (box.plot){
                    current.mar = par()$mar
                    par(mar=c(3,5,3,3))
                    if (fancy.boxplot){
                        library(grid)
                        library(gridBase)
                        plot.new()
                        vps = baseViewports()
                        pushViewport(vps$figure)
                        vp1 = plotViewport(c(0,0,0,0))
                        print(var.box.plots[[k]], vp=vp1)
                        popViewport()
                    }else{
                        with(variants, boxplot(get(s) ~ get(cluster.col.name),
                                           cex.lab=box.plot.text.size,
                                           cex.axis=box.plot.text.size,
                                           cex.main=box.plot.text.size,
                                           cex.sub=box.plot.text.size,
                                           ylab=s))
                    }

                    par(mar=current.mar)
                }
                bell.plot.center.to.top = 30
                if (disable.cell.frac){bell.plot.center.to.top=40}
                start.pos = 1
                if (disable.sample.label){start.pos=0.1}
                draw.sample.clones(m, x=start.pos+sample.pos$xstarts[k], y=0,
                                   wid=bell.plot.center.to.top,
                                   #len=7,
                                   len=sample.pos$xlens[k],
                                   wscale=7/w,
                                   clone.shape=clone.shape,
                                   bell.curve.step=bell.curve.step,
                                   clone.time.step.scale=clone.time.step.scale,
                                   bell.border.width=bell.border.width,
                                   zero.cell.frac.clone.color=zero.cell.frac.clone.color,
                                   zero.cell.frac.clone.border.color=zero.cell.frac.clone.border.color,
                                   zero.cell.frac.clone.border.width=zero.cell.frac.clone.border.width,
                                   nonzero.cell.frac.clone.border.color=nonzero.cell.frac.clone.border.color,
                                   nonzero.cell.frac.clone.border.width=nonzero.cell.frac.clone.border.width,
                                   label=lab,
                                   text.size=text.size,
                                   cell.frac.ci=cell.frac.ci,
                                   disable.cell.frac=disable.cell.frac,
                                   top.title=top.title,
                                   adjust.clone.height=adjust.clone.height,
                                   cell.frac.top.out.space=cell.frac.top.out.space,
                                   cell.frac.side.arrow.width=cell.frac.side.arrow.width,
                                   variants.to.highlight=variants.to.highlight,
                                   variant.color=bell.event.label.color,
                                   variant.angle=bell.event.label.angle,
                                   show.time.axis=show.time.axis,
                                   color.node.by.sample.group=color.node.by.sample.group,
                                   color.border.by.sample.group=color.border.by.sample.group)

                if (cell.plot){
                    current.mar = par()$mar
                    par(mar=c(0.1,0.1,0.1,0.1))
                    mx = m[(!m$excluded & !m$is.zero),]
                    cp = plot.cell.population(mx$free.mean/sum(mx$free.mean),
                        mx$color, layout=cell.layout,
                        cell.border.size=cell.border.size, cell.border.color=cell.border.color,
                        clone.grouping=clone.grouping,
                        num.cells=num.cells,
                        frame=cell.show.frame)
                    library(grid)
                    library(gridBase)
                    plot.new()
                    vps2 = baseViewports()
                    pushViewport(vps2$figure)
                    vp2 = plotViewport(c(0,0,0,0))
                    print(cp, vp=vp2)
                    popViewport()
                    par(mar=current.mar)

                }


                if (individual.sample.tree.plot){
                    gs = plot.tree(m, node.shape=tree.node.shape,
                               node.size=tree.node.size,
                               tree.node.text.size=tree.node.text.size,
                               cell.frac.ci=cell.frac.ci,
                               color.node.by.sample.group=color.node.by.sample.group,
                               color.border.by.sample.group=color.border.by.sample.group,
                               show.legend=FALSE,
                               node.prefix.to.add=paste0(s,': '),
                               out.prefix=paste0(this.out.prefix, '__', s))
                }


                # plot merged tree
                if (merged.tree.plot && k == nSamples){
                    current.mar = par()$mar
                    #par(mar=c(3,3,3,3))
                    #par(mai=c(merged.tree.distance.from.bottom,0.01,0.01,0.01))
                    par(mar=c(0,0,0,0))

                    #if (merged.tree.clone.as.branch){
                    #    plot.tree.clone.as.branch(merged.tree,
                    #        tree.rotation=mtcab.tree.rotation.angle,
                    #        text.angle=mtcab.tree.text.angle,
                    #        angle=mtcab.branch.angle,
                    #        branch.width=mtcab.branch.width,
                    #        branch.text.size=mtcab.branch.text.size,
                    #        node.size=mtcab.node.size,
                    #        node.label.size=mtcab.node.label.size,
                    #        node.text.size=mtcab.node.text.size,
                    #        event.sep.char=mtcab.event.sep.char,
                    #        tree.label = mtcab.tree.label,
                    #        show.event = mtcab.show.event
                    #    )
                    #}else{

                        # determine colors based on sample grouping
                        node.colors = NULL
                        if ('sample.group' %in% colnames(merged.tree)){
                            node.colors = merged.tree$sample.group.color
                            names(node.colors) = merged.tree$lab
                        }

                        gs2 = plot.tree(merged.tree,
                                   node.shape=tree.node.shape,
                                   node.size=tree.node.size*merged.tree.node.size.scale,
                                   tree.node.text.size=tree.node.text.size*merged.tree.node.text.size.scale,
                                   node.annotation=merged.tree.node.annotation,
                                   node.label.split.character=tree.node.label.split.character,
                                   node.num.samples.per.line=tree.node.num.samples.per.line,
                                   cell.frac.ci=merged.tree.cell.frac.ci,
                                   #title='\n\n\n\n\n\nmerged\nclonal evolution\ntree\n|\n|\nv',
                                   node.prefix.to.add=paste0(s,': '),
                                   #node.colors=node.colors,
                                   color.node.by.sample.group=color.node.by.sample.group,
                                   color.border.by.sample.group=color.border.by.sample.group,
                                   out.prefix=paste0(this.out.prefix, '__merged.tree__', s))
                    #}
                    par(mar=current.mar)
                }

                if (merged.tree.clone.as.branch && k == nSamples){
                    current.mar = par()$mar
                    #par(mar=c(3,3,3,3))
                    #par(mai=c(merged.tree.distance.from.bottom,0.01,0.01,0.01))
                    par(mar=c(0,0,0,0))

                    plot.tree.clone.as.branch(merged.tree,
                        tree.rotation=mtcab.tree.rotation.angle,
                        text.angle=mtcab.tree.text.angle,
                        angle=mtcab.branch.angle,
                        branch.width=mtcab.branch.width,
                        branch.text.size=mtcab.branch.text.size,
                        node.size=mtcab.node.size,
                        node.label.size=mtcab.node.label.size,
                        node.text.size=mtcab.node.text.size,
                        event.sep.char=mtcab.event.sep.char,
                        tree.label = mtcab.tree.label,
                        show.event = mtcab.show.event
                    )

                    par(mar=current.mar)

                }

                if (individual.sample.tree.plot){
                    if (is.null(combined.graph)){
                        combined.graph = gs
                    }else{
                        combined.graph = graph.union(combined.graph, gs,
                                                 byname=TRUE)
                        # set color for all clones, if missing in 1st sample
                        # get color in other sample
                        V(combined.graph)$color <-
                            ifelse(is.na(V(combined.graph)$color_1),
                                   V(combined.graph)$color_2,
                                   V(combined.graph)$color_1)
                   }
                }
            }
            if (out.format == 'png' || out.format == 'pdf.multi.files'){
                dev.off()
            }
            if (individual.sample.tree.plot){
                write.graph(combined.graph,
                        file=paste0(this.out.prefix, '.graphml'),
                        format='graphml')
            }
        }
        if (out.format == 'pdf'){
            #plot(combined.graph)
            dev.off()
        }

        # plot trimmed trees
        if (nrow(trimmed.trees[[1]]) == 0){#single sample, no trimmed merged tree
            trimmed.merged.tree.plot = FALSE
        }
        if (trimmed.merged.tree.plot){
            cat('Plotting pruned consensus trees...\n')
            if (is.null(trimmed.merged.tree.plot.width) ||
                is.null(trimmed.merged.tree.plot.height)){
                    trimmed.merged.tree.plot.width = w/num.plot.cols*2
                    trimmed.merged.tree.plot.height = h/nSamples*7
            }
            pdf(paste0(out.dir, '/', out.prefix, '.pruned-trees.pdf'),
                width=trimmed.merged.tree.plot.width,
                height=trimmed.merged.tree.plot.height,
                useDingbat=FALSE, title='')
            for (i in 1:length(trimmed.trees)){
                gs3 = plot.tree(trimmed.trees[[i]],
                           node.shape=tree.node.shape,
                           node.size=tree.node.size*0.5,
                           tree.node.text.size=tree.node.text.size,
                           node.annotation=merged.tree.node.annotation,
                           node.label.split.character=tree.node.label.split.character,
                           node.num.samples.per.line=tree.node.num.samples.per.line,
                           color.border.by.sample.group=color.border.by.sample.group,
                           #cell.frac.ci=cell.frac.ci,
                           cell.frac.ci=FALSE,
                           node.prefix.to.add=paste0(s,': '),
                           out.prefix=paste0(this.out.prefix,
                                             '__trimmed.merged.tree__', s))

            }
            dev.off()

        }

        # write trace file, telling what samples are eliminated from tree
        # due to its similar clonal architecture to another sample
        if (!is.null(merged.traces[[1]])){
            all.traces = NULL
            for (i in 1:nrow(matched)){
                tmp = merged.traces[[i]]
                tmp = cbind(model=i, tmp)
                if (is.null(all.traces)){all.traces = tmp}
                else{all.traces = rbind(all.traces, tmp)}
            }
            all.traces$sample = factor(all.traces$sample,
                                       levels=unique(all.traces$sample))
            all.traces = aggregate(similar.sample ~ model+sample, all.traces,
                paste, collapse=',')
            all.traces = all.traces[order(all.traces$model),]
            write.table(all.traces, paste0(out.dir, '/', out.prefix,
                '.sample-reduction.tsv'), sep='\t', row.names=FALSE,
                quote=FALSE)
        }

    }else{# of !is.null(matched$index); plot all
        # plot all models for all samples separately.
        # This will serve as a debug tool for end-user when their models
        # from different samples do not match.
        message('No consensus multi-sample models provided.
                Individual sample models will be plotted!\n')
        for (s in names(models)){
            draw.sample.clones.all(models[[s]],
                                   paste0(out.dir, '/', out.prefix, '-', s))
        }

    }
    cat(paste0('Output plots are in: ', out.dir, '\n'))

}



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

adjust.clone.vaf <- function(clone.vafs, var, cluster.col.name,
                             founding.cluster=1,
                             adjust.to.founding.cluster.only=TRUE,
                             p.value.cut=0.01){
    vaf.names = colnames(clone.vafs[2:length(colnames(clone.vafs))])
    founding.cluster.idx = which(clone.vafs$cluster == founding.cluster)
    base.clusters.idx = unique(c(founding.cluster.idx, 1:(nrow(clone.vafs)-1)))
    if (adjust.to.founding.cluster.only){
        base.clusters.idx = founding.cluster.idx
    }
    #debug
    #print(base.clusters.idx)
    for (vaf.name in vaf.names){
        #for (i in 1:(nrow(clone.vafs)-1)){
        for (i in base.clusters.idx){
            ci = clone.vafs$cluster[i]
            vaf.i = clone.vafs[clone.vafs[[cluster.col.name]]==ci, vaf.name]
            for (j in (i+1):nrow(clone.vafs)){
                cj = clone.vafs$cluster[j]
                vaf.j = clone.vafs[clone.vafs[[cluster.col.name]]==cj, vaf.name]
                if (!clone.vaf.diff(var[var[[cluster.col.name]]==ci,vaf.name],
                                    var[var[[cluster.col.name]]==cj,vaf.name])){

                    clone.vafs[clone.vafs[[cluster.col.name]]==cj, vaf.name] =
                        vaf.i
                }
            }
        }
    }
    return(clone.vafs)
}


clone.vaf.diff <- function(clone1.vafs, clone2.vafs, p.value.cut=0.05){
    tst = wilcox.test(clone1.vafs, clone2.vafs)
    #print(tst)
    if (is.na(tst$p.value) || tst$p.value <= p.value.cut){
        return(TRUE)
    }else{
        return(FALSE)
    }
}



is.poly <- function(v){
    res = FALSE
    if (nrow(v[!is.na(v$parent) & v$parent == '0' & !v$excluded,]) > 1){
        res = TRUE
    }
    return(res)
}


sum.polyclonal <- function(x){
    samples = colnames(x$matched$index)
    poly = matrix(rep(NA, length(samples)), nrow=1)
    colnames(poly) = samples
    for (i in 1:nrow(x$matched$index)){
        p = c()
        for (s in samples){
           p = c(p, is.poly(x$models[[s]][[x$matched$index[i, s]]]))
        }
        poly = rbind(poly, p)
    }
    poly = poly[!is.na(poly[,1]),]
    rownames(poly) = seq(1, nrow(poly))
    poly = as.data.frame.matrix(poly)
    return(poly)
}


merge.samples <- function(x, samples, new.sample, new.sample.group, ref.cols=NULL, var.cols=NULL){
    if (!all(samples %in% names(x$models))){
        stop('ERROR: Sample not found when merging!')
    }

    # merge trees from samples (we need to do this to make sure only clones
    # from the samples' trees are included
    x$models[[new.sample]] = list()
    for (mid in 1:x$num.matched.models){
        #print(mid)
        v = NULL
        for (i in 1:length(samples)){
            s = samples[i]
            vi = x$models[[s]][[x$matched$index[mid,s]]]
            rownames(vi) = vi$lab
            #print(vi[, c( 'lab', 'excluded', 'parent', 'free.mean')])
            if (is.null(v)){
                v = vi
            }else{
                if (nrow(v) != nrow(vi)){
                    stop('ERROR: Number of clusters/clones not equal while merging\n')
                }
                vi = vi[as.character(v$lab),]
                v$vaf = v$vaf + vi$vaf
                vi.only = !vi$excluded & v$excluded
                v$parent[vi.only] = vi$parent[vi.only]
                v$excluded[vi.only] = FALSE
                v$ancestors[vi.only] = vi$ancestors[vi.only]
            }
        }
        v$vaf = v$vaf/length(samples)
        rownames(v) = v$lab
        #print(v[, c( 'lab', 'excluded', 'parent', 'free.mean')])
        x$models[[new.sample]][[mid]] = v
    }

    # update model index in matched
    x$matched$index[[new.sample]] = seq(1,x$num.matched.models)

    # recalculate VAF using read counts combined from all samples
    # if ref and var counts available
    if (!is.null(ref.cols) && !is.null(var.cols)){
        new.ref.col = paste0(new.sample, '.ref')
        new.var.col = paste0(new.sample, '.var')
        new.depth.col = paste0(new.sample, '.depth')

        x$variants[[new.ref.col]] = rowSums(x$variants[, ref.cols], na.rm=TRUE)
        x$variants[[new.var.col]] = rowSums(x$variants[, var.cols], na.rm=TRUE)
        x$variants[[new.depth.col]] = x$variants[[new.ref.col]] +
                                        x$variants[[new.var.col]]

        # recalc vaf of variants
        x$variants[[new.sample]] = x$variants[[new.var.col]]/x$variants[[new.depth.col]]
        if (x$params$vaf.in.percent){
            x$variants[[new.sample]] = 100*x$variants[[new.sample]]
        }
    }

    # recalc center vaf of clusters
    c = estimate.clone.vaf(x$variants, x$params$cluster.col.name,
                            vaf.col.names=new.sample,
                            vaf.in.percent=x$params$vaf.in.percent,
                            method=x$params$cluster.center.method)
    tmp = make.clonal.data.frame(c[[new.sample]], c[[x$params$cluster.col.name]],
        colors=unique(x$models[[1]][[1]]$color))
    rownames(tmp) = tmp$lab
    for (mid in 1:x$num.matched.models){
        tmp = tmp[as.character(x$models[[new.sample]][[mid]]$lab),]
        x$models[[new.sample]][[mid]]$vaf = tmp$vaf
        x$models[[new.sample]][[mid]]$vaf = tmp$free
        x$models[[new.sample]][[mid]]$free.mean = 0
    }

    # update bootstrap samples for the merged sample
    boot = generate.boot(x$variants, vaf.col.names=new.sample,
                        vaf.in.percent=x$params$vaf.in.percent,
                        num.boots=x$params$num.boots,
                        bootstrap.model=x$params$bootstrap.model,
                        cluster.center.method=x$params$cluster.center.method,
                        random.seed=x$params$random.seed)
    # add new sample' bootstraps, and push zero.means down to bottom in boot list
    x$boot[[new.sample]] = boot[[new.sample]]
    tmp = x$boot[['zero.means']]
    x$boot[['zero.means']] = NULL
    x$boot[['zero.means']] = tmp
    tmp = NULL

    # reestimate CCF
    cat('Estimating CCF of clones for merged sample...\n')
    for (mid in 1:x$num.matched.models){
        v =  x$models[[new.sample]][[mid]]
        rownames(v) = v$lab
        for (cl in v$lab[!v$excluded]){
            sub.clusters = v$lab[v$parent == cl & !is.na(v$parent)]
            if(length(sub.clusters) == 0){sub.clusters = NULL}
            v = estimate.ccf(v, new.sample, which(v$lab == cl), x$boot,
                x$params$min.cluster.vaf, x$params$alpha,
                t=NULL, sub.clusters=sub.clusters)
            #cat('ccf: ', v[cl,'lab'], paste(sub.clusters, collapse=','), '\n')
        }
        clone.stat = determine.subclone(v, v$lab[!is.na(v$parent)
                                         & v$parent == '-1'])

        v$is.subclone = clone.stat$is.sub[v$lab]
        v$is.founder = clone.stat$is.founder[v$lab]
        v$is.zero = clone.stat$is.zero
        x$models[[new.sample]][[mid]] = v
    }

    # add sample group
    x$params$sample.groups[new.sample] = new.sample.group

    # cleanup: remove models of samples merged to ensure data structure
    # consistency
    x$params$vaf.col.names = c(setdiff(x$params$vaf.col.names, samples), new.sample)
    for (s in samples){
        x$params$sample.groups = x$params$sample.groups[names(x$params$sample.groups) != s]
        x$models[[s]] = NULL
        x$matched$index[[s]] = NULL
    }

    x = merge.all.matched.clone.trees(x)

    return(x)
}

merge.all.matched.clone.trees <- function(x){
    if (x$num.matched.models == 0 || is.null(x$matched)){
        message('WARN: No matched model to merge.\n')
        return(x)
    }
    samples = names(x$params$sample.groups)
    merged.trees = list()
    merged.traces = list()
    cat('Merging clonal evolution trees across samples...\n')

    #????
    #x$params$sample.groups =

    for (i in 1:x$num.matched.models){
        m = list()
        for (j in 1:length(samples)){
            # the order of samples should be the same
            m = c(m, list(x$models[[j]][[x$matched$index[i, j]]]))
        }

        zz = merge.clone.trees(m, samples=samples, x$params$sample.groups,
               merge.similar.samples=x$params$merge.similar.samples)
        mt = zz$merged.tree
        trace = zz$merged.trace
        # after merged, assign sample.group and color to individual tree
        #print(mt)
        for (j in 1:length(samples)){
            x$models[[j]][[x$matched$index[i, j]]] = merge(x$models[[j]][[x$matched$index[i, j]]],
                mt[, c('lab', 'sample.group', 'sample.group.color')], all.x=TRUE)
        }

        # if events already mapped to old merged trees, take over from one of them
        if ('events' %in% colnames(x$matched$merged.trees[[1]])){
            mt = merge(mt, x$matched$merged.trees[[1]][, c('lab', 'events')])
        }

        merged.trees = c(merged.trees, list(mt))
        merged.traces = c(merged.traces, list(trace))
    }
    x$matched$merged.trees = merged.trees
    x$matched$merged.traces = merged.traces

    return(x)
}


assign.events.to.clones.of.a.tree <- function(tree, events, samples, cutoff=0){
    rownames(tree) = tree$lab
    if(nrow(tree) == 0 || nrow(events) == 0){return(NULL)}

    # strip off sample note (eg. zero cell frac)
    tree$samples = gsub('\u00B0|\\*', '', tree$sample)

    # make binary based on vaf cutoff
    events[, samples][events[, samples] < cutoff] = 0
    events[, samples][events[, samples] >= cutoff] = 1
    events = events[apply(events[, samples] > 0, 1, sum, na.rm=TRUE) > 0,]
    rownames(events) = NULL
    tree$events = ''

    # map events to clones
    # for each event, find the 1st clone that have the max ratio of
    # samples carrying the event
    for (i in 1:nrow(events)){# each event
        e = events[i,samples] > 0
        event.samples = colnames(e)[e[1,]]
        # find clone that shared the most number of samples
        # with the event
        best.match.idx = NULL
        best.match.clone = NULL
        max.match.rate = 0
        max.match.num.samples = 0
        for (j in 1:nrow(tree)){ # each clone
            clone = tree$lab[j]
            clone.samples = unlist(strsplit(tree$samples[j], ','))
            num.match.samples = sum(clone.samples %in% event.samples)
            match.rate = num.match.samples^2/length(clone.samples)
            #match.rate = num.match.samples*length(clone.samples)
            if (match.rate > max.match.rate ){
                max.match.rate = match.rate
                best.match.clone = clone
                best.match.idx = j
            }

        }
        tree$events[best.match.idx] = paste0(tree$events[best.match.idx],
                                                events$event[i], ',')
    }
    tree$events = gsub(',$', '', tree$events)
    tree$samples = NULL
    return(tree)

}


extract.mapped.events <- function(x){
    e = x$matched$merged.trees[[1]][, c('lab', 'events')]
    ee = NULL
    for (i in 1:nrow(e)){
        if (e[i,]$events == ''){next}
        evnts = unlist(strsplit(e[i,]$events, ','))
        ei = data.frame(cluster=e[i,]$lab, clone = e[i,]$lab,
            event=evnts, stringsAsFactors=FALSE)
        if (is.null(ee)){ee = ei}else{ee=rbind(ee,ei)}
    }
    return(ee)
}



assign.events.to.clones <- function(x, events, samples, cutoff=0){
    if (x$num.matched.models > 0){
        for (i in 1:x$num.matched.models){
            x$matched$merged.trees[[i]] = assign.events.to.clones.of.a.tree(
                x$matched$merged.trees[[i]], events, samples, cutoff)
        }
        x$events = merge(extract.mapped.events(x), events)
        # TODO: think about this.
        #x$variants.with.mapped.events = merge.variants.and.events(x$variants,
        #    x$events, vaf.col.names=samples, other.col.names=c())
    }
    return(x)
}


transfer.events.to.consensus.trees <- function(x, events,
        cluster.col.name='cluster',
        event.col.name){
    events$lab = as.character(events[[cluster.col.name]])
    events = events[, c(cluster.col.name, event.col.name)]
    colnames(events) = c('lab', 'events')
    events = aggregate(events ~ lab, data=events, paste, collapse=',')

    if (x$num.matched.models > 0){
        for (i in 1:x$num.matched.models){
            mt = x$matched$merged.trees[[i]]
            mt = merge(mt, events, all.x=TRUE)
            mt$events[is.na(mt$events)] = ''
            x$matched$merged.trees[[i]] = mt
        }
    }
    return(x)
}



get.clonevol.colors <- function(num.colors, strong.color=FALSE){
    colors = c('#cccccc',
               #'#b3b3b3',
               #'#999999',
               '#a6cee3', '#b2df8a', '#cab2d6','#ff99ff', '#fdbf6f', '#fb9a99',
               #'#bf812d',
               '#bbbb77',
               '#cf8d30',
               '#41ae76', '#ff7f00',
               '#3de4c5',
               #'#aaaa55',
               '#ff1aff',
               #'#c51b7d',
               '#9933ff', '#3690c0','#8c510a', '#666633', '#54278f',
               '#e6e600',
               '#e31a1c', '#00cc00', '#0000cc', '#252525',
               #'#ef6548',
               '#fccde5',  '#d9d9d9',
               #'#33a02c', '#3f007d', '#1f78b4',
               '#f0ecd7', '#ffffb3',
               #'#ffcc00',
               #'#fca27e', '#fb8072',
               '#ffff00', rep('#e5f5f9',10000))
    if(strong.color){
        colors = c('#e41a1c', '#377eb8', '#4daf4a', '#984ea3',
            '#ff7f00', 'black', 'darkgray', rep('lightgray',10000))
        colors[1:3] = c('red', 'blue', 'green')
    }
    if (num.colors > length(colors)){
        stop('ERROR: Not enough colors!\n')
    }else{
        return(colors[1:num.colors])
    }
}



plot.clonevol.colors <- function(num.colors=40){
    library(ggplot2)
    colors = get.clonevol.colors(num.colors)
    x = data.frame(hex=colors, val=1, stringsAsFactors=FALSE)
    x$hex = paste0(sprintf('%02d', seq(1,num.colors)), '\n', x$hex)
    names(colors) = x$hex
    print(x)
    p = (ggplot(x, aes(x=hex, y=val, fill=hex))
         + geom_bar(stat='identity')
         + theme_bw()
         + scale_fill_manual(values=colors)
         + theme(legend.position='none')
         + ylab(NULL) + xlab(NULL)
         + theme(axis.text.y=element_blank())
         + theme(axis.ticks.y=element_blank())
         + ggtitle('Clonevol colors'))
    ggsave(p, file='clonevol.colors.pdf', width=num.colors*0.75, height=4)
}



insert.lf <- function(ss, n, split.char=','){
        if (!is.null(split.char)){
            num.splits = sapply(ss, function(l)
                nchar(gsub(paste0('[^', split.char, ']'), '', l)))

            # only keep the node.label.split.char in interval of node.num.samples.per.line
            # such that a block of node.num.samples.per.line samples will be grouped and
            # kept in one line
            if (!is.null(n)){
                for (i in 1:length(ss)){
                    vl = unlist(strsplit(ss[i], split.char))
                    sel = seq(min(n, length(vl)),
                        length(vl),n)
                    vl[sel] = paste0(vl[sel], split.char)
                    vl[-sel] = paste0(vl[-sel], ';')
                    ss[i] = paste(vl, collapse='')
                }
                num.splits = length(sel) + 1
            }
            extra.lf = sapply(num.splits, function(n) paste(rep('\n', n), collapse=''))
            extra.lf = ''
            ss = paste0(extra.lf, gsub(split.char, '\n',ss))
         }
         return(ss)
}




plot.cell.population <- function(cell.frac, colors, labels=NULL,
    cell.cex=2, delta=NULL, cell.border.color='black', cell.border.size=0.1,
    num.cells=200, layout='cloud', clone.grouping='random', frame=FALSE){

    # generate approximately num.cells positions
    n = round(sqrt(num.cells))
    num.cells = n*n
    x = c(); y = c()
    xpos = n:1
    for (i in 1:n){
        xpos = rev(xpos) # rev to make continous x coord in plot
        x = c(x, xpos)
        y = c(y, rep(i,n))
    }

    # sort, round cell frac, calculate number of cells for each clone
    cell.frac = 100*cell.frac
    cell.frac[cell.frac < 0] = 0
    colors = colors[order(cell.frac, decreasing=TRUE)]
    cell.frac = sort(cell.frac, decreasing=TRUE)
    cell.frac = cell.frac*num.cells/100
    cells = round(cell.frac)

    # make sure num.cells total = num.cells (100%)
    cells[1] = cells[1] - (sum(cells) - num.cells)
    # generate color vector matching with number of cells
    cols = rep('black', num.cells)
    idx = 1
    for (i in 1:length(cells)){
        ncells = cells[i]
        if (ncells > 0){
            cols[idx:(idx+ncells-1)] = colors[i]
            idx = idx+ncells
        }
    }

    current.mar = par()$mar
    par(mar=c(0.1,0.1,0.1,0.1))

    if (layout == 'cloud'){
        cells = generate.cloud.of.cells(colors=cols)
        p = plot.cloud.of.cells(cells, frame=frame,
                cell.border.color=cell.border.color,
                cell.border.size=cell.border.size,
                clone.grouping=clone.grouping)
        return(p)
    }else if (layout == 'plate'){
        # plot cells using points
        if (is.null(delta)){delta = cell.cex/5}
        plot(x, y, col=cell.border.color, bg=cols, lwd=cell.border.size,
             axes=FALSE, cex=cell.cex,
            pch=21, xlim=c(1-delta,n+delta), ylim=c(1-delta,n+delta))
    }else{
        stop(paste0('ERROR: plot cell population layout=', layout,
            ' not supported.\n'))
    }

    par(mar=current.mar)
}



generate.cloud.of.cells <- function(colors, maxiter=1000){

    n = length(colors)
    limits = c(-50, 50)
    inset = diff(limits)/3
    # scale radius to make sure cloud of cells gather approximately like a sphere
    # with 200 cells, limits = c(-50,50), radius ~ 3.3 makes sure circles
    radius = sqrt(11*200/n)
    xyr = data.frame(
      x = runif(n, min(limits) + inset, max(limits) - inset),
      y = runif(n, min(limits) + inset, max(limits) - inset),
      r = rep(radius, n))

    # Next, we use the `circleLayout` function to try to find a non-overlapping
    # arrangement, allowing the circles to occupy any part of the bounding square.
    # The returned value is a list with elements for the layout and the number
    # of iterations performed.
    library(packcircles)
    res = circleLayout(xyr, limits, limits, maxiter = maxiter)

    ## plot data for the `after` layout returned by circleLayout
    cells = circlePlotData(res$layout)

    # color the circles
    cells$color = sample(colors, nrow(cells), replace=TRUE)

    return(cells)
}


plot.cloud.of.cells <- function(cells, title='', alpha=1, frame=FALSE,
    cell.border.color='black', cell.border.size=0.1,
    clone.grouping='random', limits=c(-50,50)){

    library(ggplot2)
    library(gridExtra)
    if (clone.grouping == 'horizontal'){
        cells$color[order(cells$y, cells$x)] = sort(cells$color)
    }else if (clone.grouping == 'vertical'){
        cells$color[order(cells$x, cells$y)] = sort(cells$color)
    }

    p = (ggplot(cells) +
        geom_polygon(aes(x, y, group=id, fill=color), color=cell.border.color,
                alpha=alpha, size=cell.border.size) +
        coord_equal(xlim=limits, ylim=limits) +
        theme_bw() +
        theme(axis.text=element_blank(),
            axis.ticks=element_blank(),
            axis.title=element_blank(),
            legend.position='none',
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank()) +
        theme(plot.margin=unit(c(0,0,0,0), 'mm')) +
        #labs(title=title) +
        #scale_fill_manual(values=colors) +
        scale_fill_identity()
    )
    if (!frame){
        p = p + theme(panel.border=element_blank())
    }else{
        p = p + theme(panel.border=element_rect(linetype='dotted'))
    }
    return(p)
}



