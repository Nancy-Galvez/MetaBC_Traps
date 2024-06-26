### Functions (and example) for performing the analyses described in:
### Multimodal inference for probabilistic species-delineation in the analysis of environmental DNA sequence data 
### by Jeff R. Powell (jeffpowell2@gmail.com)

### tested on 5 November 2010
#> sessionInfo()
#R version 2.10.1 (2009-12-14) 
#i386-apple-darwin9.8.0 
#
#locale:
#[1] en_US.UTF-8/en_US.UTF-8/C/C/en_US.UTF-8/en_US.UTF-8
#
#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base     
#
#other attached packages:
# [1] gtools_2.6.1  vegan_1.17-0  igraph_0.5.3  geiger_1.3-1  ouch_2.6-1    subplex_1.1-3 msm_0.9.5     mvtnorm_0.9-9 splits_1.0-11 paran_1.4.3   MASS_7.3-6    ape_2.5-3    
#
#loaded via a namespace (and not attached):
#[1] gee_4.13-14     grid_2.10.1     lattice_0.17-26 nlme_3.1-96     splines_2.10.1  survival_2.35-7 tools_2.10.1   
#> 
#


#### example ###
#
## have to load the functions before running the example
#
#source('Powell_multimodelInferenceGMYC.R')  # file must be in working directory, loads required packages
#
#
###fit gmyc models to test.tr
#
#data(test.tr)
#test.sing<-gmyc.edit(test.tr,method='s')  # fits single-threshold model
#test.mult<-gmyc.edit(test.tr,method='m')  # fits multiple-threshold model
#
#
###estimate model scores
#
#test.scores<-model.scores(test.sing,test.mult)
#summary.scores(test.scores,deltaAICc=3)
#
#
###estimate cluster probabilities that each pair of tips co-occurs within a cluster
#
#test.probs<-cluster.probs(test.scores,deltaAICc=0.999)  # will ask if you want to proceed, type 'y'
#image(test.probs$clusterMatrix.avg)
#plot.probs(test.probs)
#
#
###generate sample matrix for sample-specific estimates
#
#samp1<-rbinom(length(test.tr$tip.label),5,0.03)
#samp2<-rbinom(length(test.tr$tip.label),5,0.3)
#test.samples<-matrix(c(samp1,samp2),ncol=2,dimnames=list(test.tr$tip.label,c('samp1','samp2')))
#rm(samp1,samp2)
#
#
###estimate richness and diversity
#
#sample.est(test.scores,deltaAICc=0.999)  # assumes single sample, will ask if you want to proceed, type 'y' --> can compare to summary.scores(test.scores) to see effect on estimates of excluding models with low weights
#sample.est(test.scores,test.samples,deltaAICc=0.999)  ## sample-specific estimates, will ask if you want to proceed, type 'y'
#
#
###END example
#



### functions ###

library(splits)
library(geiger)
library(igraph)
library(vegan)
library(gtools)


# modified version of gmyc() for estimating speciation-coalescent transitions (changes are indicated with comments)
# MUST USE THIS VERSION FOR MULTIMODAL INFERENCE:
# estimates parameters assuming all isolates belong to a unique species (Yule process)
# also added diagnostic output indicating whether convergence was achieved during optimization
gmyc.edit <- function (tr, method = "single", interval = c(0, 5), quiet = FALSE, reltol = sqrt(.Machine$double.eps)) {
    local.env <- environment()
    read.data <- function(z = 1) {
        bt <- -branching.times(tr)
        bt[bt > -1e-06] <- -1e-06
        names(bt) <- NULL
        assign("bt", bt, envir = local.env)
        assign("sb", sort(bt), envir = local.env)
        assign("numnod", length(bt), envir = local.env)
        assign("numtip", length(tr$tip.label), envir = local.env)
        assign("numall", length(bt) + length(tr$tip.label), envir = local.env)
        internod <- sb[2:numnod] - sb[1:numnod - 1]
        internod[numnod] <- 0 - sb[numnod]
        assign("internod", internod, envir = local.env)
        assign("nesting", sapply((numtip + 1):numall, nesting.nodes), 
            envir = local.env)
        assign("nested", sapply((numtip + 1):numall, nest.nodes), 
            envir = local.env)
        numnested <- array(NA, numnod)
        numnesting <- array(NA, numnod)
        des <- matrix(NA, nrow = numnod, ncol = 2)
        sis <- array(NA, numnod)
        for (j in (1:numnod)) {
            numnested[j] <- length(nested[[j]])
            numnesting[j] <- length(nesting[[j]])
            des[j, ] <- as.integer(tr$edge[, 2])[as.integer(tr$edge[, 
                1]) == j + numtip]
            if (des[j, 1] > numtip) {
                sis[des[j, 1] - numtip] <- des[j, 2]
            }
            if (des[j, 2] > numtip) {
                sis[des[j, 2] - numtip] <- des[j, 1]
            }
        }
        assign("numnested", numnested, envir = local.env)
        assign("numnesting", numnesting, envir = local.env)
        des[des <= numtip] <- NA
        assign("des", des, envir = local.env)
        assign("sis", sis, envir = local.env)
        ancs <- cbind(tr$edge[pmatch((1:numnod + numtip), tr$edge[, 
            2]), 1], (1:numnod + numtip))
        bt.ancs <- cbind(bt[ancs[, 1] - numtip], bt[ancs[, 2] - 
            numtip])
        assign("bt.ancs", bt.ancs, envir = local.env)
    }
    nest.nodes <- function(x, p = 0) {
        numtip <- length(tr$tip.label)
        nods <- array(NA, 0)
        desc <- as.integer(tr$edge[, 2][tr$edge[, 1] == x])
        if (desc[1] > numtip) {
            nods <- c(nods, desc[1], nest.nodes(desc[1]))
        }
        if (desc[2] > numtip) {
            nods <- c(nods, desc[2], nest.nodes(desc[2]))
        }
        if (length(nods) > 0) {
            return(nods)
        }
        else {
            return(NULL)
        }
    }
    nesting.nodes <- function(x, p = 0) {
        numtip <- length(tr$tip.label)
        nod <- array(NA, 0)
        if (x >= numtip + 2) {
            anc <- as.integer(tr$edge[, 1][tr$edge[, 2] == x])
        }
        else {
            anc <- 1
        }
        if (anc >= numtip + 2) {
            nod <- c(nod, anc, nesting.nodes(anc))
        }
        if (anc == numtip + 1) {
            nod <- c(nod, anc)
        }
        if (length(nod) > 0) {
            return(nod)
        }
        else {
            return(NULL)
        }
    }
    l.prep <- function() {
        assign("n", length(mrca.nodes), envir = local.env)
        i.mat <- matrix(0, ncol = numnod, nrow = (n + 1))
        s.nod <- matrix(0, ncol = numnod, nrow = (n + 1))
        assign("nod", nod.type[order(bt)], envir = local.env)
        for (i in (1:n)) {
            s.nod[i, mrca.nodes[i]] <- 2
            if (!is.null(nested[[mrca.nodes[i]]])) {
                s.nod[i, nested[[mrca.nodes[i]]] - numtip] <- 1
            }
            s.nod[i, ] <- s.nod[i, order(bt)]
            i.mat[i, ][s.nod[i, ] == 2] <- 2
            i.mat[i, ][s.nod[i, ] == 1] <- 1
            i.mat[i, ] <- cumsum(i.mat[i, ])
        }
        s.nod[s.nod == 2] <- 1
        i.mat <- i.mat * (i.mat - 1)
        s.nod[n + 1, ] <- nod == 0
        i.mat[n + 1, nod == 0] <- 1
        i.mat[n + 1, nod == 2] <- -1
        i.mat[n + 1, ] <- cumsum(i.mat[n + 1, ]) + 1
        assign("s.nod", s.nod, envir = local.env)
        assign("i.mat", i.mat, envir = local.env)
    }
    l.null <- function(p = 1) {
        i.div <- 2:(numnod + 1)
        i.div <- i.div * (i.div - 1)
        lambda.div <- numnod/sum(internod * i.div^p)
        assign("lambda.div", lambda.div, envir = local.env)
        lik <- i.div^p * lambda.div * exp(-i.div^p * lambda.div * 
            internod)
        return(sum(log(lik)))
    }
    l.null.yul <- function(p = 1) { # added (JP)
		s.yul <- 2:(numnod+1) # added (JP)
		lambda.yul <- numnod/s.yul^p %*% internod # added (JP)
        assign("lambda.yul", lambda.yul, envir = local.env) # added (JP)
        lik <- s.yul^p * lambda.yul * exp(-s.yul^p * lambda.yul * internod) # added (JP)
        return(sum(log(lik))) # added (JP)
    } # added (JP)
    l.min2 <- function(q = c(1, 1)) {
        p <- c(rep(q[1], n), q[2])
        lambda <- sum(s.nod[1:n, ])/sum(i.mat[1:n, ]^p[-(n + 
            1)] %*% internod)
        lambda <- c(rep(lambda, n), sum(s.nod[n + 1, ])/i.mat[n + 
            1, ]^p[n + 1] %*% internod)
        b <- t(i.mat^p) %*% lambda
        lik <- b * exp(-b * internod)
        lambda <- lambda
        assign("b", b, envir = local.env)
        assign("lik", lik, envir = local.env)
        assign("lambda", lambda, envir = local.env)
        return(sum(log(lik)))
    }
    gmyc.likelihood <- function() {
        l.prep()
        # l.results <- rep(NA, 7) # commented (JP)
        # x <- optim(c(1, 1), l.min2, method = "Nelder-Mead", control = list(fnscale = -1)) # commented (JP)
        # l.results[1:6] <- c(x$value, lambda[n + 1], lambda[1],  # commented (JP)
        #     x$par[2], x$par[1], as.integer(sum(s.nod[n + 1, ]) +  # commented (JP)
        #         1)) # commented (JP)
        l.results <- rep(NA, 9) # added (JP)
        x <- optim(c(1, 1), l.min2, method = "Nelder-Mead", control = list(fnscale = -1, reltol = reltol)) # added (JP)
	    l.results[c(1:6, 8:9)] <- c(x$value, lambda[n + 1], lambda[1], x$par[2], x$par[1], as.integer(sum(s.nod[n + 1, ]) + 1), x$convergence, x$counts['function']) # added (JP)
        l.results[7] <- n
        assign("temp.params", x$par, envir = local.env)
        return(l.results)
    }
    gmyc.exhaustive <- function() {
        re.combn <- function(node = numtip + 1, tr = tr, se = numtip + 
            1) {
            parents <- tr$edge[, 1]
            children <- tr$edge[, 2]
            temp.list <- list()
            ch <- children[which(parents == node)]
            se <- se[-(which(se == node))]
            se <- c(se, ch)
            temp.list <- c(temp.list, list(se))
            if (any(ch > numtip)) {
                for (cc in ch) {
                  temp.temp.list <- temp.list
                  if (cc > numtip) {
                    for (tl in temp.temp.list) {
                      temp.list <- c(temp.list, re.combn(node = cc, 
                        tr = tr, se = tl))
                    }
                  }
                }
            }
            else {
            }
            return(temp.list)
        }
        remove.tips <- function(vec) {
            return(vec[vec > numtip])
        }
        MRCA <- lapply(re.combn(tr = tr), remove.tips)
        MRCA <- MRCA[-length(MRCA)]
        if (!quiet) {
            cat("exhaustive search", "\npossible combinations of MRCAs are", 
                length(MRCA), "\n")
        }
        l.results <- c()
        l.mrca <- lapply(MRCA, "-", numtip)
        assign("temp.params", c(1, 1), envir = local.env)
        count <- 1
        for (nodes in MRCA) {
            mrca.nodes <- nodes - numtip
            assign("mrca.nodes", mrca.nodes, envir = local.env)
            nod.type <- rep(0, numnod)
            for (j in 1:length(mrca.nodes)) {
                nod.type[mrca.nodes[j]] <- 2
                nod.type[nested[[mrca.nodes[j]]] - numtip] <- 1
            }
            assign("nod.type", nod.type, envir = local.env)
            l.results <- rbind(l.results, gmyc.likelihood())
            if (!quiet) {
                cat(l.results[count, 1], "\n")
                count <- count + 1
            }
        }
        return(list(l.results, l.mrca))
    }
    gmyc.single <- function() {
        # l.results <- matrix(ncol = 7, nrow = (nthresh)) # commented (JP)
        l.results <- matrix(ncol = 9, nrow = (nthresh + 1)) # added (JP)
        l.mrca <- list()
        results <- matrix(nrow = nthresh, ncol = numnod)
        x <- optimise(l.null, interval = interval, maximum = 1)
        # l.results[1, c(6:7, 1:2, 4)] <- c(as.integer(1), as.integer(1),  # commented (JP)
        #    x$objective, lambda.div, x$maximum) # commented (JP)
        l.results[1, c(6:7, 1, 3, 5, 8:9)] <- c(as.integer(1), as.integer(1), x$objective, lambda.div, x$maximum, 0, 0)  # added (JP)
        ml <- 0
        if (!quiet) {
            cat("node\t", "T\t", "loglik\n")
        }
        stthresh <- 2
        while (sb[stthresh] == sb[1]) {
            stthresh <- stthresh + 1
        }
        assign("temp.params", c(1, 1), envir = local.env)
        for (j in (stthresh:nthresh)) {
            threshy <- sb[j]
            tmp <- (bt.ancs[, 1] < threshy) & (bt.ancs[, 2] >= 
                threshy)
            nod.type <- tmp + (bt >= threshy)
            mrca.nodes <- which(nod.type == 2)
            if (nod.type[1] == 1) 
                nod.type[1] <- 2
            mrca.nodes <- which(nod.type == 2)
            assign("mrca", mrca, envir = local.env)
            assign("mrca.nodes", mrca.nodes, envir = local.env)
            assign("nod.type", nod.type, envir = local.env)
            l.mrca[[j]] <- mrca.nodes
            l.results[j, ] <- gmyc.likelihood()
            if (!quiet) {
            #    cat(j, threshy, l.results[j, 1], "\n") # commented (JP)
                cat(j, threshy, l.results[j, c(1, 8, 9)], length(l.mrca[[j]]), "\n") # added (JP)
            }
        }
        x <- optimise(l.null.yul, interval = interval, maximum = 1) # added (JP)
        # l.results[1, c(6:7, 1:2, 4)] <- c(as.integer(1), as.integer(1),  # commented (JP)
        #    x$objective, lambda.div, x$maximum) # commented (JP)
        l.results[nthresh + 1, c(6:7, 1:2, 4, 8:9)] <- c(as.integer(1), as.integer(1), x$objective, lambda.yul, x$maximum, 0, 0)  # added (JP)
        return(list(l.results, l.mrca))
    }
    gmyc.multi <- function() {
        assign("temp.params", c(1, 1), envir = local.env)
        renew.mrca <- function(n) {
            f.renew <- function(n) {
                parents <- tr$edge[, 1]
                children <- tr$edge[, 2]
                renewed <- list()
                if (length(n) > 1) {
                  pair <- combn(n, 2)
                  for (i in 1:length(n)) {
                    parent.node <- parents[children == n[i]]
                    sibling.nodes <- children[parents == parent.node]
                    sibling <- sibling.nodes[sibling.nodes != 
                      n[i]]
                    if (sibling <= numtip) {
                      renewed <- c(renewed, list(c(n[-i], parent.node)))
                    }
                    else if (any(n == sibling)) {
                      renewed <- c(renewed, list(c(n[-c(i, which(n == 
                        sibling))], parent.node)))
                    }
                  }
                }
                return(unique(renewed))
            }
            d.renew <- function(n) {
                parents <- tr$edge[, 1]
                children <- tr$edge[, 2]
                renewed <- list()
                if (length(n) > 1) {
                  for (i in 1:length(n)) {
                    child.nodes <- children[parents == n[i]]
                    if (any(child.nodes > numtip)) {
                      renewed <- c(renewed, list(c(n[-i], child.nodes[child.nodes > 
                        numtip])))
                    }
                  }
                }
                return(renewed)
            }
            return(c(f.renew(n), d.renew(n)))
        }
        select.start.re <- function(time, ml) {
            result <- c()
            result.mrca <- list()
            m <- mean(time)
            left <- time[time < m]
            right <- time[time >= m]
            mid <- time[time >= mean(left) & time < mean(right)]
            part <- list(left, mid, right)
            start <- c(mean(left), m, mean(right))
            if (!any(is.na(start)) || length(unique(start)) == 
                1) {
                improve <- c(FALSE, FALSE, FALSE)
                num.improve <- c(0, 0, 0)
                for (i in 1:length(start)) {
                  if (!quiet) {
                    cat("start at", start[i], "\n")
                  }
                  temp <- renew.likelihood.from.thresh(start[i])
                  lik <- temp[[1]]
                  mrca <- temp[[2]]
                  if (any(lik[, 1] > ml)) {
                    if (!quiet) {
                      cat("improvement found\n")
                    }
                    result <- rbind(result, lik[which(lik[, 1] > 
                      ml), ])
                    result.mrca <- c(result.mrca, mrca[which(lik[, 
                      1] > ml)])
                    improve[i] <- TRUE
                    num.improve[i] <- length(lik[lik[, 1] > ml, 
                      1])
                  }
                }
                if (!quiet) {
                  cat(num.improve, "\n")
                }
                if (!all(num.improve == 0)) {
                  temp <- select.start.re(part[[which.max(num.improve)]], 
                    max(result[, 1]))
                  result <- rbind(result, temp[[1]])
                  result.mrca <- c(result.mrca, temp[[2]])
                }
            }
            return(list(result, result.mrca))
        }
        select.start <- function(time, ml) {
            result <- c()
            result.mrca <- list()
            m <- mean(time)
            left <- time[time < m]
            right <- time[time >= m]
            mid <- time[time >= mean(left) & time < mean(right)]
            part <- list(left, mid, right)
            start <- c(mean(left), m, mean(right))
            for (i in 1:length(start)) {
                if (!quiet) {
                  cat("start at", start[i], "\n")
                }
                temp <- renew.likelihood.from.thresh(start[i])
                lik <- temp[[1]]
                mrca <- temp[[2]]
                if (any(lik[, 1] > ml)) {
                  if (!quiet) {
                    print("improvement found\n")
                  }
                  result <- rbind(result, lik[which(lik[, 1] > 
                    ml), ])
                  result.mrca <- c(result.mrca, mrca[which(lik[, 
                    1] > ml)])
                  if (length(part[[i]]) > 1) {
                    if (!quiet) {
                      cat("recursively trying part of", start[i], 
                        "\n")
                    }
                    temp <- select.start(part[[i]], max(lik[, 
                      1]))
                    result <- rbind(result, temp[[1]])
                    result.mrca <- c(result.mrca, temp[[2]])
                  }
                }
            }
            return(list(result, result.mrca))
        }
        renew.likelihood.from.thresh <- function(start) {
            l.results <- c()
            l.mrca <- list()
            max.lik <- NA
            mrca <- array(FALSE, numnod)
            mrca[2:numnod] <- (bt[as.integer(tr$edge[, 1][tr$edge[, 
                2] > numtip]) - numtip] < start) & (bt[as.integer(tr$edge[, 
                2][tr$edge[, 2] > numtip]) - numtip] >= start)
            nod.type <- (bt >= start) + mrca
            if (nod.type[1] == 1) 
                nod.type[1] <- 2
            mrca.nodes <- which(nod.type == 2)
            assign("mrca", mrca, envir = local.env)
            assign("mrca.nodes", mrca.nodes, envir = local.env)
            assign("nod.type", nod.type, envir = local.env)
            initial.mrca <- mrca.nodes
            while (TRUE) {
                found <- FALSE
                max.mrca <- initial.mrca
                for (nodes in renew.mrca(initial.mrca + numtip)) {
                  mrca.nodes <- nodes - numtip
                  if (mrca.nodes[[1]] != 1) {
                    assign("mrca.nodes", mrca.nodes, envir = local.env)
                    nod.type <- rep(0, numnod)
                    for (j in 1:length(mrca.nodes)) {
                      nod.type[mrca.nodes[j]] <- 2
                      nod.type[nested[[mrca.nodes[j]]] - numtip] <- 1
                    }
                    assign("nod.type", nod.type, envir = local.env)
                    res <- gmyc.likelihood()
                    if (!is.na(max.lik)) {
                      if (max.lik <= res[1]) {
                        max.lik <- res[1]
                        l.results <- rbind(l.results, res)
                        max.mrca <- nodes - numtip
                        l.mrca <- c(l.mrca, list(max.mrca))
                        if (!quiet) {
                          cat(max.lik, "\n")
                        }
                        found <- TRUE
                      }
                    }
                    else {
                      max.lik <- res[1]
                      l.results <- rbind(l.results, res)
                      max.mrca <- nodes - numtip
                      l.mrca <- c(l.mrca, list(max.mrca))
                      found <- TRUE
                    }
                  }
                }
                if (found) {
                  initial.mrca <- max.mrca
                }
                else {
                  if (!quiet) {
                    cat("break\n")
                  }
                  break
                }
            }
            rownames(l.results) <- NULL
            return(list(l.results, l.mrca))
        }
        l.results <- c()
        l.mrca <- list()
        x <- optimise(l.null, interval = interval, maximum = 1)
        # l.results <- rbind(l.results, c(x$objective, lambda.div,  # commented (JP)
        #     NA, x$maximum, NA, as.integer(1), as.integer(1))) # commented (JP)
        l.results <- rbind(l.results, c(x$objective, NA, lambda.div, NA, x$maximum, as.integer(1), as.integer(1), 0, 0)) # added (JP)
        l.mrca <- c(l.mrca, list(numtip + 1))
        if (!quiet) {
            cat("null likelihood\n")
            cat(l.results[1, 1], "\n")
        }
        if (!quiet) {
            cat("\nGMYC likelihood\n")
        }
        temp <- select.start.re(sb, l.results[1, 1])
        l.results <- rbind(l.results, temp[[1]])
        l.mrca <- c(l.mrca, temp[[2]])
        x <- optimise(l.null.yul, interval = interval, maximum = 1) # added (JP)
        l.results <- rbind(l.results, c(x$objective, lambda.yul, NA, x$maximum, NA, as.integer(1), as.integer(1), 0, 0))  # added (JP)
        return(list(l.results, l.mrca))
    }
    if (!is.ultrametric(tr)) {
        stop("Your ultrametric tree is not ultrametric, please check")
    }
    if (!is.binary.tree(tr)) {
        stop("Your input tree is not fully bifurcating, please resolve with zero branch lengths")
    }
    read.data()
    ntrees <- 1
    numnod <- max(as.integer(tr$edge[, 1])) - numtip
    nthresh <- numnod
    if (method == "single" || method == "s") {
        method <- "single"
        temp <- gmyc.single()
        l.results <- temp[[1]]
        l.mrca <- temp[[2]]
    }
    else if (method == "multiple" || method == "m") {
        method <- "multiple"
        temp <- gmyc.multi()
        l.results <- temp[[1]]
        l.mrca <- temp[[2]]
    }
    else if (method == "exhaustive" || method == "e") {
        method <- "exhaustive"
        temp <- gmyc.exhaustive()
        l.results <- temp[[1]]
        l.mrca <- temp[[2]]
    }
    else {
        stop("Invalid name of optimiztion method. Only single(s) or multiple(m) and exhaustive(e) are available.")
    }
    cat("\n", date(), "\n", sep = "")
    cat("finish.\n")
    result <- list()
    result[["method"]] <- method
    result[["likelihood"]] <- l.results[, 1]
    result[["parameters"]] <- l.results[, 2:5]
    colnames(result[["parameters"]]) <- c("lambda.div", "lambda.coal", 
        "p.div", "p.coal")
    result[["entity"]] <- l.results[, 6]
    result[["cluster"]] <- l.results[, 7]
	result[["convergence"]] <- l.results[, 8] # added (JP)
    if (method == "single") {
        result[["MRCA"]] <- l.mrca
        result[["threshold.time"]] <- sb
    }
    else if (method == "multiple" || method == "exhaustive") {
        reduce.threshold <- function(mrcas) {
            parent <- tr$edge[, 1]
            child <- tr$edge[, 2]
            thresh.group <- list()
            thresh.time <- c()
            mrcas <- mrcas + numtip
            k <- 1
            while (TRUE) {
                times <- bt[mrcas - numtip]
                thresh1.time <- min(times)
                thresh1.node <- mrcas[which.min(times)]
                mrcas <- mrcas[-which.min(times)]
                if (length(mrcas) == 0) {
                  thresh.time <- c(thresh.time, thresh1.time)
                  thresh.group[[k]] <- thresh1.node
                  break
                }
                member <- thresh1.node
                del <- c()
                for (i in 1:length(mrcas)) {
                  par.nod <- parent[child == mrcas[i]]
                  t.par <- bt[par.nod - numtip]
                  if (t.par < thresh1.time) {
                    member <- c(member, mrcas[i])
                    del <- c(del, i)
                  }
                }
                thresh.time <- c(thresh.time, thresh1.time)
                thresh.group[[k]] <- member
                k <- k + 1
                if (length(del) != 0) {
                  mrcas <- mrcas[-del]
                }
                if (length(mrcas) == 0) {
                  break
                }
            }
            return(thresh.time)
        }
        result[["MRCA"]] <- l.mrca
        result[["MRCA"]] <- l.mrca
        result[["threshold.time"]] <- lapply(l.mrca, reduce.threshold)
    }
    result[["tree"]] <- tr
    class(result) <- "gmyc"
    return(result)
}


# estimates AIC scores, Akaike weights for gmyc.edit() objects ('res1','res2')
# can use on single- or multiple-threshold model results individually but gives a warning
# when including both model types, specify single-threshold model first -> currently gives an error if multiple-threshold model is specified as the first argument (need to fix; has to do with null.yule referring to res.sing during structuring of 'out' table, but res.mult during structuring of 'parameter' table when res1=mult and res2=sing)
model.scores<-function(res1,res2=NULL,deltaAICc=2){
	ifelse(res1[['method']]=='single',res.sing<-res1,ifelse(res1[['method']]=='multiple',res.mult<-res1,stop("\"res1\" is not a valid \"gmyc\" object.")))
	if(!is.null(res2)){
		ifelse(res2[['method']]=='single',res.sing<-res2,ifelse(res2[['method']]=='multiple',res.mult<-res2,stop("\"res2\" is not a valid \"gmyc\" object.")))
		} else {
		if(res1[['method']]=='single') res.mult<-NULL
		if(res1[['method']]=='multiple') res.sing<-NULL
		}
	if(!is.null(res.sing)){
		res<-res.sing
		lik<-res$likelihood
		orig.step<-1:length(lik)
		LR<-pvalue<-pval.mult.sing<-k<-AICc<-numeric(length(lik))
		LR[1]<-pvalue[1]<-pval.mult.sing[1]<-NA
		n<-length(res$tree$tip.label)
		k[1]<-2
		AICc[1]<-2*k[1]-2*lik[1]+((2*k[1]*(k[1]+1))/(n-k[1]-1))
		for(i in 2:(length(lik)-1)){
			LR[i]<-2*(lik[i]-lik[1])
			pvalue[i]<-1-pchisq(LR[i],3)
			pval.mult.sing[i]<-NA
			k[i]<-3+2
			AICc[i]<-2*k[i]-2*lik[i]+((2*k[i]*(k[i]+1))/(n-k[i]-1))
			}
		LR[length(lik)]<-pvalue[length(lik)]<-pval.mult.sing[length(lik)]<-NA
		k[length(lik)]<-2
		AICc[length(lik)]<-2*k[length(lik)]-2*lik[length(lik)]+((2*k[length(lik)]*(k[length(lik)]+1))/(n-k[length(lik)]-1))
		out.sing<-data.frame(step=orig.step,likelihood=lik,LR=LR,k=k,pvalue=pvalue,pval.mult.sing=pval.mult.sing,AICc=AICc)
		out.sing$method<-'single'
		out.sing$method[1]<-'null.coal'
		out.sing$method[length(lik)]<-'null.yule'
		}
	if(!is.null(res.mult)){
		res<-res.mult
		lik<-res$likelihood
		orig.step<-1:length(lik)
		if(!is.null(res.sing)) sing.lik<-max(res.sing$likelihood)
		LR<-pvalue<-pval.mult.sing<-k<-AICc<-numeric(length(lik))
		LR[1]<-pvalue[1]<-pval.mult.sing[1]<-NA
		n<-length(res$tree$tip.label)
		k[1]<-2
		AICc[1]<-2*k[1]-2*lik[1]+((2*k[1]*(k[1]+1))/(n-k[1]-1))
		for(i in 2:(length(lik)-1)){
			LR[i]<-2*(lik[i]-lik[1])
			pvalue[i]<-1-pchisq(LR[i],3+length(res$threshold.time[[i]])-1)
			ifelse(!is.null(res.sing),pval.mult.sing[i]<-1-pchisq(2*(lik[i]-sing.lik),length(res$threshold.time[[i]])),pval.mult.sing[i]<-NA)  #if 1 parameter per threshold
			#ifelse(!is.null(res.sing),pval.mult.sing[i]<-1-pchisq(2*(lik[i]-sing.lik),3*(length(res$threshold.time[[i]])-1)),pval.mult.sing[i]<-NA)  #if 3 parameters per threshold
			k[i]<-3+length(res$threshold.time[[i]])-1+2  #if 1 parameter per threshold
			#k[i]<-3*length(res$threshold.time[[i]])+2  #if 3 parameters per threshold
			AICc[i]<-2*k[i]-2*lik[i]+((2*k[i]*(k[i]+1))/(n-k[i]-1))
			}
		LR[length(lik)]<-pvalue[length(lik)]<-pval.mult.sing[length(lik)]<-NA
		k[length(lik)]<-2
		AICc[length(lik)]<-2*k[length(lik)]-2*lik[length(lik)]+((2*k[length(lik)]*(k[length(lik)]+1))/(n-k[length(lik)]-1))
		out.mult<-data.frame(step=orig.step,likelihood=lik,LR=LR,k=k,pvalue=pvalue,pval.mult.sing=pval.mult.sing,AICc=AICc)
		out.mult$method<-'multiple'
		out.mult$method[1]<-'null.coal'
		out.mult$method[length(lik)]<-'null.yule'
		out.mult<-out.mult[which(out.mult$k!=5),] # eliminates cases where multiple-threshold model finds single-likelihood solution by chance
		}
	if(!is.null(res.sing)&!is.null(res.mult)){
		out<-rbind(out.sing[c(1,nrow(out.sing),2:(nrow(out.sing)-1)),],out.mult[2:(nrow(out.mult)-1),])
		sing.mrcas<-res.sing[['MRCA']][out$step[which(out$method!='multiple')]]
		mult.mrcas<-res.mult[['MRCA']][out$step[which(out$method=='multiple')]]
		#mrcas<-list(sing.mrcas=sing.mrcas,mult.mrcas=mult.mrcas)
		mrcas<-c(sing.mrcas,mult.mrcas)
		}
	if(!is.null(res.sing)&is.null(res.mult)){
		out<-out.sing[c(1,nrow(out.sing),2:(nrow(out.sing)-1)),]
		mrcas<-res.sing[['MRCA']][out$step]
		}
	if(is.null(res.sing)&!is.null(res.mult)){
		out<-out.mult[c(1,nrow(out.mult),2:(nrow(out.mult)-1)),]
		mrcas<-res.mult[['MRCA']][out$step]
		}
	mrcas[[1]]<-1
	rownames(out)<-1:nrow(out)
	
	nclust<-nent<-numeric(0)
	param<-res1$parameters[0,]
	for(i in 1:nrow(out)){
		if(out[i,'method']=='null.coal'){
			nclust<-c(nclust,1)
			nent<-c(nent,1)
			param<-rbind(param,res1$parameters[out[i,'step'],])
			}
		if(out[i,'method']=='null.yule'){
			nclust<-c(nclust,0)
			nent<-c(nent,length(res$tree$tip.label))
			param<-rbind(param,res1$parameters[out[i,'step'],])
			}
		if(out[i,'method']=='single'){
			nclust<-c(nclust,res.sing$cluster[out[i,'step']])
			nent<-c(nent,res.sing$entity[out[i,'step']])
			param<-rbind(param,res.sing$parameters[out[i,'step'],])
			}
		if(out[i,'method']=='multiple'){
			nclust<-c(nclust,res.mult$cluster[out[i,'step']])
			nent<-c(nent,res.mult$entity[out[i,'step']])
			param<-rbind(param,res.mult$parameters[out[i,'step'],])
			}
		}
	if(is.null(dimnames(param))) dimnames(param)[[2]]<-c('lambda.div','lambda.coal','p.div','p.coal') # for gmyc runs from original code
	
	del<-out$AICc-out$AICc[which.min(out$AICc)]
	w<-exp(-0.5*del)/sum(exp(-0.5*del))
	aic.wts<-list(delta.AICc=del,Akaike.weights=w)
	
	species.est<-function(parameter,del=del,w=w){
		p<-parameter
		p[is.na(p)]<-0
		av.the<-sum(w*p)
		var.the<-sum(w*sqrt((p-av.the)^2))
		data.frame(estimate=av.the,variance=var.the)
		}
	species.tab<-data.frame(estimate=numeric(0),variance=numeric(0))
	for(i in list(nclust,nent)) species.tab<-rbind(species.tab,species.est(i,del,w))
	row.names(species.tab)<-c('clusters','entities')
	
	param.est<-function(parameter,del=del,w=w){
		p<-param[,parameter]
		p[is.na(p)]<-0
		av.the<-sum(w*p)
		var.the<-sum(w*sqrt((p-av.the)^2))
		data.frame(estimate=av.the,variance=var.the)
		}
	param.tab<-data.frame(estimate=numeric(0),variance=numeric(0))
	for(i in c('lambda.div','lambda.coal','p.div','p.coal')) param.tab<-rbind(param.tab,param.est(i,del,w))
	row.names(param.tab)<-c('lambda.div','lambda.coal','p.div','p.coal')
	
	df.tab<-data.frame(estimate=res1$tree$Nnode+1-sum(w*out$k),variance=sum(w*sqrt((out$k-sum(w*out$k))^2)))
	row.names(df.tab)<-'df'
	modelAv.tab<-rbind(species.tab,param.tab,df.tab)

	warns<-'No warnings'
	if(is.null(res.sing)) warns<-'Model scores only calculated for multiple-threshold model'
	if(is.null(res.mult)) warns<-'Model scores only calculated for single-threshold model'

	cat('\nSummary of models (model averaged):\n')
	print(modelAv.tab)
	cat('\n',warns,'\n')
	
	out.sig<-out[which(out$fit!=''),]
	all.out<-list(scores=out,aic.weights=aic.wts,clusters=nclust,entities=nent,parameters=param,model.averages=modelAv.tab,tree=res$tree,MRCA=mrcas,deltaAICc=deltaAICc,noSigModels=nrow(out.sig),noModels=nrow(out),warnings=warns)
	all.out
	}



# for visualizing a summary of the model.scores() result
# deltaAICc argument allows for visualizing models within a certain distance of the best-fit model
summary.scores<-function(model.scores.result, deltaAICc=2){
	x<-model.scores.result
	out<-x$scores
	out$fit<-''
	out$fit[which(out$AICc<=min(out$AICc)+deltaAICc)]<-'*'
	out$fit[which.min(out$AICc)]<-'**'
	modelAv.tab<-x$model.averages
	
	warns<-x$warnings

	tab<-out[which(out$fit!=''),c('method','step','likelihood','k','AICc')]
	tab$clusters<-x$clusters[which(out$fit!='')]
	tab$entities<-x$entities[which(out$fit!='')]
	tab$delta.AICc<-x$aic.weights$delta.AICc[which(out$fit!='')]
	tab$weights<-x$aic.weights$Akaike.weights[which(out$fit!='')]

	cat('\nBest-fit and other models (within delta AICc =',deltaAICc,') model statistics:\n\n')
	print(tab[order(tab$delta.AICc),])
	cat('\n\nSummary of models (model averaged):\n')
	print(modelAv.tab)
	cat('\n',warns,'\n')
	}


# estimate probabilities that tips co-occur within a cluster, using Akaike weights
# can use 'deltaAICc' argument to exclude models with low Akaike weights (e.g. 0.99 includes those models with cumulative weight up to 0.99) -> significantly speeds up function for large trees
# first element of resulting list is the tree
# second element ('clusterMatrix.avg') contains the probability matrix
# also estimates lower and upper range of confidence interval, stores each in a separate matrix
cluster.probs<-function(model.scores.result,deltaAICc=NULL){
	res<-model.scores.result
	if(is.null(deltaAICc)) keep<-1:nrow(res$scores)
	if(deltaAICc>=1) keep<-which(res[['aic.weights']][['delta.AICc']]<deltaAICc)
	if(deltaAICc<1) {
		wts.ord<-order(res[['aic.weights']][['Akaike.weights']],decreasing=T)
		cum.sum<-res[['aic.weights']][['Akaike.weights']][wts.ord[1]]
		for(i in 2:length(wts.ord)) cum.sum<-c(cum.sum,cum.sum[i-1]+res[['aic.weights']][['Akaike.weights']][wts.ord[i]])
		keep<-wts.ord[which(cum.sum<=deltaAICc)]
		}
	cat('Calculating probabilities from',length(keep),'of',nrow(res$scores),'models.\n')
	ask<-ask(msg='Do you want to continue? (y/n): ')
	if(ask=='n') cat('Result not generated.')
	if(ask=='y'){
		cat('Calculating clade probabilities (this may take awhile) ...\n')
		out<-res$scores[keep,]
		del<-out$AICc-out$AICc[which.min(out$AICc)]
		w<-exp(-0.5*del)/sum(exp(-0.5*del))
		mat.clust1<-matrix(0,nrow=length(res$tree$tip.label),ncol=length(res$tree$tip.label),dimnames=list(res$tree$tip.label,res$tree$tip.label))
		out.clust<-vector('list',length(keep))
		for(i in 1:length(keep)){
			if(out[i,'method']=='null.yule'){
				cluster.names<-vector('list',length=length(res$tree$tip.label))
				for(j in 1:length(cluster.names)) cluster.names[[j]]<-res$tree$tip.label[j]
				out.clust[[i]]<-cluster.names
				} else {
				cluster.names<-vector('list',length=length(res$MRCA[[keep[i]]]))
				for(j in 1:length(cluster.names)) cluster.names[[j]]<-node.leaves(res$tree,res$MRCA[[keep[i]]][j]+length(res$tree$tip.label))
				out.clust[[i]]<-cluster.names
				}
			}
		mat.clust<-mat.clust1
		for(i in 1:length(out.clust)){
			clust<-out.clust[[i]]
			z1<-mat.clust1
			for(j in 1:length(clust)){
				group<-clust[[j]]
				z1[group[1:length(group)],group[1:length(group)]]<-1
				}
			mat.clust<-mat.clust+w[i]*z1
			}
		mat.clust[which(mat.clust>1)]<-1 # issue with approximation? --> mat.clust[which(mat.clust>1)]   [1] 1 1 1 1 1 ...
		mat.clust.avg<-mat.clust
		ci95<-1.96*sqrt((mat.clust.avg*(1-mat.clust.avg))/nrow(out))
		mat.clust.upper<-mat.clust.avg+ci95
		mat.clust.lower<-mat.clust.avg-ci95
		cat('\nFinished calculating clade probabilities.\n')
		all.out<-list(tree=res$tree,clusterMatrix.avg=mat.clust.avg,clusterMatrix.lowerCI=mat.clust.lower,clusterMatrix.upperCI=mat.clust.upper)
		all.out
		}
	}



# for plotting phylogeny and cluster probabilities from cluster.probs() result to the left of each node
# 'clade' is a vector of tip labels that can be used to limit the estimate to a subsample of isolates in the tree
# 'CI' is a logical element indicating whether the confidence interval should included below the node
# 'main' and 'sub' are character elements specifying the heading and subheading
# works best if plotted to file with large page dimensions, for example with pdf()
plot.probs<-function(cluster.probs.result,clade=NULL,CI=F,main=NULL,sub=NULL){
	x<-cluster.probs.result
	plot(x$tree)
	title(main=main,sub=sub)
	prob.ave<-x$clusterMatrix.avg
	node.ave<-numeric(0)
	if(CI){ 
		prob.low<-x$clusterMatrix.lowerCI
		prob.upp<-x$clusterMatrix.upperCI
		node.low<-node.upp<-numeric(0)
		}
	node.lab<-character(0)
	for(i in 1:x$tree$Nnode){
		node<-i+length(x$tree$tip.label)
		tips<-node.leaves(x$tree,node)
		ave1<-round(min(prob.ave[tips,tips]),3)
		if(CI){
			low1<-round(min(prob.low[tips,tips]),3)
			upp1<-round(min(prob.upp[tips,tips]),3)
			}
		ifelse(CI,node.lab<-c(node.lab,paste(ave1,'\n(',low1,',',upp1,')',sep='')),node.lab<-c(node.lab,paste(ave1,'\n',sep='')))
		}
	nodelabels(node.lab,adj=1,frame='none',bg='none',cex=1)
	}



# estimate richness, Shannon diversity from model.scores result, with or without a sample matrix ('samples.result')
# can use 'deltaAICc' argument to exclude models with low Akaike weights (e.g. 0.99 includes those models with cumulative weight up to 0.99) -> significantly speeds up function for large trees
# 'haplotypes' is used in special cases to link duplicated sequences that are represented in the sample matrix but are not represented in the tree 
# 'clade' is a vector of tip labels that can be used to limit the estimate to a subsample of isolates in the tree
sample.est<-function(model.scores.result,samples.result=NULL,haplotypes=F,deltaAICc=NULL,clade=NULL){
	#functions for extracting cluster information from the gmyc output
	clusters<-function(){ # for extracting a list of tips belonging to each gmyc cluster
		cluster.names<-vector('list',length=length(res$MRCA[[i]]))
		for(j in 1:length(cluster.names)) cluster.names[[j]]<-node.leaves(res$tree,res$MRCA[[i]][j]+length(res$tree$tip.label))
		cluster.names
		}
	singles<-function(){ # for extracting a vector of singletons from the gmyc model (each of which may refer to multiple haplotypes)
		cluster.tips<-unlist(clus1)
		all.tips<-res$tree$tip.label
		single.tips<-res$tree$tip.label[-c(match(cluster.tips,all.tips))]
		single.tips
		}

	x<-model.scores.result
	if(is.null(deltaAICc)) keep<-1:nrow(x$scores)
	if(!is.null(deltaAICc)){
		if(deltaAICc>=1) keep<-which(x[['aic.weights']][['delta.AICc']]<deltaAICc)
		if(deltaAICc<1) {
			wts.ord<-order(x[['aic.weights']][['Akaike.weights']],decreasing=T)
			cum.sum<-x[['aic.weights']][['Akaike.weights']][wts.ord[1]]
			for(i in 2:length(wts.ord)) cum.sum<-c(cum.sum,cum.sum[i-1]+x[['aic.weights']][['Akaike.weights']][wts.ord[i]])
			keep<-wts.ord[which(cum.sum<=deltaAICc)]
			}
		}
	cat('Calculating estimates from',length(keep),'of',nrow(x$scores),'models.\n')
	ask<-ask(msg='Do you want to continue? (y/n): ')
	if(ask=='n') cat('Result not generated.')
	if(ask=='y'){
		cat('Calculating estimates (this may take awhile) ...\n')
		if(!is.null(clade)) { # prunes clade from tree and creates new MRCA list corresponding to nodes within the pruned tree
			z<-mrca(x$tree)
			z<-z[match(clade,x$tree$tip.label),match(clade,x$tree$tip.label)]
			diag(z)<-max(z)+1
			clade.tree<-extract.clade(x$tree,min(z))
			z<-z-length(x$tree$tip.label)
			z<-sort(unique(as.numeric(z)))
			z<-z[1:(length(z)-1)]
			clade.MRCA<-vector('list',keep)
			for(i in 1:length(keep)){
				zz<-match(z,x$MRCA[[keep[i]]],nomatch=0)
				clade.MRCA[[i]]<-x$MRCA[[keep[i]]][zz]-(min(z)-1)
				if(length(clade.MRCA[[i]])==0) clade.MRCA[[i]]<-1
				}
			res<-list(tree=clade.tree,MRCA=clade.MRCA)
			}
		if(is.null(clade)) res<-list(tree=x$tree,MRCA=x$MRCA[keep])
		if(is.null(samples.result)) samples.result<-matrix(1,ncol=1,nrow=length(res$tree$tip.label),dimnames=list(res$tree$tip.label,'phylo'))
		clus.mat<-ents.mat<-shan.mat<-matrix(0,nrow=length(keep),ncol=ncol(samples.result),dimnames=list(keep,colnames(samples.result)))
		# abun.list<-vector('list',length(keep))  # if want to keep abundance
		for(i in 1:length(keep)){
			if(x$scores$method[keep[i]]=='null.yule'){
				ents1<-as.list(res$tree$tip.label)
				# abun.list[[i]]<-vector('list',ncol(samples.result))  # if want to keep abundance
				for(j in 1:ncol(samples.result)){
					abun.vec<-numeric(0)
					for(k in 1:length(ents1)){
						if(haplotypes) no.k<-sum(samples.result[match(ents1[[k]],rownames(samples.result)),j])
						if(!haplotypes){
							spp.vec<-samples.result[match(ents1[[k]],rownames(samples.result)),j]
							no.k<-length(spp.vec[which(spp.vec>0)])
							}
						if(no.k>=1) ents.mat[i,j]<-ents.mat[i,j]+1
						abun.vec<-c(abun.vec,no.k)
						}
					# abun.list[[i]][[j]]<-abun.vec  # if want to keep abundance
					shan.mat[i,j]<-diversity(matrix(abun.vec,nrow=1),index='shannon')
					}
				}
			else {
				clus1<-clusters()
				ents1<-c(clus1,singles())
				for(j in 1:ncol(samples.result)) for(k in 1:length(clus1)) if(sum(samples.result[match(clus1[[k]],rownames(samples.result)),j])>=1) clus.mat[i,j]<-clus.mat[i,j]+1
				for(j in 1:ncol(samples.result)){
					abun.vec<-numeric(0)
					for(k in 1:length(ents1)){
						if(haplotypes) no.k<-sum(samples.result[match(ents1[[k]],rownames(samples.result)),j])
						if(!haplotypes){
							spp.vec<-samples.result[match(ents1[[k]],rownames(samples.result)),j]
							no.k<-length(spp.vec[which(spp.vec>0)])
							}
						if(no.k>=1) ents.mat[i,j]<-ents.mat[i,j]+1
						abun.vec<-c(abun.vec,no.k)
						}
					# abun.list[[i]][[j]]<-abun.vec  # if want to keep abundance
					shan.mat[i,j]<-diversity(matrix(abun.vec,nrow=1),index='shannon')
					}
				}
			}

		del<-x[['aic.weights']][['delta.AICc']][keep]
		w<-exp(-0.5*del)/sum(exp(-0.5*del))

		species.est<-function(parameter,del=del,w=w){
			p<-parameter
			p[is.na(p)]<-0
			av.the<-sum(w*p)
			var.the<-sum(w*sqrt((p-av.the)^2)) 
			c(estimate=av.the,variance=var.the)
			}
		species.mat<-matrix(0,nrow=ncol(samples.result),ncol=6)
		for(i in 1:ncol(samples.result)) species.mat[i,]<-c(species.est(clus.mat[,i],del,w),species.est(ents.mat[,i],del,w),species.est(shan.mat[,i],del,w))
		dimnames(species.mat)<-list(colnames(samples.result),c('clusters.est','clusters.var','entities.est','entities.var','shannon.est','shannon.var'))
		cat('\nFinished calculating estimates.\n')
		species.mat
		}
	}



##END functions


