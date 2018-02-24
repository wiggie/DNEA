# ---
#   title: "CRIC Differential Network"
# author: "Jing Ma"
# date: "Oct 30,2017"

rm(list=ls())

## Loading packages
require(gdata)
require(zoo)
require(igraph)
require(glasso)
library(glmnet)
library(corpcor)
library(GSA)

library(parallel)
library(foreach)
library(doParallel)

source("lib/preprocess_lib.R")
source("lib/JSEM.R")
source("lib/netgsa_complex.R")

set.seed(1)

####################################  I/O #######################################
results_dir =
input_data =

results_dir = gsub("/$", "", results_dir)
if (!file.exists(results_dir))
   {
   cat("\n")
   stop(paste0("Results directory :", results_dir, " does not exist\n\n"))
   }

out_edgelist_full = paste(results_dir, "edgelist_full.csv", sep="/")
out_edgelist = paste(results_dir, "edgelist.csv", sep="/")
out_nodelist = paste(results_dir, "nodelist.csv", sep="/")
out_netgsa = paste(results_dir, "results.csv", sep="/")

out_clusters_pdf = paste(results_dir, "consensus_clusters.pdf", sep="/")

out_filtered_rdata = paste(results_dir, "filtered.RData", sep = "/")
out_tuning_data = paste(results_dir, "BIC_tuning.RData", sep = "/")
out_stable_networks = paste(results_dir, "stable_networks.RData", sep = "/")
out_netgsa_rda = paste(results_dir, "results_netgsa.RData", sep="/")



################################# Analysis -- workflow ###################################
# In the first step, select the tuning parameter **BIC=T, SS=F, runNetGSA=F **.
# In the second step, run stability selection by setting **BIC=F, SS=T, runNetGSA=F **.
# In the last step, run NetGSA and get graph by setting **BIC=F, SS=F, runNetGSA=T **.

BIC <- T          
SS <- F           
runNetGSA <- F    
savePlots <- T           # create PDF plots
pearsonFiltering <- F    # prefilter by Pearson correlations

################################# Analysis -- workflow ###################################
fdr.cutoff <- 0.2       # FDR threshold for defining whether a compound is differential
cor.threshold <- 0.7    # Threshold for pearson correlations
nCores <- 10            # Number of iterations per core
nreps <- 10             # Number of cores needed in the stability selection procedure
tau0 <- 0.5             # Threshold for getting consensus matrix
eps <- 1e-06            # Threshold for deciding zeros in the partial correlation matrix


################################### File checks  ########################################

if (!file.exists(results_dir))
   stop(paste0("Results directory does not exist ", results_dir))

if (!BIC && !SS && !runNetGSA)
  {
  cat("\n\n")
  stop("=>\tMust choose at least one analysis option\n\n")
  }

input_ok = TRUE

if (! isTRUE(runNetGSA))
    {
    cat("\n\n\nData preparation step...\n\n")
    
    this_input_ok = file.exists(input_data)
    if (! this_input_ok)
        stop(paste0("\t=>\tInput data " , input_data, " is missing\n"))
    else
        cat("\t=>\tInput data ", input_data, "is available\n")
    
    input_ok = (input_ok && this_input_ok)
    if (input_ok)
        {
        cat("\t=>\tFiltered data will be written to", out_filtered_rdata, "\n\n")
        }
    }

if (BIC)
    {
    if (SS || runNetGSA)
      stop("=>\tCannot choose two analysis options ")

    cat("\n\nTuning parameter calculation...  ")
    if (input_ok)
        cat("\n\t=>\tTuning data will be written to ", out_tuning_data, "\n\n")
    }

if (SS)
    {
    if (BIC || runNetGSA)
  	stop("=>\tCannot choose two analysis options ")
    
    cat("\n\nStability selection...\n\n")
    
    this_input_ok = file.exists(out_tuning_data)
    if (! this_input_ok)
        stop("=>\tInput data ", out_tuning_data, " is missing", "\n")
    else
        cat("Input data is ", out_tuning_data, "\n")
    
    input_ok = (input_ok && this_input_ok)
    if (input_ok)
        cat("\t=>\tResults for stability selection will be written to ", out_stable_networks, "\n\n")
    }


if(runNetGSA)
    {
    if (SS || BIC)
       stop("=>\tCannot choose two analysis options ")
    cat("\nRunning NetGSA\n")
    
    this_input_ok = file.exists(out_stable_networks)
    if (! this_input_ok)
        stop(paste0("=>\tInput data ",	out_tuning_data, " is missing", "\n\n"))
    else
        cat("\t=>\tInput data is ", out_stable_networks, "\n")
    
    input_ok = (input_ok && this_input_ok)
    
    this_input_ok	= file.exists(out_filtered_rdata)
    if (! this_input_ok)
        stop(paste0(" =>\tInput data ", out_filtered_rdata, " is missing", "\n\n"))
    else
        cat("\t=>\tInput data is ", out_filtered_rdata, "\n")
    
    input_ok = (input_ok && this_input_ok)
    
    if (input_ok)
        {
        cat("\n\t=>\tThresholded edgelist will be written to ", out_edgelist, "\n")
        cat("\t=>\tFull edgelist will be written to ", out_edgelist_full, "\n")
        cat("\t=>\tNodelist will be written to ", out_nodelist, "\n")
        cat("\t=>\tClusters pdf will be called  ", out_clusters_pdf, "\n")
        cat("\t=>\tNetGSA rdata will be written to ", out_netgsa_rda, "\n")
        cat("\t=>\tNetGSA csv will be written to ", out_netgsa, "\n")
        }
    }

cat("\n\n\n\n")


#########################  Start analysis ##################################

if (SS || BIC)
  {
  ## Load data and prepare for joint estimation. Note in this file, data is organized such that the first column is group information, the column names represent variables/CRIC
  cat("Input data is " , input_data, "\n")
  dat <- read.csv(input_data, row.names = 1, check.names = FALSE, as.is = TRUE, na.strings = c("NA", ""))
  
  if (! "GROUP" %in% names(dat))
      stop("Your data does not contain a column with header 'GROUP' indicating group membership")

   sample_info <- dat$GROUP
   metabs <- colnames(dat)[-1]


   if (length(metabs) != length(unique(metabs)))
      stop("One or more metabolites are duplicated in your input file. Metabolites must be unique")

   if (length(unique(sample_info)) < 2)
      stop("Your input file only contains 1 group -- differential analysis requires two groups")

   if (length(unique(sample_info)) > 2)
      cat("Warning : More than two groups are defined in your input file - only the first two will be analysed")
  
  dat <- dat[,-1]
  dat <- as.matrix(dat)
  rownames(dat) <- NULL
  colnames(dat) <-NULL
  
  dataset <- list()
  dataset$metab_info <- data.frame("Compound"=metabs, stringsAsFactors=FALSE)
  dataset$sample_info <- data.frame("ID"=seq(1:length(sample_info)),"GROUP"=sample_info, stringsAsFactors=FALSE)
  dataset$dat <- t(dat)
  temp3 <- sapply(dataset$metab_info$Compound, function(v) unlist(strsplit(as.character(v)[1], split="___"))[[1]])
  names(temp3) <- NULL
  dataset$metab_info$ShortName <- temp3
  
  table(dataset$sample_info$GROUP)
  p <- nrow(dataset$metab_info)
  ncond <- length(unique(dataset$sample_info$GROUP))
  cond_names <- names(table(dataset$sample_info$GROUP))
  
  dat <- vector("list", ncond)
  for (loop_cond in 1:ncond)
     {
     dat[[loop_cond]] <- dataset$dat[,(dataset$sample_info$GROUP==cond_names[loop_cond])]
     }
  
  #if (pearsonFiltering)
  #{
  # ## Next check correlations among variables. Remove those that have correlations lower than the specified threshold.
  # print("Screening correlation...")
  # dat <- lapply(dat, function(d) t(scale(t(d))))
  # r <- lapply(dat, function(v) cor(t(v)))
  # drop_index <- vector("list", ncond)
  # for (loop_cond in 1:ncond)
  #     {
  #     print(loop_cond)
  #     r[[loop_cond]][which(abs(r[[loop_cond]]) < cor.threshold)] <- 0
  #     r[[loop_cond]][which(abs(r[[loop_cond]]) >= cor.threshold)] <- 1
  #     diag(r[[loop_cond]]) <- 0
  #     drop_index[[loop_cond]] <- (apply(r[[loop_cond]], 1, sum) == 0)
  #     }
  
  #    #Drop variables if needed
  #   if ( max(Reduce("+", drop_index)) == ncond)
  #    {
  #    dataset <- drop.metabs(dataset, (Reduce("+", drop_index) == ncond))
  #    # correction/change TBD 02/21/18
  #    dim(dataset$dat)
  #    }
  
  #  for (loop_cond in 1:ncond)
  #    {
  #    dat[[loop_cond]] <- dataset$dat[,(dataset$sample_info$GROUP==cond_names[loop_cond])]
  #    }
  #p <-nrow(dataset$metab_info)
  #cat('Filtering by correlation complete\n')
  #}

  
  dataset$metab_info$foldchange <- rowMeans(dat[[2]]) - rowMeans(dat[[1]])
  dataset$metab_info$fcdirection <- sapply(1:p, function(i) ifelse(dataset$metab_info$foldchange[i]>0, "Up", "Down"))
  dataset$metab_info$fc.notes <- "Group 2 over Group 1"
  dataset$metab_info$statistic <- sapply(1:p, function(i) t.test(dat[[2]][i,], dat[[1]][i,], var.equal=FALSE)$statistic)
  dataset$metab_info$pvalue <- sapply(1:p, function(i) t.test(dat[[2]][i,], dat[[1]][i,], var.equal=FALSE)$p.value)
  dataset$metab_info$qvalue <- p.adjust(dataset$metab_info$pvalue, "BH")
  dataset$metab_info$DEstatus <- sapply(1:p, function(i) ifelse(abs(dataset$metab_info$qvalue[i])>=0.05, FALSE, TRUE))

  #IO : Filtered data to be used with NetGSA
  save(dataset, file=out_filtered_rdata)
  
  ## Joint estimation
  dat <- lapply(dat, function(d) t(scale(t(d)))) 
  metab_info <- dataset$metab_info
  n4cov <- max(sapply(dat, ncol))
  trainX <- t(do.call(cbind, dat))
  trainY <- c(rep(1, sum((dataset$sample_info$GROUP==cond_names[1]))), rep(2, sum((dataset$sample_info$GROUP==cond_names[2]))))
  p <- ncol(trainX)
  
  # Tuning parameter -- select (if (BIC)) or load (if SS))
  lambda.guo <- seq(0.01, 0.3, 0.02)*sqrt(log(p)/n4cov)
  
  if (BIC)
    {
    cat("BIC using Guo et al ... \n")
  
    cl <- makeCluster(nCores)
    registerDoParallel(cl)
    bic.guo <- foreach(i=1:length(lambda.guo), .packages = c("MASS", "glasso")) %dopar%
     CGM_AHP_tune(trainX, testX=trainX, model=trainY, lambda=lambda.guo[i], BIC=TRUE, eta=0.1)
    stopCluster(cl)
    
    lastar.guo <- lambda.guo[which.min(sapply(bic.guo, function(a) a$BIC))]
    #IO: Save tuning data
    save(bic.guo, file=out_tuning_data)
    }
  else if (SS)
    {
    #IO: Load tuning data
    load(out_tuning_data)
    
    tmp <- sapply(bic.guo, function(a) a$BIC)
    if (max(is.infinite(tmp))==1)
      {
      bic.guo <- bic.guo[is.finite(tmp)]
      lambda.guo <- lambda.guo[is.finite(tmp)]
      lastar.guo <- lambda.guo[which.min(sapply(bic.guo, function(a) a$BIC))]
      }
    else
       {
       lastar.guo <- lambda.guo[which.min(sapply(bic.guo, function(a) a$BIC))]
       }
    }
  
  
  if (SS)
    {
    listX <- lapply(dat, t)
    
    my.iter <- function(iter, seed.base)
        {
        #Run subsampling and estimation 5 times.
        fit <- CGM_AHP_stabsel(X=listX, cnt=nreps, lastar=lastar.guo, seed.base=seed.base, eta=0.1)
        return(fit)
        }
    
    cat("Stability selection with Guo et al ... \n")
    #Use multiple cores to run the function my.iter nreps = 10 times. In total we get 10*10 subsampling for stability selection
    cl <- makeCluster(nCores)
    registerDoParallel(cl)
    stab_guo = foreach(i = 1:nCores,.packages = c("MASS", "glasso")) %dopar% my.iter(i,i*100)
    stopCluster(cl)
    
    #IO: Save stable networks
    save(stab_guo, file=out_stable_networks)
    }
  }

## Get the partial correlations
if (runNetGSA)
    {
    #IO: Load stable networks and filtered data
    load(out_stable_networks)
    load(out_filtered_rdata)

    ncond <- length(unique(dataset$sample_info$GROUP))
    sel_mat <- vector("list", ncond)
    for (k in 1:ncond)
        {
        sel_mat[[k]] <- lapply(stab_guo, function(r) r$mat[[k]])
        sel_mat[[k]] <- Reduce("+", sel_mat[[k]])
        sel_mat[[k]] <- sel_mat[[k]]/(2 * nCores * nreps)
        }
  
    ###*********************************************###
    ## Estimate the partial correlation matrix
    ###*********************************************###
    n <- ncol(dataset$dat)
    p <- nrow(dataset$dat)
    x <- dataset$dat
    y <- dataset$sample_info$GROUP
    xx <- vector("list", ncond)
    xx[[1]] <- x[, which(y == 0)]
    xx[[2]] <- x[, which(y == 1)]
    Ip <- diag(rep(1,p))
 
    ## Model selection is done via adjusted DGlasso, where the inverse frequency weighted graphical lasso is applied.
    ## FDR 20% is applied to select the significant partial correlations
    wAdj <- vector("list", ncond)
    Qmat <- vector("list", ncond)
    pCorMat <- vector("list", ncond)
    for (k in 1:ncond)
        {
        cat('Estimating model ...', k, '...\n')
        fit <- adjDGlasso(t(xx[[k]]), weights=1/(0.1 + sel_mat[[k]]), alpha = fdr.cutoff)
        pCorMat[[k]] <- fit$Theta
        Qmat[[k]] <- fit$qvalue
        wAdj[[k]] <- fit$Theta * fit$qvalue.fdr
        }
  
    ###*********************************************###
    ##              Output full edge list
    ###*********************************************###
    pairs <- combn(as.character(dataset$metab_info$Compound), 2, simplify=FALSE)
    
    df <- data.frame(Metabolite.A=rep(0,length(pairs)), Metabolite.B=rep(0,length(pairs)),
       pcor.0=rep(0,length(pairs)), pcor.1=rep(0,length(pairs)), qval.0=rep(0,length(pairs)),
       qval.1=rep(0,length(pairs)))
    df[,1:2] <- do.call(rbind, pairs)
    df[,3] <- lowerTriangle(pCorMat[[1]])
    df[,4] <- lowerTriangle(pCorMat[[2]])
    df[,5] <- lowerTriangle(Qmat[[1]])
    df[,6] <- lowerTriangle(Qmat[[2]])
  
    #IO: Full edgelist
    write.csv(df, file=out_edgelist_full, row.names=FALSE, quote = FALSE)

    ## Get the adjacency matrix by thresholding the partial correlations
    Ahat <- NULL
    for (k in 1:ncond)
        {
        Ahat[[k]] <- abs(wAdj[[k]]) >= matrix(rep(eps, p^2), p, p)
        }
  
    cat("Number of edges in Group 1: ", sum(Ahat[[1]])/2, "\n")
    cat("Number of edges in Group 2: ", sum(Ahat[[2]])/2, "\n")
  
    ###*********************************************###
    ##      Output FDR thresholded edge list
    ###*********************************************###
    pairs <- combn(as.character(dataset$metab_info$Compound), 2, simplify=FALSE)
    df <- data.frame(Metabolite.A=rep(0,length(pairs)), Metabolite.B=rep(0,length(pairs)),
           pcor.0=rep(0,length(pairs)), pcor.1=rep(0,length(pairs)), qval.0=rep(0,length(pairs)),
           qval.1=rep(0,length(pairs)))
    df[,1:2] <- do.call(rbind, pairs)
    df[,3] <- lowerTriangle(wAdj[[1]])
    df[,4] <- lowerTriangle(wAdj[[2]])
    df[,5] <- lowerTriangle(Qmat[[1]])
    df[,6] <- lowerTriangle(Qmat[[2]])
    df$edge <- rep(-99, length(pairs))#non-edge
    df$edge[which((abs(df$pcor.0) >= eps)*(abs(df$pcor.1) >= eps)==1)] <- "Both" #common edge
    df$edge[which((abs(df$pcor.0) >= eps)*(abs(df$pcor.1)  < eps)==1)] <- "Group 1"
    df$edge[which((abs(df$pcor.0) <  eps)*(abs(df$pcor.1) >= eps)==1)] <- "Group 2"
    df <- df[(df$edge!=-99),]
    rownames(df) <- NULL
  
    # IO: Thresholded edgelist
    write.csv(df, file=out_edgelist, row.names=FALSE, quote = FALSE)
  
    ## Join the two networks
    myGraph <- vector("list", length(Ahat))
    for (loop_el in 1:length(Ahat))
        {
        g <- graph_from_adjacency_matrix(wAdj[[loop_el]], mode="undirected", weighted = TRUE)
        V(g)$name <- as.character(dataset$metab_info$ShortName)
        myGraph[[loop_el]] <- g
        }
  
    jointGraph <- igraph::union(myGraph[[1]], myGraph[[2]])
    jointLayout <- layout_nicely(jointGraph)
    E(jointGraph)$lty <- 1
    E(jointGraph)$color <- "black"
    E(jointGraph)$lty[is.na(E(jointGraph)$weight_2)] <- 2 ## 0: non-progressor
    E(jointGraph)$lty[is.na(E(jointGraph)$weight_1)] <- 3 ## 1: progressor
    E(jointGraph)$color[is.na(E(jointGraph)$weight_2)] <- "green" ## 0: non-progressor
    E(jointGraph)$color[is.na(E(jointGraph)$weight_1)] <- "red" ## 1: progressor
    V(jointGraph)$color <- ifelse(dataset$metab_info$DEstatus=="TRUE", "purple", "white")
    V(jointGraph)$DE <- dataset$metab_info$DEstatus
  
    if (runNetGSA)
        {
        ###*********************************************###
        ###        Ensemble community detection         ###
        ###*********************************************###
        fit <- run_consensus_cluster(jointGraph,tau=tau0,method="ensemble")
        consensus_membership <- fit$dcl
        B <- matrix(0, nrow=length(unique(consensus_membership)), p)
        rownames(B) <- paste0("Cluster",1:length(unique(consensus_membership)))
        for (j in 1:nrow(B))
            {
            B[j,which(consensus_membership==j)] <- 1
            }
        if (length(which(rowSums(B)<5))>0)
            {
            B <- B[-which(rowSums(B)<5),]
            }
        npath <- nrow(B)
    
        summary_list <- list()
        for (loop_cluster in 1:nrow(B) )
            {
            cluster_c <- induced.subgraph(jointGraph, V(jointGraph)$name[(B[loop_cluster,]==1)])
            summary_list[[loop_cluster]] <- data.frame("number.of.nodes"=length(V(cluster_c)),
              "number.of.edges"=length(E(cluster_c)),
              "number.of.DE.nodes"=sum(as.numeric(table(V(cluster_c)$DE)[names(table(V(cluster_c)$DE))==TRUE])),
              "number.of.DE.edges"=sum(as.numeric(table(E(cluster_c)$color)[names(table(E(cluster_c)$color)) %in% c("red", "green")])))
            }
        
        summary_stat <- data.frame("cluster.name"= rownames(B), do.call(rbind, summary_list))
    
    
        ###*********************************************###
        ##                   Run NetGSA
        ###*********************************************###
        out.netgsa <- NetGSA(wAdj, x = cbind(xx[[1]], xx[[2]]), y = c(rep(1, ncol(xx[[1]])), rep(2, ncol(xx[[2]]))), B = B, lklMethod = "REML")
        dataset$metab_info$meanchange <- out.netgsa$beta[[2]] - out.netgsa$beta[[1]]
        dataset$metab_info$mc.notes <- "Group 2 over Group 1"
    
    
        ###****************************************************************###
        ##        Output enrichment results -- clusters sorted by significance
        ###****************************************************************###

        res <- data.frame(summary_stat, "NetGSA-pval"=out.netgsa$p.value, "NetGSA-pFDR"=p.adjust(out.netgsa$p.value, "BH"))
        res <- res[order(res$NetGSA.pFDR),]
        rownames(res) <- 1:nrow(res)
        
        dataset$metab_info$membership <- consensus_membership
        dataset$metab_info$membership[!(dataset$metab_info$membership %in% gsub('Cluster','',res$cluster.name))] <- NA
        dataset$metab_info$membership <- rownames(res)[match(dataset$metab_info$membership, as.numeric(gsub('Cluster','',res$cluster.name)))]
        dataset$metab_info$membership <- as.numeric(dataset$metab_info$membership)
        res$cluster.name <- paste0("Cluster",rownames(res))
    
        # IO: netgsa results (csv + binary) and nodelist
        write.csv(res, file=out_netgsa, row.names=FALSE)
        write.csv(dataset$metab_info, file=out_nodelist, row.names=FALSE, quote = FALSE)
        save.image(file=out_netgsa_rda)

        ## Output the clusters
        if (savePlots)
            {
            ## Size nodes by mean change
            V(jointGraph)$size <- abs(dataset$metab_info$meanchange)*30
            
            #IO: PDF graphs
            pdf(out_clusters_pdf, height = 10, width = 10)
            
            for (loop_cluster in 1:nrow(B) )
                {
                cluster_c <- induced.subgraph(jointGraph, V(jointGraph)$name[which(dataset$metab_info$membership==loop_cluster)])
                plot(cluster_c, vertex.label = V(cluster_c)$name, vertex.label.cex = 1,
                  layout = layout.fruchterman.reingold(cluster_c),
                  main = paste0(res$cluster.name[loop_cluster]," ( qvalue - ",round(res$NetGSA.pFDR[loop_cluster],2),")"))
                }
            dev.off()
            }
        }
    }



