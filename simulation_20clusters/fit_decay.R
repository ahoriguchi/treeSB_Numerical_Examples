print("start of fit_decay.R")
rm(list=ls())


# --------------------------
# User-set settings --------
# --------------------------

to_debug <- TRUE  # Want to run quickly?

taskID <- 10  # should be one of 10,11,60,61
seed <- taskID
tstr <- ifelse(floor(taskID/10)%%10 ==6, "BT", "LT")
add_cov_eff <- (taskID%%2 == 1)  # Which dataset to use: covariate-independent or covariate-dependent?
mysig <- 10
tag <- paste0(tstr, "sig", mysig, "_covdep", add_cov_eff, "_ID", taskID)
print(tag)

set.seed(seed)

# --------------------------
# Read data ----------------
# --------------------------
raw_data <- readRDS("dat_decay_v3.rds")
str(raw_data)

Y <- raw_data$Y
n <- nrow(Y)
C <- rep(1, n)
J <- length(unique(C))
Y <- apply(Y, 2, function(cc) cc / max(abs(cc)) * 10)
X <- raw_data$Xindep
if (add_cov_eff) {
    X <- raw_data$Xdep
}
covX <- X  # X is either 0 or 1 already


# --------------------------
# Train model data ---------
# --------------------------

pmc <- list(npart=10, nburn=100000, nsave=10000, nskip=50, ndisplay=100) 
if (to_debug) {  # to debug
    pmc <- list(npart=10, nburn=50, nsave=200, nskip=1, ndisplay=10)  
}
# Model parameters:
prior <- list(K=32, merge_par=0, zeta=1.0, merge_step=F, use_skew=T) 
psiX <- matrix(1, nrow=n, ncol=1)  
psiX <- cbind(psiX, covX) 
R <- ncol(psiX)
prior$gam_Sig <- Matrix::diag(nrow=R) * mysig  
prior$treestr <- ifelse(tstr=="BT", 1, 0)  # 0 is LT; 1 is BT

try({
    # Fit model:
    state <- list(t=stats::kmeans(Y, prior$K/2, iter.max = 100)$cluster - 1)
    res <- treeSB::fit_treeSB(Y, psiX, C, prior = prior, pmc = pmc, state=state)  # COMIX::
    saveRDS(res, paste0("res_", tag, "_K", prior$K, ".rds"))
    res$data$Y <- Y
    res$data$psiX <- psiX 
    res$data$C <- C
    
    # Handle label switching:
    resRelab <- treeSB::relabelChain(res)  # COMIX::
    saveRDS(resRelab, paste0("rr_", tag, "_K", prior$K, ".rds"))
    
    # Compute Jaccard distances
    source(file.path("..", "jaccard.R"))
    jm <- JaccardMat(resRelab$chain$t, raw_data$Z)
    saveRDS(jm, paste0("jm_", taskID, ".rds"))
    
})  # end try

print("end of fit_decay.R")