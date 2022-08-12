if (FALSE) {
    Rcpp::compileAttributes("~/Dropbox/Academic/Papers/004_FCM/Code/treeSB-master/", verbose = TRUE)
    remove.packages("treeSB")
    unloadNamespace("treeSB")
    # 3. Build the changes via R CMD build treeSB-master/ 
    #    which will result in a new archive file named treeSB_0.1.0.tar.gz 
    # 4. Install your modified archive via R CMD INSTALL treeSB_0.1.0.tar.gz 
    # 5. Close and open RStudio. 
}

print("start of example.R")
rm(list=ls())

tstr <- "BT"

n <- 100
Y <- rnorm(n)
C <- rep(1, n)
J <- length(unique(C))
covX <- rep(c(0,1), n/2)
psiX <- matrix(1, nrow=n, ncol=1)  
psiX <- cbind(psiX, covX) 
R <- ncol(psiX)


# --------------------------
# Train model data ---------
# --------------------------

pmc <- list(npart=10, nburn=5, nsave=200, nskip=1, ndisplay=10)  
# Model parameters:
prior <- list(K=8, merge_par=0, zeta=1.0, merge_step=F, use_skew=T) 
prior$gam_Sig <- Matrix::diag(nrow=R) * 100  
prior$treestr <- ifelse(tstr=="BT", 1, 0)  # 0 is LT; 1 is BT

tag <- paste0("n", n, "_", tstr, "_K", prior$K)

try({
    # Fit model:
    state <- list(t=stats::kmeans(Y, prior$K/2, iter.max = 100)$cluster - 1)
    res <- treeSB::fit_treeSB(Y, psiX, C, prior=prior, pmc=pmc, state=state) 
    saveRDS(res, paste0("res_", tag, ".rds"))
    res$data$Y <- Y
    res$data$psiX <- psiX 
    res$data$C <- C
    
    # Handle label switching:
    resRelab <- treeSB::relabelChain(res) 
    saveRDS(resRelab, paste0("rr_", tag, ".rds"))
    
})  # end try

print("end of example.R")