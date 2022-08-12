rm(list=ls())

# --------------------------
# User-set settings --------
# --------------------------

# Set to TRUE to save memory so R doesn't crash as often
save_memory <- TRUE
# Use covariate information?
use_age  <- TRUE
use_sex  <- TRUE
use_race <- TRUE
use_bmi  <- FALSE

to_debug <- TRUE  # Want to run quickly?
tstr <- "BT"  # or "LT"
zeta_seq <- 0.1  # a vector of values between 0 and 1
myK <- 64  # usually 64, but need to increase if (zeta > 0.15 and n=443520)

tag <- "grif_manysubj_"
if (use_age)  tag <- paste0(tag, "age_")   # age group 18-29, 30-39, 40-49, 50-69
if (use_sex)  tag <- paste0(tag, "sex_")   # sex - male 0
if (use_race) tag <- paste0(tag, "race_")  # race - white 0
if (use_bmi)  tag <- paste0(tag, "bmi_")   # BMI group 18-24, 25-29, 30-35

tag <- paste0(tag, tstr, "_mp", 0.0, "_zeta", zeta_seq)




# --------------------------
# Read data ----------------
# --------------------------

# 40320 201600 443520
raw_data = readRDS("data_grifols_manysubj_n201600.rds")
str(raw_data)
n <- nrow(raw_data$Y)
inds_tr <- 1:n
tag <- paste0(tag, "_ntr", n)

Y <- raw_data$Y[inds_tr, ]
scaleY <- apply(abs(Y), 2, max) * 0.3  # To avoid numerical instability
p <- ncol(Y)
for (i in 1:p) {
    Y[, i] <- Y[, i] / scaleY[i]
}
X <- raw_data$X[inds_tr, ]
covX <- raw_data$covX[inds_tr, ]
C <- raw_data$C[inds_tr]
C_uniq <- unique(C)
J <- length(C_uniq)
C <- sapply(C, function(x) which(C_uniq == x))  # shift elements in C to 1:whatever

inds_covX <- c(rep(use_age, 3), use_race, use_sex)
covX <- covX[, inds_covX]







# --------------------------
# Train model data ---------
# --------------------------

# library(COMIX)

# Range of tuning parameter zeta:
# Fit model for each zeta:
for (zeta in zeta_seq) {
    print(paste0("zeta = ", zeta))
    
    if (to_debug) {
        pmc <- list(npart=10, nburn=10, nsave=5, nskip=2, ndisplay=10)  
    } else {
        pmc <- list(npart=10, nburn=10000, nsave=1000, nskip=50, ndisplay=100)  
    }
    
    # Model parameters:
    prior <- list(scaleY=scaleY, K=myK, merge_par=0, zeta=zeta, merge_step=F, use_skew=T) 
    
    psiX <- matrix(1, nrow=n, ncol=1)  
    psiX <- cbind(psiX, covX) 
    if (save_memory) {
        rm(covX)
        rm(X)
        gc()
    }
    R <- ncol(psiX)
    prior$gam_Sig <- Matrix::diag(nrow=R) * 100  # 10
    prior$treestr <- ifelse(tstr=="BT", 1, 0)  
    seed <- as.numeric(zeta * 101)  # zeta * 100
    set.seed(seed)
    try({
        # Fit model:
        print("before treeSB::fit_treeSB()")
        # browser()
        state <- list(t=stats::kmeans(Y, prior$K/2, iter.max = 100)$cluster - 1)
        start_time <- Sys.time()
        res <- treeSB::fit_treeSB(Y, psiX, C, prior = prior, pmc = pmc, state=state)
        end_time <- Sys.time()
        print("end_time - start_time")
        print(end_time - start_time)
        # saveRDS(res, paste0("res_", tag, ".rds"))
        res$data$Y <- Y
        res$data$psiX <- psiX      
        res$data$C <- C
        
        # Handle label switching:
        print("before treeSB::relabelChain()")  
        resRelab <- treeSB::relabelChain(res)  
        if (save_memory) { 
            rm(res)
            gc()
        }
        resRelab$inds_tr <- inds_tr
        saveRDS(resRelab, paste0("rr_", tag, ".rds"))
        # source("../count_clusters.R")
    })
    
}