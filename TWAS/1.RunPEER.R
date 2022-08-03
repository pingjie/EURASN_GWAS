set.seed(1024)
library(peer)

args <- commandArgs(trailingOnly = T)

categ <- args[1]
tissue <- args[2]

workDir <- paste0("exp/", categ, "/", tissue)
covDir <- paste0("cov/", categ, "/")

### function load Rdata
loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
}

print(paste0("Working on ", categ, ", ", tissue))

if (file.exists(workDir)) {
    expTab <- read.table(paste0(workDir, "/", tissue, "_RNASeq_RMA_Inverse.txt"), sep = "\t", header = T, row.names = 1, as.is = T)
    covTab <- read.table(paste0(covDir, tissue, "_Cov.", categ, ".txt"), sep = "\t", header = T, row.names = 1, as.is = T)

    model <- PEER()
    PEER_setPhenoMean(model, as.matrix(expTab))
    PEER_setCovariates(model, as.matrix(covTab))

    peerN <- 15

    # the num here need to be decided per no. of subjects in the tissue of interest, see https://www.gtexportal.org/home/documentationPage
    # for eQTL: 15 factors for N<150, 30 factors for 150<= N <250, 45 factors for 250<= N <350, and 60 factors for N>=350
    # as a result of optimizing for the number of eGenes discovered

    if ( nrow(PEER_getPhenoMean(model)) < 150) {
        peerN <- 15
    } else if ( nrow(PEER_getPhenoMean(model)) < 250 & nrow(PEER_getPhenoMean(model)) >= 150 ) {
        peerN <- 30
    } else if ( nrow(PEER_getPhenoMean(model)) < 350 & nrow(PEER_getPhenoMean(model)) >= 250 ) {
        peerN <- 45
    } else if ( nrow(PEER_getPhenoMean(model)) >= 350 ) {
        peerN <- 60
    }

    PEER_setNk(model, peerN)

    PEER_update(model)  # 100+ iterations; take 50+ hours
    factors <- PEER_getX(model)
    residuals <- PEER_getResiduals(model)  # This is the Residuals after adjusting for PEER factors and PC variables, which are to be used for following analyses

    save(factors, file = paste0(workDir, "/", tissue, "_", categ, "_PEER_factors.rda"))
    write.table(factors, file = paste0(workDir, "/", tissue, "_", categ, "_PEER_factors.txt"), sep = "\t", row.names = T, quote = F)

    pdf(paste0(workDir, "/diagnostics_peer.GTEx.", tissue, ".", categ, ".pdf"))
    PEER_plotModel(model)
    dev.off()

    rownames(residuals) <- rownames(expTab)
    colnames(residuals) <- colnames(expTab)

    save(residuals, file = paste0(workDir, "/", tissue, "_", categ, "_residuals_gtex.rda"))
    write.table(residuals, file = paste0(workDir, "/", tissue, "_", categ, "_RNASeq_Normal_Inverse_PEER_Corrected_Exp.txt"), sep = "\t", row.names = T, quote = F)

    # Run Inverse quantile transformation again
    residuals_update <- apply(residuals, 2, function(x) qnorm(rank(x, ties.method = "r")/(length(x) + 1)))

    save(residuals_update, file = paste0(workDir, "/", tissue, "_", categ, "_Residuals_GTEx_Inversed.rda"))

    write.table(residuals_update, file = paste0(workDir, "/", tissue, "_", categ, "_RNASeq_Normal_Inverse_PEER_Inversed_Corrected_Exp.txt"), sep = "\t", row.names = T, quote = F)
} else {
    print(paste0("No data of ", tissue))
}

