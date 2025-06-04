rm(list = ls())
set.seed(2024)

setwd("C:/Trang/KIProjects/ComprehensionDR/ResultPaper")

# Install required packages if not already installed
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("reshape2")) install.packages("reshape2")
if (!require("pheatmap")) install.packages("pheatmap")

# Load libraries
library(ggplot2)
library(reshape2)
library(pheatmap)


# 1. NMFS
NMF_HDVD_AUC  = read.csv("NMF/NMFS_AUC_HDVD_25runs_rerun.csv" , header = FALSE, skip = 1)
NMF_LAGCN_AUC  = read.csv("NMF/NMFS_AUC_LAGCN_25runs_rerun.csv" , header = FALSE, skip = 1)
NMF_Fdata_AUC  = read.csv("NMF/NMFS_AUC_Fdata_25runs_rerun.csv" , header = FALSE, skip = 1)
NMF_Cdata_AUC  = read.csv("NMF/NMFS_AUC_Cdata_25runs_rerun.csv" , header = FALSE, skip = 1)
NMF_LRSSL_AUC  = read.csv("NMF/NMFS_AUC_LRSSL_25runs_rerun.csv" , header = FALSE, skip = 1)
NMF_Ydata_AUC  = read.csv("NMF/NMFS_AUC_Ydata_25runs_rerun.csv" , header = FALSE, skip = 1)
NMF_oMat_AUC  = read.csv("NMF/NMFS_AUC_oMat_25runs_rerun_neworder.csv" , header = FALSE, skip = 1)
NMF_hsdn_AUC  = read.csv("NMF/NMFS_AUC_hsdn_25runs_rerun_neworder.csv" , header = FALSE, skip = 1)

NMF_HDVD_AUPR  = read.csv("NMF/NMFS_AUPR_HDVD_25runs_rerun.csv" , header = FALSE, skip = 1)
NMF_LAGCN_AUPR  = read.csv("NMF/NMFS_AUPR_LAGCN_25runs_rerun.csv" , header = FALSE, skip = 1)
NMF_Fdata_AUPR  = read.csv("NMF/NMFS_AUPR_Fdata_25runs_rerun.csv" , header = FALSE, skip = 1)
NMF_Cdata_AUPR  = read.csv("NMF/NMFS_AUPR_Cdata_25runs_rerun.csv" , header = FALSE, skip = 1)
NMF_LRSSL_AUPR  = read.csv("NMF/NMFS_AUPR_LRSSL_25runs_rerun.csv" , header = FALSE, skip = 1)
NMF_Ydata_AUPR  = read.csv("NMF/NMFS_AUPR_Ydata_25runs_rerun.csv" , header = FALSE, skip = 1)
NMF_oMat_AUPR  = read.csv("NMF/NMFS_AUPR_oMat_25runs_rerun_neworder.csv" , header = FALSE, skip = 1)
NMF_hsdn_AUPR  = read.csv("NMF/NMFS_AUPR_hsdn_25runs_rerun_neworder.csv" , header = FALSE, skip = 1)

#2. NMFPer
NMFPer_HDVD_AUC = read.csv("NMF-PDR/NMFPer_AUC_HDVD_Wiltest_25.csv", header = FALSE, skip = 1)
NMFPer_LAGCN_AUC_24 = read.csv("NMF-PDR/NMFPer_AUC_wiltest_LAGCN_25runs.csv", header = FALSE, skip = 1)
NMFPer_LAGCN_AUC = rbind(NMFPer_LAGCN_AUC_24, 0.7073)

NMFPer_HDVD_AUPR = read.csv("NMF-PDR/NMFPer_AUPR_HDVD_Wiltest_25.csv", header = FALSE, skip = 2)
NMFPer_LAGCN_AUPR_24 = read.csv("NMF-PDR/LAGCN/NMFPer_AUPR_wiltest_LAGCN_25runs.csv", header = FALSE, skip = 1)
NMFPer_LAGCN_AUPR = rbind(NMFPer_LAGCN_AUPR_24, 0.2381)


load("NMF-PDR/NMFPer_AUC_25_wil_Fdata.RData")
NMFPer_Fdata_AUC = AUC
load("NMF-PDR/NMFPer_AUC_25_wil_Cdata.RData")
NMFPer_Cdata_AUC = AUC
load("NMF-PDR/NMFPer_AUC_25_wil_LRSSL.RData")
NMFPer_LRSSL_AUC = AUC
load("NMF-PDR/NMFPer_AUC_25_wil_Ydata.RData")
NMFPer_Ydata_AUC = AUC
load("NMF-PDR/NMFPer_AUC_25_wil_oMat_neworder.RData")
NMFPer_oMat_AUC = AUC
load("NMF-PDR/NMFPer_AUC_25_wil_hsdn_neworder.RData")
NMFPer_hsdn_AUC = AUC

load("NMF-PDR/NMFPer_AUPR_25_wil_Fdata.RData")
NMFPer_Fdata_AUPR = AUPR
load("NMF-PDR/NMFPer_AUPR_25_wil_Cdata.RData")
NMFPer_Cdata_AUPR = AUPR
load("NMF-PDR/NMFPer_AUPR_25_wil_LRSSL.RData")
NMFPer_LRSSL_AUPR = AUPR
load("NMF-PDR/NMFPer_AUPR_25_wil_Ydata.RData")
NMFPer_Ydata_AUPR = AUPR
load("NMF-PDR/NMFPer_AUPR_25_wil_oMat_neworder.RData")
NMFPer_oMat_AUPR = AUPR
load("NMF-PDR/NMFPer_AUPR_25_wil_hsdn_neworder.RData")
NMFPer_hsdn_AUPR = AUPR

# 3. BNNR
BNNR_HDVD_AUC = read.csv("BNNR/BNNR_HDVD_AUC_25runs_rerun.csv", header = FALSE)
BNNR_LAGCN_AUC = read.csv("BNNR/BNNR_LAGCN_AUC_25runs_rerun.csv", header = FALSE )
BNNR_Fdata_AUC = read.csv("BNNR/BNNR_Fdata_AUC_25runs_rerun.csv", header = FALSE )
BNNR_Cdata_AUC = read.csv("BNNR/BNNR_Cdata_AUC_25runs_rerun.csv", header = FALSE )
BNNR_LRSSL_AUC = cbind(read.csv("BNNR/BNNR_LRSSL_AUC_9runs_rerun.csv", header = FALSE ),
                       read.csv("BNNR/BNNR_LRSSL_AUC_16runs_rerun.csv", header = FALSE ))
BNNR_Ydata_AUC = cbind(read.csv("BNNR/BNNR_Ydata_AUC_5runs_rerun.csv", header = FALSE ),
                       read.csv("BNNR/BNNR_Ydata_AUC_25runs_rerun.csv", header = FALSE ))
BNNR_oMat_AUC = read.csv("BNNR/BNNR_oMat_AUC_25runs_rerun_neworder.csv", header = FALSE )
BNNR_hsdn_AUC = read.csv("BNNR/25runs/BNNR_hsdn_AUC_25runs_rerun_neworder.csv", header = FALSE )



BNNR_HDVD_AUPR = read.csv("BNNR/BNNR_HDVD_AUPR_25runs_rerun.csv", header = FALSE )
BNNR_LAGCN_AUPR = read.csv("BNNR/BNNR_LAGCN_AUPR_25runs_rerun.csv", header = FALSE )
BNNR_Fdata_AUPR = read.csv("BNNR/BNNR_Fdata_AUPR_25runs_rerun.csv", header = FALSE )
BNNR_Cdata_AUPR = read.csv("BNNR/BNNR_Cdata_AUPR_25runs_rerun.csv", header = FALSE )
BNNR_LRSSL_AUPR = cbind(read.csv("BNNR/BNNR_LRSSL_AUPR_9runs_rerun.csv", header = FALSE ),
                        read.csv("BNNR/BNNR_LRSSL_AUPR_16runs_rerun.csv", header = FALSE ))
BNNR_Ydata_AUPR = cbind(read.csv("BNNR/BNNR_Ydata_AUPR_5runs_rerun.csv", header = FALSE ),
                        read.csv("BNNR/BNNR_Ydata_AUPR_25runs_rerun.csv", header = FALSE ))
BNNR_oMat_AUPR = read.csv("BNNR/BNNR_oMat_AUPR_25runs_rerun_neworder.csv", header = FALSE )
BNNR_hsdn_AUPR = read.csv("BNNR/BNNR_hsdn_AUPR_25runs_rerun_neworder.csv", header = FALSE )

## OMC
OMC_HDVD_AUC = read.csv("OMC/OMC_HDVD_AUC_25runs_rerun.csv", header = FALSE )
OMC_LAGCN_AUC = read.csv("OMC/OMC_LAGCN_AUC_25runs_rerun.csv", header = FALSE )
OMC_Fdata_AUC = read.csv("OMC/OMC_Fdata_AUC_25runs_rerun.csv", header = FALSE )
OMC_Cdata_AUC = read.csv("OMC/OMC_Cdata_AUC_25runs_rerun.csv", header = FALSE )
OMC_LRSSL_AUC = read.csv("OMC/OMC_LRSSL_AUC_25runs_rerun.csv", header = FALSE )
OMC_Ydata_AUC = cbind(read.csv("OMC/OMC_Ydata_AUC_25runs_rerun.csv", header = FALSE ),
                      read.csv("OMC/OMC_Ydata_AUC_5runs_rerun.csv", header = FALSE ))
OMC_oMat_AUC = read.csv("OMC/OMC_oMat_AUC_25runs_rerun_neworder.csv", header = FALSE )
OMC_hsdn_AUC = read.csv("OMC/OMC_hsdn_AUC_25runs_rerun_neworder.csv", header = FALSE )


OMC_HDVD_AUPR = read.csv("OMC/OMC_HDVD_AUPR_25runs_rerun.csv", header = FALSE )
OMC_LAGCN_AUPR = read.csv("OMC/OMC_LAGCN_AUPR_25runs_rerun.csv", header = FALSE )
OMC_Fdata_AUPR = read.csv("OMC/OMC_Fdata_AUPR_25runs_rerun.csv", header = FALSE )
OMC_Cdata_AUPR = read.csv("OMC/OMC_Cdata_AUPR_25runs_rerun.csv", header = FALSE )
OMC_LRSSL_AUPR = read.csv("OMC/OMC_LRSSL_AUPR_25runs_rerun.csv", header = FALSE )
OMC_Ydata_AUPR = cbind(read.csv("OMC/OMC_Ydata_AUPR_25runs_rerun.csv", header = FALSE ),
                       read.csv("OMC/OMC_Ydata_AUPR_5runs_rerun.csv", header = FALSE ))

OMC_oMat_AUPR = read.csv("OMC/OMC_oMat_AUPR_25runs_rerun_neworder.csv", header = FALSE )
OMC_hsdn_AUPR = read.csv("OMC/OMC_hsdn_AUPR_25runs_rerun_neworder.csv", header = FALSE )



# 5. HGIMC
HGIMC_HDVD_AUC = read.csv("HGIMC/HGIMC_HDVD_AUC_25runs_rerun.csv",header = FALSE)
HGIMC_LAGCN_AUC = read.csv("HGIMC/HGIMC_LAGCN_AUC_25runs_rerun.csv",header = FALSE)
HGIMC_Fdata_AUC = read.csv("HGIMC/HGIMC_Fdata_AUC_25runs_rerun.csv",header = FALSE)
HGIMC_Cdata_AUC = read.csv("HGIMC/HGIMC_Cdata_AUC_25runs_rerun.csv",header = FALSE)
HGIMC_LRSSL_AUC = read.csv("HGIMC/HGIMC_LRSSL_AUC_25runs_rerun.csv",header = FALSE)
HGIMC_Ydata_AUC = read.csv("HGIMC/HGIMC_Ydata_AUC_25runs_rerun.csv",header = FALSE)
HGIMC_oMat_AUC = read.csv("HGIMC/HGIMC_oMat_AUC_25runs_neworder.csv",header = FALSE)
HGIMC_hsdn_AUC = read.csv("HGIMC/HGIMC_hsdn_AUC_25runs_neworder.csv",header = FALSE)


HGIMC_HDVD_AUPR = read.csv("HGIMC/HGIMC_HDVD_AUPR_25runs_rerun.csv",header = FALSE)
HGIMC_LAGCN_AUPR = read.csv("HGIMC/HGIMC_LAGCN_AUPR_25runs_rerun.csv",header = FALSE)
HGIMC_Fdata_AUPR = read.csv("HGIMC/HGIMC_Fdata_AUPR_25runs_rerun.csv",header = FALSE)
HGIMC_Cdata_AUPR = read.csv("HGIMC/HGIMC_Cdata_AUPR_25runs_rerun.csv",header = FALSE)
HGIMC_LRSSL_AUPR = read.csv("HGIMC/HGIMC_LRSSL_AUPR_25runs_rerun.csv",header = FALSE)
HGIMC_Ydata_AUPR = read.csv("HGIMC/HGIMC_Ydata_AUPR_25runs_rerun.csv",header = FALSE)
HGIMC_oMat_AUPR = read.csv("HGIMC/HGIMC_oMat_AUPR_25runs_neworder.csv",header = FALSE)
HGIMC_hsdn_AUPR = read.csv("HGIMC/HGIMC_hsdn_AUPR_25runs_neworder.csv",header = FALSE)

# 6. VDA
VDA_HDVD_AUC = read.csv("VDA-GKSBMF/VDA_HDVD_AUC_25runs_rerun.csv", header = FALSE)
VDA_LAGCN_AUC = read.csv("VDA-GKSBMF/VDA_LAGCN_AUC_25runs_rerun.csv", header = FALSE)
VDA_Fdata_AUC = read.csv("VDA-GKSBMF/VDA_Fdata_AUC_25runs_rerun.csv", header = FALSE)
VDA_Cdata_AUC = read.csv("VDA-GKSBMF/VDA_Cdata_AUC_25runs_rerun.csv", header = FALSE)
VDA_LRSSL_AUC = read.csv("VDA-GKSBMF/VDA_LRSSL_AUC_25runs_rerun.csv", header = FALSE)
VDA_Ydata_AUC = cbind(read.csv("VDA-GKSBMF/VDA_Ydata_AUC_20runs_rerun.csv", header = FALSE),
                      read.csv("VDA-GKSBMF/VDA_Ydata_AUC_5runs_rerun.csv", header = FALSE))
VDA_oMat_AUC = read.csv("VDA-GKSBMF/VDA_oMat_AUC_25runs_rerun_neworder.csv", header = FALSE)
VDA_hsdn_AUC = read.csv("VDA-GKSBMF/VDA_hsdn_AUC_25runs_rerun_neworder.csv", header = FALSE)


VDA_HDVD_AUPR = read.csv("VDA-GKSBMF/VDA_HDVD_AUPR_25runs_rerun.csv", header = FALSE)
VDA_LAGCN_AUPR = read.csv("VDA-GKSBMF/VDA_LAGCN_AUPR_25runs_rerun.csv", header = FALSE)
VDA_Fdata_AUPR = read.csv("VDA-GKSBMF/VDA_Fdata_AUPR_25runs_rerun.csv", header = FALSE)
VDA_Cdata_AUPR = read.csv("VDA-GKSBMF/VDA_Cdata_AUPR_25runs_rerun.csv", header = FALSE)
VDA_LRSSL_AUPR = read.csv("VDA-GKSBMF/VDA_LRSSL_AUPR_25runs_rerun.csv", header = FALSE)
VDA_Ydata_AUPR = cbind(read.csv("VDA-GKSBMF/VDA_Ydata_AUPR_20runs_rerun.csv", header = FALSE),
                       read.csv("VDA-GKSBMF/VDA_Ydata_AUPR_5runs_rerun.csv", header = FALSE))
VDA_oMat_AUPR = read.csv("VDA-GKSBMF/VDA_oMat_AUPR_25runs_rerun_neworder.csv", header = FALSE)
VDA_hsdn_AUPR = read.csv("VDA-GKSBMF/VDA_hsdn_AUPR_25runs_rerun_neworder.csv", header = FALSE)


# 7. NMFDR

NMFDR_HDVD_AUC = read.csv("NMFDR/NMFDR_HDVD_AUC_25runs_rerun.csv", header = FALSE)
NMFDR_LAGCN_AUC = read.csv("NMFDR/NMFDR_LAGCN_AUC_25runs_rerun.csv", header = FALSE)
NMFDR_Fdata_AUC = read.csv("NMFDR/NMFDR_Fdata_AUC_25runs_rerun.csv", header = FALSE)
NMFDR_Cdata_AUC = read.csv("NMFDR/NMFDR_Cdata_AUC_25runs_rerun.csv", header = FALSE)
NMFDR_LRSSL_AUC = read.csv("NMFDR/NMFDR_LRSSL_AUC_25runs_rerun.csv", header = FALSE)
NMFDR_Ydata_AUC = read.csv("NMFDR/NMFDR_Ydata_AUC_25runs_rerun.csv", header = FALSE)
NMFDR_oMat_AUC = read.csv("NMFDR/NMFDR_oMat_AUC_25runs_rerun_neworder.csv", header = FALSE)
NMFDR_hsdn_AUC = cbind(read.csv("NMFDR/NMFDR_hsdn_AUC_16runs_rerun_neworder.csv", header = FALSE),
                       read.csv("NMFDR/NMFDR_hsdn_AUC_9runs_rerun_neworder.csv", header = FALSE),
                       read.csv("NMFDR/NMFDR_hsdn_AUC_2runs_rerun_neworder.csv", header = FALSE))

NMFDR_HDVD_AUPR = read.csv("NMFDR/NMFDR_HDVD_AUPR_25runs_rerun.csv", header = FALSE)
NMFDR_LAGCN_AUPR = read.csv("NMFDR/NMFDR_LAGCN_AUPR_25runs_rerun.csv", header = FALSE)
NMFDR_Fdata_AUPR = read.csv("NMFDR/NMFDR_Fdata_AUPR_25runs_rerun.csv", header = FALSE)
NMFDR_Cdata_AUPR = read.csv("NMFDR/NMFDR_Cdata_AUPR_25runs_rerun.csv", header = FALSE)
NMFDR_LRSSL_AUPR = read.csv("NMFDR/NMFDR_LRSSL_AUPR_25runs_rerun.csv", header = FALSE)
NMFDR_Ydata_AUPR = read.csv("NMFDR/NMFDR_Ydata_AUPR_25runs_rerun.csv", header = FALSE)
NMFDR_oMat_AUPR = read.csv("NMFDR/NMFDR_oMat_AUPR_25runs_rerun_neworder.csv", header = FALSE)
NMFDR_hsdn_AUPR = cbind(read.csv("NMFDR/NMFDR_hsdn_AUPR_16runs_rerun_neworder.csv", header = FALSE),
                        read.csv("NMFDR/NMFDR_hsdn_AUPR_9runs_rerun_neworder.csv", header = FALSE),
                        read.csv("NMFDR/NMFDR_hsdn_AUPR_2runs_rerun_neworder.csv", header = FALSE))


# 8. IBCF
IBCF_LAGCN_AUC = read.csv("IBCF/RS_IBCF_LAGCN_AUC.csv", header = FALSE, skip = 1)
IBCF_HDVD_AUC = read.csv("IBCF/RS_IBCF_HDVD_AUC.csv", header = FALSE, skip = 1)
IBCF_Fdata_AUC = read.csv("IBCF/RS_IBCF_Fdata_AUC.csv", header = FALSE, skip = 1)
IBCF_Cdata_AUC = read.csv("IBCF/RS_IBCF_Cdata_AUC.csv", header = FALSE, skip = 1)
IBCF_LRSSL_AUC = rbind(read.csv("IBCF/RS_IBCF_LRSSL_AUC_10.csv", header = FALSE, skip = 1),
                       read.csv("IBCF/RS_IBCF_LRSSL_AUC_15.csv", header = FALSE, skip = 1))
IBCF_Ydata_AUC = read.csv("IBCF/RS_IBCF_Ydata_AUC.csv", header = FALSE, skip = 1)
IBCF_oMat_AUC = read.csv("IBCF/25runs/RS_IBCF_oMat_AUC_neworder.csv", header = FALSE, skip = 1)
IBCF_hsdn_AUC = read.csv("IBCF/RS_IBCF_hsdn_AUC_neworder.csv", header = FALSE, skip = 1)

IBCF_LAGCN_AUPR = read.csv("IBCF/RS_IBCF_LAGCN_AUPR.csv", header = FALSE, skip = 1)
IBCF_HDVD_AUPR = read.csv("IBCF/RS_IBCF_HDVD_AUPR.csv", header = FALSE, skip = 1)
IBCF_Fdata_AUPR = read.csv("IBCF/RS_IBCF_Fdata_AUPR.csv", header = FALSE, skip = 1)
IBCF_Cdata_AUPR = read.csv("IBCF/RS_IBCF_Cdata_AUPR.csv", header = FALSE, skip = 1)
IBCF_LRSSL_AUPR = rbind(read.csv("IBCF/RS_IBCF_LRSSL_AUPR_10.csv", header = FALSE, skip = 1),
                        read.csv("IBCF/RS_IBCF_LRSSL_AUPR_15.csv", header = FALSE, skip = 1))
IBCF_Ydata_AUPR = read.csv("IBCF/RS_IBCF_Ydata_AUPR.csv", header = FALSE, skip = 1)
IBCF_oMat_AUPR = read.csv("IBCF/RS_IBCF_oMat_AUPR_neworder.csv", header = FALSE, skip = 1)
IBCF_hsdn_AUPR = read.csv("IBCF/RS_IBCF_hsdn_AUPR_neworder.csv", header = FALSE, skip = 1)



# 9. LIBCF
LIBMF_LAGCN_AUC = read.csv("LIBCF/RS_LIBCF_LAGCN_AUC.csv", header = FALSE, skip = 1)
LIBMF_HDVD_AUC = read.csv("LIBCF/RS_LIBCF_HDVD_AUC.csv", header = FALSE, skip = 1)
LIBMF_Fdata_AUC = read.csv("LIBCF/RS_LIBCF_Fdata_AUC.csv", header = FALSE, skip = 1)
LIBMF_Cdata_AUC = read.csv("LIBCF/RS_LIBCF_Cdata_AUC.csv", header = FALSE, skip = 1)
LIBMF_LRSSL_AUC = read.csv("LIBCF/RS_LIBCF_LRSSL_AUC.csv", header = FALSE, skip = 1)
LIBMF_Ydata_AUC = read.csv("LIBCF/RS_LIBCF_Ydata_AUC.csv", header = FALSE, skip = 1)
LIBMF_oMat_AUC = read.csv("LIBCF/RS_LIBMF_oMat_AUC_neworder.csv", header = FALSE, skip = 1)
LIBMF_hsdn_AUC = read.csv("LIBCF/RS_LIBMF_hsdn_AUC_neworder.csv", header = FALSE, skip = 1)

LIBMF_LAGCN_AUPR = read.csv("LIBCF/RS_LIBCF_LAGCN_AUPR.csv", header = FALSE, skip = 1)
LIBMF_HDVD_AUPR = read.csv("LIBCF/RS_LIBCF_HDVD_AUPR.csv", header = FALSE, skip = 1)
LIBMF_Fdata_AUPR = read.csv("LIBCF/RS_LIBCF_Fdata_AUPR.csv", header = FALSE, skip = 1)
LIBMF_Cdata_AUPR = read.csv("LIBCF/RS_LIBCF_Cdata_AUPR.csv", header = FALSE, skip = 1)
LIBMF_LRSSL_AUPR = read.csv("LIBCF/RS_LIBCF_LRSSL_AUPR.csv", header = FALSE, skip = 1)
LIBMF_Ydata_AUPR = read.csv("LIBCF/RS_LIBCF_Ydata_AUPR.csv", header = FALSE, skip = 1)
LIBMF_oMat_AUPR = read.csv("LIBCF/RS_LIBMF_oMat_AUPR_neworder.csv", header = FALSE, skip = 1)
LIBMF_hsdn_AUPR = read.csv("LIBCF/RS_LIBMF_hsdn_AUPR_neworder.csv", header = FALSE, skip = 1)


# 1 HDVD
median_AUC_HDVD = c(median(as.numeric(NMF_HDVD_AUC$V2)),median(as.numeric(NMFPer_HDVD_AUC$V2)),
                    median(as.numeric(NMFDR_HDVD_AUC)),median(as.numeric(VDA_HDVD_AUC)),
                    median(as.numeric(BNNR_HDVD_AUC)),
                    median(as.numeric(OMC_HDVD_AUC)),median(as.numeric(HGIMC_HDVD_AUC)),
                    median(IBCF_HDVD_AUC$V2), median(LIBMF_HDVD_AUC$V2))

sd_AUC_HDVD = c(sd(as.numeric(NMF_HDVD_AUC$V2)),sd(as.numeric(NMFPer_HDVD_AUC$V2)),sd(as.numeric(NMFDR_HDVD_AUC)),sd(as.numeric(VDA_HDVD_AUC)),
                sd(as.numeric(BNNR_HDVD_AUC)),
                sd(as.numeric(OMC_HDVD_AUC)),sd(as.numeric(HGIMC_HDVD_AUC)), sd(IBCF_HDVD_AUC$V2), sd(LIBMF_HDVD_AUC$V2))

max_AUC_HDVD = c(max(as.numeric(NMF_HDVD_AUC$V2)),max(as.numeric(NMFPer_HDVD_AUC$V2)),max(as.numeric(NMFDR_HDVD_AUC)),max(as.numeric(VDA_HDVD_AUC)),
                 max(as.numeric(BNNR_HDVD_AUC)),
                 max(as.numeric(OMC_HDVD_AUC)),max(as.numeric(HGIMC_HDVD_AUC)), max(IBCF_HDVD_AUC$V2), max(LIBMF_HDVD_AUC$V2))


median_AUPR_HDVD = c(median(as.numeric(NMF_HDVD_AUPR$V2)),median(as.numeric(NMFPer_HDVD_AUPR$V1)),median(as.numeric(NMFDR_HDVD_AUPR)),median(as.numeric(VDA_HDVD_AUPR)),
                     median(as.numeric(BNNR_HDVD_AUPR)),
                     median(as.numeric(OMC_HDVD_AUPR)),median(as.numeric(HGIMC_HDVD_AUPR)), median(IBCF_HDVD_AUPR$V2), median(LIBMF_HDVD_AUPR$V2))

sd_AUPR_HDVD = c(sd(as.numeric(NMF_HDVD_AUPR$V2)),sd(as.numeric(NMFPer_HDVD_AUPR$V1)),sd(as.numeric(NMFDR_HDVD_AUPR)),sd(as.numeric(VDA_HDVD_AUPR)),
                 sd(as.numeric(BNNR_HDVD_AUPR)),
                 sd(as.numeric(OMC_HDVD_AUPR)),sd(as.numeric(HGIMC_HDVD_AUPR)), sd(IBCF_HDVD_AUPR$V2), sd(LIBMF_HDVD_AUPR$V2))


## 2. LAGCN
median_AUC_LAGCN = c(median(as.numeric(NMF_LAGCN_AUC$V2)),
                     median(as.numeric(NMFPer_LAGCN_AUC$V1)),
                     median(as.numeric(NMFDR_LAGCN_AUC)),
                     median(as.numeric(VDA_LAGCN_AUC)),
                     median(as.numeric(BNNR_LAGCN_AUC)),
                     median(as.numeric(OMC_LAGCN_AUC)),
                     median(as.numeric(HGIMC_LAGCN_AUC)),
                     median(IBCF_LAGCN_AUC$V2), median(LIBMF_LAGCN_AUC$V2))


sd_AUC_LAGCN = c(sd(as.numeric(NMF_LAGCN_AUC$V2)),
                 sd(as.numeric(NMFPer_LAGCN_AUC$V1)),
                 sd(as.numeric(NMFDR_LAGCN_AUC)),
                 sd(as.numeric(VDA_LAGCN_AUC)),
                 sd(as.numeric(BNNR_LAGCN_AUC)),
                 sd(as.numeric(OMC_LAGCN_AUC)),
                 sd(as.numeric(HGIMC_LAGCN_AUC)),
                 sd(IBCF_LAGCN_AUC$V2), sd(LIBMF_LAGCN_AUC$V2))

median_AUPR_LAGCN = c(median(as.numeric(NMF_LAGCN_AUPR$V2)),
                      median(as.numeric(NMFPer_LAGCN_AUPR$V1)),
                      median(as.numeric(NMFDR_LAGCN_AUPR)),
                      median(as.numeric(VDA_LAGCN_AUPR)),
                      median(as.numeric(BNNR_LAGCN_AUPR)),
                      median(as.numeric(OMC_LAGCN_AUPR)),
                      median(as.numeric(HGIMC_LAGCN_AUPR)),
                      median(IBCF_LAGCN_AUPR$V2), median(LIBMF_LAGCN_AUPR$V2))


sd_AUPR_LAGCN = c(sd(as.numeric(NMF_LAGCN_AUPR$V2)),
                  sd(as.numeric(NMFPer_LAGCN_AUPR$V1)),
                  sd(as.numeric(NMFDR_LAGCN_AUPR)),
                  sd(as.numeric(VDA_LAGCN_AUPR)),
                  sd(as.numeric(BNNR_LAGCN_AUPR)),
                  sd(as.numeric(OMC_LAGCN_AUPR)),
                  sd(as.numeric(HGIMC_LAGCN_AUPR)),
                  sd(IBCF_LAGCN_AUPR$V2), sd(LIBMF_LAGCN_AUPR$V2))



## Fdata
median_AUC_Fdata = c(median(as.numeric(NMF_Fdata_AUC$V2)),median(as.numeric(NMFPer_Fdata_AUC)), median(as.numeric(NMFDR_Fdata_AUC)),
                     median(as.numeric(VDA_Fdata_AUC)),median(as.numeric(BNNR_Fdata_AUC)),
                     median(as.numeric(OMC_Fdata_AUC)),median(as.numeric(HGIMC_Fdata_AUC)), median(IBCF_Fdata_AUC$V2), median(LIBMF_Fdata_AUC$V2))

sd_AUC_Fdata = c(sd(as.numeric(NMF_Fdata_AUC$V2)),sd(as.numeric(NMFPer_Fdata_AUC)), sd(as.numeric(NMFDR_Fdata_AUC)),
                 sd(as.numeric(VDA_Fdata_AUC)),sd(as.numeric(BNNR_Fdata_AUC)),
                 sd(as.numeric(OMC_Fdata_AUC)),sd(as.numeric(HGIMC_Fdata_AUC)), sd(IBCF_Fdata_AUC$V2), sd(LIBMF_Fdata_AUC$V2))


median_AUPR_Fdata = c(median(as.numeric(NMF_Fdata_AUPR$V2)),median(as.numeric(NMFPer_Fdata_AUPR)), median(as.numeric(NMFDR_Fdata_AUPR)),
                      median(as.numeric(VDA_Fdata_AUPR)),median(as.numeric(BNNR_Fdata_AUPR)),
                      median(as.numeric(OMC_Fdata_AUPR)),median(as.numeric(HGIMC_Fdata_AUPR)), median(IBCF_Fdata_AUPR$V2), median(LIBMF_Fdata_AUPR$V2))

sd_AUPR_Fdata = c(sd(as.numeric(NMF_Fdata_AUPR$V2)),sd(as.numeric(NMFPer_Fdata_AUPR)), sd(as.numeric(NMFDR_Fdata_AUPR)),
                  sd(as.numeric(VDA_Fdata_AUPR)),sd(as.numeric(BNNR_Fdata_AUPR)),
                  sd(as.numeric(OMC_Fdata_AUPR)),sd(as.numeric(HGIMC_Fdata_AUPR)), sd(IBCF_Fdata_AUPR$V2), sd(LIBMF_Fdata_AUPR$V2))




## 4. Cdataset
median_AUC_Cdata =c(median(as.numeric(NMF_Cdata_AUC$V2)),
                    median(as.numeric(NMFPer_Cdata_AUC)),median(as.numeric(NMFDR_Cdata_AUC)),
                    median(as.numeric(VDA_Cdata_AUC)),
                    median(as.numeric(BNNR_Cdata_AUC)),
                    median(as.numeric(OMC_Cdata_AUC)),median(as.numeric(HGIMC_Cdata_AUC)), median(IBCF_Cdata_AUC$V2), 
                    median(LIBMF_Cdata_AUC$V2))

sd_AUC_Cdata =c(sd(as.numeric(NMF_Cdata_AUC$V2)),
                sd(as.numeric(NMFPer_Cdata_AUC)),sd(as.numeric(NMFDR_Cdata_AUC)),
                sd(as.numeric(VDA_Cdata_AUC)),
                sd(as.numeric(BNNR_Cdata_AUC)),
                sd(as.numeric(OMC_Cdata_AUC)),sd(as.numeric(HGIMC_Cdata_AUC)), sd(IBCF_Cdata_AUC$V2), sd(LIBMF_Cdata_AUC$V2))

median_AUPR_Cdata =c(median(as.numeric(NMF_Cdata_AUPR$V2)),
                     median(as.numeric(NMFPer_Cdata_AUPR)),median(as.numeric(NMFDR_Cdata_AUPR)),
                     median(as.numeric(VDA_Cdata_AUPR)),
                     median(as.numeric(BNNR_Cdata_AUPR)),
                     median(as.numeric(OMC_Cdata_AUPR)),median(as.numeric(HGIMC_Cdata_AUPR)), median(IBCF_Cdata_AUPR$V2), 
                     median(LIBMF_Cdata_AUPR$V2))

sd_AUPR_Cdata =c(sd(as.numeric(NMF_Cdata_AUPR$V2)),
                 sd(as.numeric(NMFPer_Cdata_AUPR)),sd(as.numeric(NMFDR_Cdata_AUPR)),
                 sd(as.numeric(VDA_Cdata_AUPR)),
                 sd(as.numeric(BNNR_Cdata_AUPR)),
                 sd(as.numeric(OMC_Cdata_AUPR)),sd(as.numeric(HGIMC_Cdata_AUPR)), sd(IBCF_Cdata_AUPR$V2), sd(LIBMF_Cdata_AUPR$V2))



## 5. LRSSL
median_AUC_LRSSL = c(median(as.numeric(NMF_LRSSL_AUC$V2)),
                     median(as.numeric(NMFPer_LRSSL_AUC)),
                     median(as.numeric(NMFDR_LRSSL_AUC)),
                     median(as.numeric(VDA_LRSSL_AUC)),
                     median(as.numeric(BNNR_LRSSL_AUC)),
                     median(as.numeric(OMC_LRSSL_AUC)),
                     median(as.numeric(HGIMC_LRSSL_AUC)),
                     median(IBCF_LRSSL_AUC$V2),
                     median(LIBMF_LRSSL_AUC$V2))

sd_AUC_LRSSL = c(sd(as.numeric(NMF_LRSSL_AUC$V2)),
                 sd(as.numeric(NMFPer_LRSSL_AUC)),
                 sd(as.numeric(NMFDR_LRSSL_AUC)),
                 sd(as.numeric(VDA_LRSSL_AUC)),
                 sd(as.numeric(BNNR_LRSSL_AUC)),
                 sd(as.numeric(OMC_LRSSL_AUC)),
                 sd(as.numeric(HGIMC_LRSSL_AUC)),
                 sd(IBCF_LRSSL_AUC$V2),
                 sd(LIBMF_LRSSL_AUC$V2))

median_AUPR_LRSSL = c(median(as.numeric(NMF_LRSSL_AUPR$V2)),
                      median(as.numeric(NMFPer_LRSSL_AUPR)),
                      median(as.numeric(NMFDR_LRSSL_AUPR)),
                      median(as.numeric(VDA_LRSSL_AUPR)),
                      median(as.numeric(BNNR_LRSSL_AUPR)),
                      median(as.numeric(OMC_LRSSL_AUPR)),
                      median(as.numeric(HGIMC_LRSSL_AUPR)),
                      median(IBCF_LRSSL_AUPR$V2),
                      median(LIBMF_LRSSL_AUPR$V2))

sd_AUPR_LRSSL = c(sd(as.numeric(NMF_LRSSL_AUPR$V2)),
                  sd(as.numeric(NMFPer_LRSSL_AUPR)),
                  sd(as.numeric(NMFDR_LRSSL_AUPR)),
                  sd(as.numeric(VDA_LRSSL_AUPR)),
                  sd(as.numeric(BNNR_LRSSL_AUPR)),
                  sd(as.numeric(OMC_LRSSL_AUPR)),
                  sd(as.numeric(HGIMC_LRSSL_AUPR)),
                  sd(IBCF_LRSSL_AUPR$V2),
                  sd(LIBMF_LRSSL_AUPR$V2))



## 6. Ydataset 
median_AUC_Ydata = c(median(as.numeric(NMF_Ydata_AUC$V2)),
                     median(as.numeric(NMFPer_Ydata_AUC)),
                     median(as.numeric(NMFDR_Ydata_AUC)),
                     median(as.numeric(VDA_Ydata_AUC)),
                     median(as.numeric(BNNR_Ydata_AUC)),
                     median(as.numeric(OMC_Ydata_AUC)),
                     median(as.numeric(HGIMC_Ydata_AUC)),
                     median(IBCF_Ydata_AUC$V2),median(LIBMF_Ydata_AUC$V2))
sd_AUC_Ydata = c(sd(as.numeric(NMF_Ydata_AUC$V2)),
                 sd(as.numeric(NMFPer_Ydata_AUC)),
                 sd(as.numeric(NMFDR_Ydata_AUC)),
                 sd(as.numeric(VDA_Ydata_AUC)),
                 sd(as.numeric(BNNR_Ydata_AUC)),
                 sd(as.numeric(OMC_Ydata_AUC)),
                 sd(as.numeric(HGIMC_Ydata_AUC)),
                 sd(IBCF_Ydata_AUC$V2),sd(LIBMF_Ydata_AUC$V2))

median_AUPR_Ydata = c(median(as.numeric(NMF_Ydata_AUPR$V2)),
                      median(as.numeric(NMFPer_Ydata_AUPR)),
                      median(as.numeric(NMFDR_Ydata_AUPR)),
                      median(as.numeric(VDA_Ydata_AUPR)),
                      median(as.numeric(BNNR_Ydata_AUPR)),
                      median(as.numeric(OMC_Ydata_AUPR)),
                      median(as.numeric(HGIMC_Ydata_AUPR)),
                      median(IBCF_Ydata_AUPR$V2),median(LIBMF_Ydata_AUPR$V2))
sd_AUPR_Ydata = c(sd(as.numeric(NMF_Ydata_AUPR$V2)),
                  sd(as.numeric(NMFPer_Ydata_AUPR)),
                  sd(as.numeric(NMFDR_Ydata_AUPR)),
                  sd(as.numeric(VDA_Ydata_AUPR)),
                  sd(as.numeric(BNNR_Ydata_AUPR)),
                  sd(as.numeric(OMC_Ydata_AUPR)),
                  sd(as.numeric(HGIMC_Ydata_AUPR)),
                  sd(IBCF_Ydata_AUPR$V2),sd(LIBMF_Ydata_AUPR$V2))

## 11. oMat data
median_AUC_oMat = c(median(as.numeric(NMF_oMat_AUC$V2)),median(as.numeric(NMFPer_oMat_AUC)),
                    median(as.numeric(NMFDR_oMat_AUC)),median(as.numeric(VDA_oMat_AUC)),
                    median(as.numeric(BNNR_oMat_AUC)),
                    median(as.numeric(OMC_oMat_AUC)),median(as.numeric(HGIMC_oMat_AUC)),
                    median(IBCF_oMat_AUC$V2), median(LIBMF_oMat_AUC$V2))


sd_AUC_oMat = c(sd(as.numeric(NMF_oMat_AUC$V2)),sd(as.numeric(NMFPer_oMat_AUC)),
                sd(as.numeric(NMFDR_oMat_AUC)),sd(as.numeric(VDA_oMat_AUC)),
                sd(as.numeric(BNNR_oMat_AUC)),
                sd(as.numeric(OMC_oMat_AUC)),sd(as.numeric(HGIMC_oMat_AUC)),
                sd(IBCF_oMat_AUC$V2), sd(LIBMF_oMat_AUC$V2))


median_AUPR_oMat = c(median(as.numeric(NMF_oMat_AUPR$V2)),median(as.numeric(NMFPer_oMat_AUPR)),
                     median(as.numeric(NMFDR_oMat_AUPR)),median(as.numeric(VDA_oMat_AUPR)),
                     median(as.numeric(BNNR_oMat_AUPR)),
                     median(as.numeric(OMC_oMat_AUPR)),median(as.numeric(HGIMC_oMat_AUPR)),
                     median(IBCF_oMat_AUPR$V2), median(LIBMF_oMat_AUPR$V2))


sd_AUPR_oMat = c(sd(as.numeric(NMF_oMat_AUPR$V2)),sd(as.numeric(NMFPer_oMat_AUPR)),
                 sd(as.numeric(NMFDR_oMat_AUPR)),sd(as.numeric(VDA_oMat_AUPR)),
                 sd(as.numeric(BNNR_oMat_AUPR)),
                 sd(as.numeric(OMC_oMat_AUPR)),sd(as.numeric(HGIMC_oMat_AUPR)),
                 sd(IBCF_oMat_AUPR$V2), sd(LIBMF_oMat_AUPR$V2))



## 12. hsdn
median_AUC_hsdn = c(median(as.numeric(NMF_hsdn_AUC$V2)),median(as.numeric(NMFPer_hsdn_AUC)), 
                    median(as.numeric(NMFDR_hsdn_AUC)),median(as.numeric(VDA_hsdn_AUC)),
                    median(as.numeric(BNNR_hsdn_AUC)),
                    median(as.numeric(OMC_hsdn_AUC)),median(as.numeric(HGIMC_hsdn_AUC)),
                    median(IBCF_hsdn_AUC$V2), median(LIBMF_hsdn_AUC$V2))


sd_AUC_hsdn = c(sd(as.numeric(NMF_hsdn_AUC$V2)),sd(as.numeric(NMFPer_hsdn_AUC)), 
                sd(as.numeric(NMFDR_hsdn_AUC)),sd(as.numeric(VDA_hsdn_AUC)),
                sd(as.numeric(BNNR_hsdn_AUC)),
                sd(as.numeric(OMC_hsdn_AUC)),sd(as.numeric(HGIMC_hsdn_AUC)),
                sd(IBCF_hsdn_AUC$V2), sd(LIBMF_hsdn_AUC$V2))

median_AUPR_hsdn = c(median(as.numeric(NMF_hsdn_AUPR$V2)),median(as.numeric(NMFPer_hsdn_AUPR)), 
                     median(as.numeric(NMFDR_hsdn_AUPR)),median(as.numeric(VDA_hsdn_AUPR)),
                     median(as.numeric(BNNR_hsdn_AUPR)),
                     median(as.numeric(OMC_hsdn_AUPR)),median(as.numeric(HGIMC_hsdn_AUPR)),
                     median(IBCF_hsdn_AUPR$V2), median(LIBMF_hsdn_AUPR$V2))


sd_AUPR_hsdn = c(sd(as.numeric(NMF_hsdn_AUPR$V2)),sd(as.numeric(NMFPer_hsdn_AUPR)), 
                 sd(as.numeric(NMFDR_hsdn_AUPR)),sd(as.numeric(VDA_hsdn_AUPR)),
                 sd(as.numeric(BNNR_hsdn_AUPR)),
                 sd(as.numeric(OMC_hsdn_AUPR)),sd(as.numeric(HGIMC_hsdn_AUPR)),
                 sd(IBCF_hsdn_AUPR$V2), sd(LIBMF_hsdn_AUPR$V2))





# Install and load pheatmap if not already installed
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
library(pheatmap)
AUC_median_matrix <- cbind(median_AUC_HDVD,
                           median_AUC_LAGCN,
                           median_AUC_Fdata,
                           median_AUC_Cdata,
                           median_AUC_LRSSL,
                           median_AUC_Ydata,
                           median_AUC_oMat,
                           median_AUC_hsdn)


AUC_sd_matrix <- cbind(sd_AUC_HDVD,
                       sd_AUC_LAGCN,
                       sd_AUC_Fdata,
                       sd_AUC_Cdata,
                       sd_AUC_LRSSL,
                       sd_AUC_Ydata,
                       sd_AUC_oMat,
                       sd_AUC_hsdn)

AUPR_median_matrix <- cbind(median_AUPR_HDVD,
                            median_AUPR_LAGCN,
                            median_AUPR_Fdata,
                            median_AUPR_Cdata,
                            median_AUPR_LRSSL,
                            median_AUPR_Ydata,
                            median_AUPR_oMat,
                            median_AUPR_hsdn)


AUPR_sd_matrix <- cbind(sd_AUPR_HDVD,
                        sd_AUPR_LAGCN,
                        sd_AUPR_Fdata,
                        sd_AUPR_Cdata,
                        sd_AUPR_LRSSL,
                        sd_AUPR_Ydata,
                        sd_AUPR_oMat,
                        sd_AUPR_hsdn)


# save(AUC_median_matrix, file = "AUC_median_matrix_rerun_neworder.RData")
# save(AUC_sd_matrix, file = "AUC_sd_matrix_rerun_neworder.RData")
# save(AUPR_median_matrix, file = "AUPR_median_matrix_rerun_neworder.RData")
# save(AUPR_sd_matrix, file = "AUPR_sd_matrix_rerun_neworder.RData")


# load("AUC_median_matrix_rerun_neworder.RData")
# load("AUC_sd_matrix_rerun_neworder.RData")
# load("AUPR_median_matrix_rerun_neworder.RData")
# load("AUPR_sd_matrix_rerun_neworder.RData")
### Plot AUC and AUPR in the heatmap


if (!require("ggplot2")) install.packages("ggplot2")
if (!require("reshape2")) install.packages("reshape2")
if (!require("scales")) install.packages("scales")
if (!require("dplyr")) install.packages("dplyr") 


library(ggplot2)
library(reshape2)
library(scales)
library(dplyr)


# Define row and column names
colnames(AUC_median_matrix) <- c("HDVD", "LAGCN", "Fdata", "Cdata", "LRSSL", "Ydata", "oMat-MechDB", "HSDN-MechDB")
rownames(AUC_median_matrix) <- c("NMF", "NMF-PDR", "NMF-DR", "VDA-GKSBMF", "BNNR", "OMC", "HGIMC", "IBCF", "LIBMF")

colnames(AUPR_median_matrix) <- c("HDVD", "LAGCN", "Fdata", "Cdata", "LRSSL", "Ydata", "oMat-MechDB", "HSDN-MechDB")
rownames(AUPR_median_matrix) <- c("NMF", "NMF-PDR", "NMF-DR", "VDA-GKSBMF", "BNNR", "OMC", "HGIMC", "IBCF", "LIBMF")


# Create a data frame from the AUC matrix
df <- as.data.frame(as.table(AUC_median_matrix)) %>%
  rename(method = Var1, dataset = Var2, median = Freq)

# Add standard deviation (AUPR) values to the data frame
df$AUPR <- as.vector(AUPR_median_matrix)

# Inspect the data frame
head(df)

# Plot the bubble heatmap with smaller circle sizes
ggplot(df, aes(x = dataset, y = method)) +
  geom_point(aes(size = AUPR, color = median)) +  # Use 'median' for color and 'sd' for size
  scale_color_gradient2(low = "red", mid = "white", high = "darkblue", midpoint = 0.5) +
  
  
  
  scale_size(
    range = c(4, 12),  # Adjust the overall range of circle sizes
    breaks = seq(0.05, 0.3, by = 0.05),  # Specify custom increments
    guide = guide_legend(title = "AUPR (Size)",
                         override.aes = list(shape = 1, color = "black", stroke = 1.5))
    
  ) +
  theme_minimal() +
  theme(
    
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
    
    #plot.title = element_text(size = 16, hjust = 0.5),  
    
    axis.title.x = element_text(size = 16, color = "black",hjust = 0.5),  # x-axis label font, face = "bold",
    axis.title.y = element_text(size = 16, color = "black",hjust = 0.5),  # y-axis label font
    
    
    axis.text.x = element_text(size = 14,  angle = 45, hjust = 0.8, color = "black", face ="plain"),  # x-axis ticks font
    axis.text.y = element_text(size = 14, angle = 0, hjust = 0.8, color = "black",face ="plain"),   # y-axis ticks font
    
    legend.title = element_text(size = 14),  # Legend title size
    legend.text = element_text(size = 14)    # Legend text size
    # legend.position = "none"
  )+
  labs(
    # title = "Heatmap of Median AUC (Color) and AUPR (Size)",
    x = "Datasets",
    y = "Methods",
    color = "AUC (Color)"
  ) +
  coord_fixed()


### plot sd values

# Define row and column names
colnames(AUC_sd_matrix) <- c("HDVD", "LAGCN", "Fdata", "Cdata", "LRSSL", "Ydata", "oMat-MechDB", "HSDN-MechDB")
rownames(AUC_sd_matrix) <- c("NMF", "NMF-PDR", "NMF-DR", "VDA-GKSBMF", "BNNR", "OMC", "HGIMC", "IBCF", "LIBMF")

colnames(AUPR_sd_matrix) <- c("HDVD", "LAGCN", "Fdata", "Cdata", "LRSSL", "Ydata", "oMat-MechDB", "HSDN-MechDB")
rownames(AUPR_sd_matrix) <- c("NMF", "NMF-PDR", "NMF-DR", "VDA-GKSBMF", "BNNR", "OMC", "HGIMC", "IBCF", "LIBMF")


# Create a data frame from the AUC matrix
df <- as.data.frame(as.table(AUC_sd_matrix)) %>%
  rename(method = Var1, dataset = Var2, sd_AUC = Freq)

# Add standard deviation (AUPR) values to the data frame
df$sd_AUPR <- as.vector(AUPR_sd_matrix)

# Inspect the data frame
head(df)
x0 = min(AUC_sd_matrix, AUPR_sd_matrix) # = 0.0002413462
x1 = max(AUC_sd_matrix, AUPR_sd_matrix) # 0.01242173

ggplot(df, aes(x = dataset, y = method)) +
  geom_point(aes(size = sd_AUPR, color = sd_AUC)) +  # Use 'median' for color and 'sd' for size
  # scale_color_gradient2(low = "red", mid = "white", high = "darkblue", midpoint = 0.5) +
  scale_fill_gradient(
    low = "red", high = "darkblue", name = "SD (AUC)", 
    limits = c(min(df$sd_AUC, na.rm = TRUE), max(df$sd_AUC, na.rm = TRUE)))+
  scale_size(
    range = c(1, 8),  # Adjust the overall range of circle sizes
    breaks = seq(x0, x1, length.out = 5),  # Specify custom increments
    labels = scales::scientific_format(digits = 3),
    # guide = guide_legend(title = "SD (AUPR) (Size)", color = NA),
    # override.aes = list(shape = 1, color = "black", stroke = 1.5)
    guide = guide_legend(
      title = "SD (AUPR) (Size)",
      override.aes = list(shape = 1, color = "black", stroke = 1.5)
    )
  ) +
  theme_minimal() +
  theme(
    
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
    
    plot.title = element_text(size = 16, hjust = 0.5),  
    
    axis.title.x = element_text(size = 16, color = "black",hjust = 0.5),  # x-axis label font, face = "bold",
    axis.title.y = element_text(size = 16, color = "black",hjust = 0.5),  # y-axis label font
    
    
    axis.text.x = element_text(size = 14,  angle = 45, hjust = 0.8, color = "black", face ="plain"),  # x-axis ticks font
    axis.text.y = element_text(size = 14, angle = 0, hjust = 0.8, color = "black",face ="plain"),   # y-axis ticks font
    
    legend.title = element_text(size = 14),  # Legend title size
    legend.text = element_text(size = 14)    # Legend text size
    # legend.position = "none"
    
  ) +
  labs(
    #title = "Heatmap of SD values of median of AUC and AUPR",
    x = "Datasets",
    y = "Methods",
    color = "SD (AUC) (Color)"
  ) +
  coord_fixed()

#####

median_matrix = AUC_median_matrix
sd_matrix = AUC_sd_matrix

colnames(median_matrix) <- c("HDVD", "LAGCN","Fdata", "Cdata", "LRSSL", "Ydata", "oMat-MechDB", "hsdn-MechDB")
rownames(median_matrix) <- c("NMF", "NMFPer","NMFDR", "VDA-GKSBMF", "BNNR", "OMC", "HGIMC", "IBCF", "LIBMF")

colnames(sd_matrix) <- c("HDVD", "LAGCN","Fdata", "Cdata", "LRSSL", "Ydata", "oMat-MechDB", "hsdn-MechDB")
rownames(sd_matrix) <- c("NMF", "NMFPer","NMFDR", "VDA-GKSBMF", "BNNR", "OMC", "HGIMC", "IBCF", "LIBMF")




### boxplot of AUC values

par(mfrow = c(2,4))
boxplot(list(as.numeric(NMF_HDVD_AUC$V2),as.numeric(NMFPer_HDVD_AUC$V2),as.numeric(NMFDR_HDVD_AUC),as.numeric(VDA_HDVD_AUC), as.numeric(BNNR_HDVD_AUC),
             as.numeric(OMC_HDVD_AUC),as.numeric(HGIMC_HDVD_AUC),IBCF_HDVD_AUC$V2, LIBMF_HDVD_AUC$V2),
        main = "HDVD", ylab = "AUC Values", 
        cex.lab = 1.2,
        cex.axis = 1,
        font.axis = 2,
        boxwex = 1,
        las = 2,
        names = c("NMF", "NMFPer","NMFDR", "VDA-GKSBMF", "BNNR", "OMC", "HGIMC", "IBCF", "LIBMF"),
        col = c("lightblue", "lightgreen", "lightpink", "lightgray", "lightyellow", "gold", "deeppink", "blue", "coral"))

boxplot(list(as.numeric(NMF_LAGCN_AUC$V2),as.numeric(NMFPer_LAGCN_AUC$V1), as.numeric(NMFDR_LAGCN_AUC), as.numeric(VDA_LAGCN_AUC),as.numeric(BNNR_LAGCN_AUC),
             as.numeric(OMC_LAGCN_AUC),as.numeric(HGIMC_LAGCN_AUC), IBCF_LAGCN_AUC$V2,LIBMF_LAGCN_AUC$V2 ),
        main = "LAGCN", ylab = "AUC Values", 
        cex.lab = 1.2,
        cex.axis = 1,
        font.axis = 2,
        boxwex = 1,
        las = 2,
        names = c("NMF", "NMFPer","NMFDR", "VDA-GKSBMF", "BNNR", "OMC", "HGIMC", "IBCF", "LIBMF"),
        col = c("lightblue", "lightgreen", "lightpink", "lightgray", "lightyellow", "gold", "deeppink",  "blue", "coral"))

boxplot(list(as.numeric(NMF_Fdata_AUC$V2),as.numeric(NMFPer_Fdata_AUC), as.numeric(NMFDR_Fdata_AUC),as.numeric(VDA_Fdata_AUC),as.numeric(BNNR_Fdata_AUC),
             as.numeric(OMC_Fdata_AUC),as.numeric(HGIMC_Fdata_AUC), IBCF_Fdata_AUC$V2, LIBMF_Fdata_AUC$V2),
        main = "Fdata", ylab = "AUC Values", 
        cex.lab = 1.2,
        cex.axis = 1,
        font.axis = 2,
        boxwex = 1,
        las = 2,
        names = c("NMF", "NMFPer","NMFDR", "VDA-GKSBMF", "BNNR", "OMC", "HGIMC", "IBCF", "LIBMF"),
        col = c("lightblue", "lightgreen", "lightpink", "lightgray", "lightyellow", "gold", "deeppink",  "blue", "coral"))

boxplot(list(as.numeric(NMF_Cdata_AUC$V2),as.numeric(NMFPer_Cdata_AUC),as.numeric(NMFDR_Cdata_AUC), as.numeric(VDA_Cdata_AUC),
             as.numeric(BNNR_Cdata_AUC),
             as.numeric(OMC_Cdata_AUC),as.numeric(HGIMC_Cdata_AUC), IBCF_Cdata_AUC$V2, LIBMF_Cdata_AUC$V2),
        main = "Cdata", ylab = "AUC Values", 
        cex.lab = 1.2,
        cex.axis = 1,
        font.axis = 2,
        boxwex = 1,
        las = 2,
        names = c("NMF", "NMFPer","NMFDR", "VDA-GKSBMF", "BNNR", "OMC", "HGIMC", "IBCF", "LIBMF"),
        col = c("lightblue", "lightgreen", "lightpink", "lightgray", "lightyellow", "gold", "deeppink",  "blue", "coral"))


boxplot(list(as.numeric(NMF_LRSSL_AUC$V2),as.numeric(NMFPer_LRSSL_AUC),as.numeric(NMFDR_LRSSL_AUC),
             as.numeric(VDA_LRSSL_AUC), as.numeric(BNNR_LRSSL_AUC),
             as.numeric(OMC_LRSSL_AUC),as.numeric(HGIMC_LRSSL_AUC), IBCF_LRSSL_AUC$V2,LIBMF_LRSSL_AUC$V2 ),
        main = "LRSSL", ylab = "AUC Values", 
        cex.lab = 1.2,
        cex.axis = 1,
        font.axis = 2,
        boxwex = 1,
        las = 2,
        names = c("NMF", "NMFPer","NMFDR", "VDA-GKSBMF", "BNNR", "OMC", "HGIMC", "IBCF", "LIBMF"),
        col = c("lightblue", "lightgreen", "lightpink", "lightgray", "lightyellow", "gold", "deeppink",  "blue", "coral"))


boxplot(list(as.numeric(NMF_Ydata_AUC$V2),as.numeric(NMFPer_Ydata_AUC),as.numeric(NMFDR_Ydata_AUC),
             as.numeric(VDA_Ydata_AUC), as.numeric(BNNR_Ydata_AUC),
             as.numeric(OMC_Ydata_AUC),as.numeric(HGIMC_Ydata_AUC), IBCF_Ydata_AUC$V2,LIBMF_Ydata_AUC$V2 ),
        main = "Ydata", ylab = "AUC Values", 
        cex.lab = 1.2,
        cex.axis = 1,
        font.axis = 2,
        boxwex = 1,
        las = 2,
        names = c("NMF", "NMFPer","NMFDR", "VDA-GKSBMF", "BNNR", "OMC", "HGIMC", "IBCF", "LIBMF"),
        col = c("lightblue", "lightgreen", "lightpink", "lightgray", "lightyellow", "gold", "deeppink",  "blue", "coral"))


boxplot(list(as.numeric(NMF_oMat_AUC$V2),as.numeric(NMFPer_oMat_AUC),as.numeric(NMFDR_oMat_AUC),as.numeric(VDA_oMat_AUC), as.numeric(BNNR_oMat_AUC),
             as.numeric(OMC_oMat_AUC),as.numeric(HGIMC_oMat_AUC), IBCF_oMat_AUC$V2, LIBMF_oMat_AUC$V2),
        main = "OMat-MechDB", ylab = "AUC Values", 
        cex.lab = 1.2,
        cex.axis = 1,
        font.axis = 2,
        boxwex = 1,
        las = 2,
        names = c("NMF", "NMFPer","NMFDR", "VDA-GKSBMF", "BNNR", "OMC", "HGIMC", "IBCF", "LIBMF"),
        col = c("lightblue", "lightgreen", "lightpink", "lightgray", "lightyellow", "gold", "deeppink",  "blue", "coral"))

boxplot(list(as.numeric(NMF_hsdn_AUC$V2),as.numeric(NMFPer_hsdn_AUC), as.numeric(NMFDR_hsdn_AUC),as.numeric(VDA_hsdn_AUC), 
             as.numeric(BNNR_hsdn_AUC),
             as.numeric(OMC_hsdn_AUC),as.numeric(HGIMC_hsdn_AUC), 
             IBCF_hsdn_AUC$V2, LIBMF_hsdn_AUC$V2),
        main = "hsdn-MechDB", ylab = "AUC Values", 
        cex.lab = 1.2,
        cex.axis = 1,
        font.axis = 2,
        boxwex = 1,
        las = 2,
        names = c("NMF", "NMFPer","NMFDR", "VDA-GKSBMF", "BNNR", "OMC", "HGIMC", "IBCF", "LIBMF"),
        col = c("lightblue", "lightgreen", "lightpink", "lightgray", "lightyellow", "gold", "deeppink",  "blue", "coral"))


# dev.off()
par(mfrow = c(1,1))


#### Plot boxplot AUPR
par(mfrow = c(2, 4))
boxplot(list(as.numeric(NMF_HDVD_AUPR$V2),as.numeric(NMFPer_HDVD_AUPR$V1), as.numeric(NMFDR_HDVD_AUPR), 
             as.numeric(VDA_HDVD_AUPR), as.numeric(BNNR_HDVD_AUPR$V1),as.numeric(OMC_HDVD_AUPR),
             as.numeric(HGIMC_HDVD_AUPR), IBCF_HDVD_AUPR$V2, LIBMF_HDVD_AUPR$V2),
        main = "HDVD", ylab = "AUPR Values", 
        cex.lab = 1.2,
        cex.axis = 1,
        font.axis = 2,
        boxwex = 1,
        las = 2,
        names = c("NMF", "NMFPer","NMFDR", "VDA-GKSBMF", "BNNR", "OMC", "HGIMC", "IBCF", "LIBMF"),
        col = c("lightblue", "lightgreen", "lightpink", "lightgray", "lightyellow", "gold", "deeppink",  "blue", "coral"))


boxplot(list(as.numeric(NMF_LAGCN_AUPR$V2),as.numeric(NMFPer_LAGCN_AUPR$V1),as.numeric(NMFDR_LAGCN_AUPR),as.numeric(VDA_LAGCN_AUPR), as.numeric(BNNR_LAGCN_AUPR),
             as.numeric(OMC_LAGCN_AUPR),as.numeric(HGIMC_LAGCN_AUPR), IBCF_LAGCN_AUPR$V2,LIBMF_LAGCN_AUPR$V2),
        main = "LAGCN", ylab = "AUPR Values", 
        cex.lab = 1.2,
        cex.axis = 1,
        font.axis = 2,
        boxwex = 1,
        las = 2,
        names = c("NMF", "NMFPer","NMFDR", "VDA-GKSBMF", "BNNR", "OMC", "HGIMC", "IBCF", "LIBMF"),
        col = c("lightblue", "lightgreen", "lightpink", "lightgray", "lightyellow", "gold", "deeppink",  "blue", "coral"))

boxplot(list(as.numeric(NMF_Fdata_AUPR$V2),as.numeric(NMFPer_Fdata_AUPR), as.numeric(NMFDR_Fdata_AUPR), as.numeric(VDA_Fdata_AUPR),
             as.numeric(BNNR_Fdata_AUPR),as.numeric(OMC_Fdata_AUPR),
             as.numeric(HGIMC_Fdata_AUPR), IBCF_Fdata_AUPR$V2, LIBMF_Fdata_AUPR$V2),
        main = "Fdata", ylab = "AUPR Values", 
        cex.lab = 1.2,
        cex.axis = 1,
        font.axis = 2,
        boxwex = 1,
        las = 2,
        names = c("NMF", "NMFPer","NMFDR", "VDA-GKSBMF", "BNNR", "OMC", "HGIMC", "IBCF", "LIBMF"),
        col = c("lightblue", "lightgreen", "lightpink", "lightgray", "lightyellow", "gold", "deeppink",  "blue", "coral"))


boxplot(list(as.numeric(NMF_Cdata_AUPR$V2),as.numeric(NMFPer_Cdata_AUPR),as.numeric(NMFDR_Cdata_AUPR), as.numeric(VDA_Cdata_AUPR),as.numeric(BNNR_Cdata_AUPR),
             as.numeric(OMC_Cdata_AUPR),as.numeric(HGIMC_Cdata_AUPR), IBCF_Cdata_AUPR$V2, LIBMF_Cdata_AUPR$V2),
        main = "Cdata", ylab = "AUPR Values", 
        cex.lab = 1.2,
        cex.axis = 1,
        font.axis = 2,
        boxwex = 1,
        las = 2,
        names = c("NMF", "NMFPer","NMFDR", "VDA-GKSBMF", "BNNR", "OMC", "HGIMC", "IBCF", "LIBMF"),
        col = c("lightblue", "lightgreen", "lightpink", "lightgray", "lightyellow", "gold", "deeppink",  "blue", "coral"))


boxplot(list(as.numeric(NMF_LRSSL_AUPR$V2),as.numeric(NMFPer_LRSSL_AUPR),as.numeric(NMFDR_LRSSL_AUPR),
             as.numeric(VDA_LRSSL_AUPR), as.numeric(BNNR_LRSSL_AUPR),
             as.numeric(OMC_LRSSL_AUPR),as.numeric(HGIMC_LRSSL_AUPR), IBCF_LRSSL_AUPR$V2, LIBMF_LRSSL_AUPR$V2),
        main = "LRSSL", ylab = "AUPR Values", 
        cex.lab = 1.2,
        cex.axis = 1,
        font.axis = 2,
        boxwex = 1,
        las = 2,
        names = c("NMF", "NMFPer","NMFDR", "VDA-GKSBMF", "BNNR", "OMC", "HGIMC", "IBCF", "LIBMF"),
        col = c("lightblue", "lightgreen", "lightpink", "lightgray", "lightyellow", "gold", "deeppink",  "blue", "coral"))


boxplot(list(as.numeric(NMF_Ydata_AUPR$V2),as.numeric(NMFPer_Ydata_AUPR),as.numeric(NMFDR_Ydata_AUPR),
             as.numeric(VDA_Ydata_AUPR), as.numeric(BNNR_Ydata_AUPR),
             as.numeric(OMC_Ydata_AUPR),as.numeric(HGIMC_Ydata_AUPR), IBCF_Ydata_AUPR$V2, LIBMF_Ydata_AUPR$V2),
        main = "Ydata", ylab = "AUPR Values", 
        cex.lab = 1.2,
        cex.axis = 1,
        font.axis = 2,
        boxwex = 1,
        las = 2,
        names = c("NMF", "NMFPer","NMFDR", "VDA-GKSBMF", "BNNR", "OMC", "HGIMC", "IBCF", "LIBMF"),
        col = c("lightblue", "lightgreen", "lightpink", "lightgray", "lightyellow", "gold", "deeppink",  "blue", "coral"))


boxplot(list(as.numeric(NMF_oMat_AUPR$V2),as.numeric(NMFPer_oMat_AUPR), as.numeric(NMFDR_oMat_AUPR),as.numeric(VDA_oMat_AUPR),as.numeric(BNNR_oMat_AUPR),
             as.numeric(OMC_oMat_AUPR),as.numeric(HGIMC_oMat_AUPR), IBCF_oMat_AUPR$V2, LIBMF_oMat_AUPR$V2),
        main = "oMat-MechDB", ylab = "AUPR Values", 
        cex.lab = 1.2,
        cex.axis = 1,
        font.axis = 2,
        boxwex = 1,
        las = 2,
        names = c("NMF", "NMFPer","NMFDR", "VDA-GKSBMF", "BNNR", "OMC", "HGIMC", "IBCF", "LIBMF"),
        col = c("lightblue", "lightgreen", "lightpink", "lightgray", "lightyellow", "gold", "deeppink",  "blue", "coral"))

boxplot(list(as.numeric(NMF_hsdn_AUPR$V2),as.numeric(NMFPer_hsdn_AUPR),as.numeric(NMFDR_hsdn_AUPR),as.numeric(VDA_hsdn_AUPR), as.numeric(BNNR_hsdn_AUPR),
             as.numeric(OMC_hsdn_AUPR),as.numeric(HGIMC_hsdn_AUPR), IBCF_hsdn_AUPR$V2, LIBMF_hsdn_AUPR$V2),
        main = "hsdn-MechDB", ylab = "AUPR Values", 
        cex.lab = 1.2,
        cex.axis = 1,
        font.axis = 2,
        boxwex = 1,
        las = 2,
        names = c("NMF", "NMFPer","NMFDR", "VDA-GKSBMF", "BNNR", "OMC", "HGIMC", "IBCF", "LIBMF"),
        col = c("lightblue", "lightgreen", "lightpink", "lightgray", "lightyellow", "gold", "deeppink",  "blue", "coral"))


# dev.off()
par(mfrow = c(1,2))

## overall AUC and AUPR for methods
overall_AUC_matrix = rbind( median_AUC_HDVD,  median_AUC_LAGCN  ,
                            median_AUC_Fdata,  median_AUC_Cdata  ,
                            median_AUC_LRSSL,  median_AUC_Ydata  ,
                            median_AUC_oMat,  median_AUC_hsdn )
overall_AUC = colMeans(overall_AUC_matrix)

overall_AUC_SD = apply(overall_AUC_matrix, 2, sd)
overall_AUPR_matrix = rbind( median_AUPR_HDVD,  median_AUPR_LAGCN  ,
                             median_AUPR_Fdata,  median_AUPR_Cdata  ,
                             median_AUPR_LRSSL,  median_AUPR_Ydata  ,
                             median_AUPR_oMat,  median_AUPR_hsdn )
overall_AUPR = colMeans(overall_AUPR_matrix)
overall_AUPR_SD = apply(overall_AUPR_matrix, 2, sd)



labels = c("NMF", "NMF-PDR","NMF-DR", "VDA-GKSBMF", "BNNR", "OMC", "HGIMC", "IBCF", "LIBMF")
Method = c("OMC", "BNNR" ,"HGIMC","VDA-GKSBMF", "NMF-PDR",  "IBCF", "NMF", "NMF-DR",  "LIBMF")

ordered_indices = match(Method, labels)  # Get the indices of the new order
overall_AUC_reordered = overall_AUC[ordered_indices]
overall_AUPR_reordered = overall_AUPR[ordered_indices]
overall_AUC_SD_reordered = overall_AUC_SD[ordered_indices]
overall_AUPR_SD_reordered = overall_AUPR_SD[ordered_indices]

# Data for the table
data <- data.frame(
  Method = Method,
  AUC = overall_AUC_reordered,
  AUPR = overall_AUPR_reordered,
  SD_AUC = overall_AUC_SD_reordered ,
  SD_AUPR = overall_AUPR_SD_reordered 
)


data_melted <- melt(data, id.vars = c("Method","SD_AUC", "SD_AUPR"), variable.name = "Metric", value.name = "Value")

data_melted$Method <- factor(data_melted$Method, levels = unique(data$Method))


# Create grouped bar plot
ggplot(data_melted, aes(x = Method, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = Value - ifelse(Metric == "AUC", SD_AUC, SD_AUPR),
                    ymax = Value + ifelse(Metric == "AUC", SD_AUC, SD_AUPR)), 
                position = position_dodge(width = 0.8), 
                width = 0.25) +
  # Add text labels above bars
  # geom_text(aes(label = round(Value, 2)), 
  #           position = position_dodge(width = 1), 
  #           vjust = -0.5, size = 4) +
  labs(#title = "Overall evaluation values (AUC and AUPR) across methods",
    x = "Methods", cex = 0.8,
    y = "Performance Metric",
    fill = "Metric") +
  theme_minimal() +
  
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),  
    
    axis.title.x = element_text(size = 18, color = "black",hjust = 0.5),  # x-axis label font, face = "bold",
    axis.title.y = element_text(size = 18, color = "black",hjust = 0.5),  # y-axis label font
    
    
    axis.text.x = element_text(size = 16,  angle = 45, hjust = 0.8, color = "black", face ="plain"),  # x-axis ticks font
    axis.text.y = element_text(size = 16, angle = 0, hjust = 0.8, color = "black",face ="plain"),   # y-axis ticks font
    
    legend.title = element_text(size = 16),  # Legend title size
    legend.text = element_text(size = 14)    # Legend text size
    # legend.position = "none"
  )




## overall AUC and AUPR for datasets
labels = c("HDVD", "LAGCN","Fdata", "Cdata", "LRSSL", "Ydata", "oMat-MechDB", "HSDN-MechDB")


overall_AUC_data = c(mean(median_AUC_HDVD) , mean(median_AUC_LAGCN) ,
                     mean(median_AUC_Fdata) , mean(median_AUC_Cdata) ,
                     mean(median_AUC_LRSSL) , mean(median_AUC_Ydata) ,
                     mean(median_AUC_oMat) , mean(median_AUC_hsdn))

overall_AUPR_data = c(mean(median_AUPR_HDVD) , mean(median_AUPR_LAGCN) ,
                      mean(median_AUPR_Fdata) , mean(median_AUPR_Cdata) ,
                      mean(median_AUPR_LRSSL) , mean(median_AUPR_Ydata) ,
                      mean(median_AUPR_oMat) , mean(median_AUPR_hsdn))

SD_AUC_data = c(sd(median_AUC_HDVD), sd(median_AUC_LAGCN), sd(median_AUC_Fdata),
                sd(median_AUC_Cdata), sd(median_AUC_LRSSL), sd(median_AUC_Ydata),
                sd(median_AUC_oMat), sd(median_AUC_hsdn))

SD_AUPR_data = c(sd(median_AUPR_HDVD), sd(median_AUPR_LAGCN), sd(median_AUPR_Fdata),
                 sd(median_AUPR_Cdata), sd(median_AUPR_LRSSL), sd(median_AUPR_Ydata),
                 sd(median_AUPR_oMat), sd(median_AUPR_hsdn))


overall_AUPR_reordered_data =sort(overall_AUPR_data, decreasing = TRUE)
AUPR_order = match(overall_AUPR_reordered_data, overall_AUPR_data)
levels = labels[AUPR_order]
overall_AUC_reordered_data = overall_AUC_data[AUPR_order]
SD_AUC_reordered_data = SD_AUC_data[AUPR_order]
SD_AUPR_reordered_data = SD_AUPR_data[AUPR_order]

# Data for the table
data <- data.frame(
  Datasets = levels,
  AUC = overall_AUC_reordered_data,
  AUPR = overall_AUPR_reordered_data,
  SD_AUC = SD_AUC_reordered_data ,
  SD_AUPR = SD_AUPR_reordered_data 
)


# Melt the data for grouped bar plot
data_melted <- melt(data, id.vars = c("Datasets","SD_AUC", "SD_AUPR"), variable.name = "Metric", value.name = "Value")

# Ensure the order of Method matches the table
data_melted$Datasets <- factor(data_melted$Datasets, levels = unique(data$Datasets))


# Create grouped bar plot
ggplot(data_melted, aes(x = Datasets, y = Value, fill = Metric)) +
  
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  
  geom_errorbar(aes(ymin = Value - ifelse(Metric == "AUC", SD_AUC, SD_AUPR),
                    ymax = Value + ifelse(Metric == "AUC", SD_AUC, SD_AUPR)), 
                position = position_dodge(width = 0.8), 
                width = 0.25) +
  # Add text labels above bars
  # # geom_text(aes(label = round(Value, 2)), 
  #           position = position_dodge(width = 1), 
  #           vjust = ifelse(data_melted$Metric == "AUC", 
  #                          ifelse(data_melted$SD_AUC > 0.085,-4, -2)), 
  #                          #ifelse(data_melted$SD_AUC <= 0.085,-3.5, -2)),
  #           #nudge_y = ifelse(data_melted$Metric == "AUC", 0.05, 0.05),  
  #           size = 5) +
  labs(#title = "Overall evaluation values (AUC and AUPR) across datasets",
    x = "Datasets", cex = 0.8,
    y = "Performance Metric",
    fill = "Metric") +
  theme_minimal() +
  
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),  
    
    axis.title.x = element_text(size = 18, color = "black",hjust = 0.5),  # x-axis label font, face = "bold",
    axis.title.y = element_text(size = 18, color = "black",hjust = 0.5),  # y-axis label font
    
    
    axis.text.x = element_text(size = 16,  angle = 45, hjust = 0.8, color = "black", face ="plain"),  # x-axis ticks font
    axis.text.y = element_text(size = 16, angle = 0, hjust = 0.8, color = "black",face ="plain"),   # y-axis ticks font
    
    legend.title = element_text(size = 16),  # Legend title size
    legend.text = element_text(size = 14)    # Legend text size
    # legend.position = "none"
  )


### Part 4: Compare the different way of cross-validation for HGIMC
# idrug: single = 0.87602, multiple = 0.9307
HGIMC_AUC_multiple_Fdata = read.csv("HGIMC_MultipleSimimarityData/HGIMC_Multiple_Fdata_AUC_25runs.csv",   header = FALSE)
HGIMC_AUC_multiple_Cdata = read.csv("HGIMC_MultipleSimimarityData/HGIMC_Multiple_Cdata_AUC_25runs.csv",   header = FALSE)
HGIMC_AUC_multiple_Ydata = read.csv("HGIMC_MultipleSimimarityData/HGIMC_Multiple_Ydata_AUC_25runs.csv",   header = FALSE)

HGIMC_AUPR_multiple_Fdata = read.csv("HGIMC_MultipleSimimarityData/HGIMC_Multiple_Fdata_AUPR_25runs.csv",   header = FALSE)
HGIMC_AUPR_multiple_Cdata = read.csv("HGIMC_MultipleSimimarityData/HGIMC_Multiple_Cdata_AUPR_25runs.csv",   header = FALSE)
HGIMC_AUPR_multiple_Ydata = read.csv("HGIMC_MultipleSimimarityData/HGIMC_Multiple_Ydata_AUPR_25runs.csv",   header = FALSE)



### Compare HGIMC for single and multiple similarity measures
HGIMC_single_AUC = c(0.8874, 0.9135 ,0.9092 ) # median for 25 runs
HGIMC_multiple_AUC = round(c(median(as.numeric(HGIMC_AUC_multiple_Fdata)), 
                             median(as.numeric(HGIMC_AUC_multiple_Cdata)), 
                             median(as.numeric(HGIMC_AUC_multiple_Ydata))),4)




par(mfrow = c(1,1))

# Combine the vectors into a matrix for grouped bar plotting
auc_matrix <- rbind(HGIMC_single_AUC, HGIMC_multiple_AUC)

# Define custom labels for the x-axis
labels <- c("Fdata", "Cdata", "Ydata")

# Create a bar plot with bars side-by-side and x-axis labels
bar_positions <- barplot(auc_matrix, beside = TRUE, col = c("deeppink", "skyblue"), names.arg = labels,
                         ylim = c(0.5, 1),xpd=FALSE, ylab = "AUC Value",
                         cex.lab = 1.2, cex.axis = 1.2)
# Add the values on top of each bar
text(x = bar_positions, y = auc_matrix, labels = round(auc_matrix, 3), pos = 3, cex = 1)


# Vectors for AUPR values
HGIMC_single_AUPR <- c(0.1750, 0.1969, 0.1998 )
HGIMC_multiple_AUPR = round(c(median(as.numeric(HGIMC_AUPR_multiple_Fdata)), 
                              median(as.numeric(HGIMC_AUPR_multiple_Cdata)), 
                              median(as.numeric(HGIMC_AUPR_multiple_Ydata))),4)

# Combine the vectors into a matrix for grouped bar plotting
aupr_matrix <- rbind(HGIMC_single_AUPR, HGIMC_multiple_AUPR)

# Define custom labels for the x-axis
labels <- c("Fdata", "Cdata", "Ydata")

# Create a bar plot with bars side-by-side and x-axis labels
bar_positions <- barplot(aupr_matrix, beside = TRUE, col = c("deeppink", "skyblue"), names.arg = labels,
                         ylim = c(0, 0.25),  ylab = "AUPR Value",
                         cex.lab = 1.2, cex.axis = 1.2)
# Add the values on top of each bar
text(x = bar_positions, y = aupr_matrix, labels = round(aupr_matrix, 3), pos = 3, cex = 1)

# # Add a custom legend at the center
# par(mar = c(0, 0, 0, 0))  # Remove margins in the legend area
# plot.new()  # Create a blank canvas
# legend("center", legend = c("Single Similarity", "Multiple Similarities"), 
#        fill = c("deeppink", "skyblue"))





# # Combine the vectors into a matrix for grouped bar plotting
# auc_matrix = rbind(HGIMC_single_AUPR, HGIMC_multiple_AUPR)
# 
# # Define custom labels for the x-axis
# labels = c("Fdata", "Cdata", "Ydata", "iDrug")
# 
# # Create a bar plot with the bars side-by-side and x-axis labels
# bar_positions=barplot(auc_matrix, beside = TRUE, col = c("deeppink", "skyblue"), names.arg = labels,
#                       ylim = c(0, 0.4), main = "Comparison of AUPR Values", ylab = "The AUPR values")
# # Add the values on top of each bar
# text(x = bar_positions, y = auc_matrix, labels = round(auc_matrix, 4), pos = 3, cex = 0.8)
# 
# par(mar = c(0, 0, 0, 0))  # Remove margins in legend area
# plot.new()  # Create a blank canvas
# legend("center", legend = c("Single Similarity", "Multiple Similarities"), col = c("deeppink", "skyblue"),
#        fill = c("deeppink", "skyblue"),ncol = 2)






#### Sparsity 
## Overall data


# Load required packages
library(ggplot2)
library(tidyr)
library(dplyr)

## for each data and each method


# save(AUC_median_matrix, file = "AUC_median_matrix_rerun_neworder.RData")
# save(AUC_sd_matrix, file = "AUC_sd_matrix_rerun_neworder.RData")
# save(AUPR_median_matrix, file = "AUPR_median_matrix_rerun_neworder.RData")
# save(AUPR_sd_matrix, file = "AUPR_sd_matrix_rerun_neworder.RData")



AUC_median_matrix_sparsity <- cbind(median_AUC_LAGCN,
                                    median_AUC_HDVD,
                                    median_AUC_oMat,
                                    
                                    median_AUC_Fdata,
                                    median_AUC_Ydata,
                                    median_AUC_Cdata,
                                    median_AUC_LRSSL,
                                    
                                    median_AUC_hsdn)



# AUC_median_matrix_sparsity <- cbind(AUC_median_matrix[colnames(AUC_median_matrix) == "LAGCN"],
#                                     AUC_median_matrix[colnames(AUC_median_matrix) == "HDVD"],
#                                     AUC_median_matrix[colnames(AUC_median_matrix) == "oMat-MechDB"],
#                                     AUC_median_matrix[colnames(AUC_median_matrix) == "Fdata"],
#                                     AUC_median_matrix[colnames(AUC_median_matrix) == "Ydata"],
#                                     AUC_median_matrix[colnames(AUC_median_matrix) == "Cdata"],
#                                     AUC_median_matrix[colnames(AUC_median_matrix) == "LRSSL"],
#                                     AUC_median_matrix[colnames(AUC_median_matrix) == "HSDN-MechDB"]
#                                     )
colnames(AUC_median_matrix_sparsity) = c("LAGCN", "HDVD", "oMat-MechDB", "Fdata", "Ydata", "Cdata", "LRSSL", "HSDN-MechDB")



AUPR_median_matrix_sparsity <- cbind(median_AUPR_LAGCN,
                                     median_AUPR_HDVD,
                                     median_AUPR_oMat,
                                     median_AUPR_Fdata,
                                     median_AUPR_Ydata,
                                     median_AUPR_Cdata,
                                     median_AUPR_LRSSL,
                                     median_AUPR_hsdn)




rownames(AUPR_median_matrix_sparsity) <- c("NMF", "NMF-PDR","NMF-DR", "VDA-GKSBMF", "BNNR", "OMC", "HGIMC", "IBCF", "LIBMF")
rownames(AUC_median_matrix_sparsity) <- c("NMF", "NMF-PDR","NMF-DR", "VDA-GKSBMF", "BNNR", "OMC", "HGIMC", "IBCF", "LIBMF")
data_name = c("HDVD", "LAGCN", "Fdata", "Cdata", "LRSSL", "Ydata", "oMat-MechDB", "HSDN-MechDB")

sparsity = c(0.9389, 0.8855, 0.9895, 0.9913, 0.9941, 0.9912, 0.979, 0.9952)


sorted_indices = order(sparsity, decreasing = FALSE)

data_name_order = data_name[sorted_indices]
sparsity_order = sparsity[sorted_indices]
x_sparsity = c(1:length(sparsity_order))

# Convert matrix to dataframe and reshape to long format
AUPR_df <- as.data.frame(AUPR_median_matrix_sparsity)
AUPR_df$Method <- rownames(AUPR_df)

AUPR_long <- AUPR_df %>%
  pivot_longer(cols = -Method, names_to = "Dataset", values_to = "AUPR")

# Add sparsity values corresponding to datasets
sparsity_df <- data.frame(Dataset = colnames(AUPR_median_matrix_sparsity), Sparsity = x_sparsity)
AUPR_long <- left_join(AUPR_long, sparsity_df, by = "Dataset")

AUPR_long$Method <- factor(AUPR_long$Method, levels = rownames(AUPR_median_matrix_sparsity))


AUPR_summary <- AUPR_long %>%
  group_by(Sparsity) %>%
  summarise(
    mean_AUPR = mean(AUPR),
    se_AUPR = sd(AUPR) / sqrt(n()),  # Standard error
    lower_CI = mean_AUPR - 1.96 * se_AUPR,  # 95% CI lower bound
    upper_CI = mean_AUPR + 1.96 * se_AUPR ,  # 95% CI upper bound
    dif = upper_CI - lower_CI
  )


data_name_order[data_name_order == "oMat-MechDB"] = "oMat"
data_name_order[data_name_order == "HSDN-MechDB"] = "HSDN"

ggplot(AUPR_long, aes(x = Sparsity, y = AUPR, color = Method, group = Method)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = AUPR_summary, aes(x = Sparsity, ymin = lower_CI, ymax = upper_CI),
              inherit.aes = FALSE, fill = "gray70", alpha = 0.3) +
  geom_point(size = 3, alpha = 0.8) +
  scale_x_continuous(
    # scale_x_log10(
    breaks = x_sparsity,
    labels = paste0(sparsity_order, "\n", data_name_order)
    # labels = sparsity_order,
    # guide = guide_axis(n.dodge = 3)  
  ) +
  #scale_color_brewer(palette = "Set2")+
  scale_color_manual(values = c(
    #"black", "red", "green", "darkgreen", "deeppink", "blue", "skyblue", "yellow", "purple"
    "black", "purple", "deeppink", "red", "#ff7f0e", "#1f77b4", "#2ca02c", "#8c564b", "#bcbd22"
  )) + 
  labs(
    # title = "AUPR vs Sparsity for Different Methods",
    x = "Sparsity",
    y = "AUPR Value"
  ) +
  theme_minimal()+
  # geom_smooth(method="lm", se = TRUE, color = NA)+
  # facet_wrap(~ Method)+
  theme(
    # plot.title = element_text(size = 16, hjust = 0.5),  
    
    axis.title.x = element_text(size = 18, color = "black",hjust = 0.5,face = "bold"),  # x-axis label font, face = "bold",
    axis.title.y = element_text(size = 18, color = "black",hjust = 0.5,face = "bold"),  # y-axis label font
    
    
    axis.text.x = element_text(size = 16,  angle = 90, hjust = 0,face = "bold"),  # x-axis ticks font
    axis.text.y = element_text(size = 16, angle = 0, hjust = 0,face = "bold"),   # y-axis ticks font
    
    # legend.title = element_text(size = 16),  # Legend title size
    # legend.text = element_text(size = 14),    # Legend text size
    # legend.position = "bottom"
  )


# Convert matrix to dataframe and reshape to long format
AUC_df <- as.data.frame(AUC_median_matrix_sparsity)
AUC_df$Method <- rownames(AUC_df)

AUC_long <- AUC_df %>%
  pivot_longer(cols = -Method, names_to = "Dataset", values_to = "AUC")

# Add sparsity values corresponding to datasets
sparsity_df <- data.frame(Dataset = colnames(AUC_median_matrix_sparsity), Sparsity = x_sparsity)
AUC_long <- left_join(AUC_long, sparsity_df, by = "Dataset")
AUC_long$Method <- factor(AUC_long$Method, levels = rownames(AUC_median_matrix_sparsity))

AUC_summary <- AUC_long %>%
  group_by(Sparsity) %>%
  summarise(
    mean_AUC = mean(AUC),
    se_AUC = sd(AUC) / sqrt(n()),  # Standard error
    lower_CI = mean_AUC - 1.96 * se_AUC,  # 95% CI lower bound
    upper_CI = mean_AUC + 1.96 * se_AUC   # 95% CI upper bound
  )



ggplot(AUC_long, aes(x = Sparsity, y = AUC, color = Method, group = Method)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = AUC_summary, aes(x = Sparsity, ymin = lower_CI, ymax = upper_CI),
              inherit.aes = FALSE, fill = "gray70", alpha = 0.3) +
  geom_point(size = 3, alpha = 0.8) +
  scale_x_continuous(
    # scale_x_log10(
    breaks = x_sparsity,
    labels = paste0(sparsity_order, "\n", data_name_order)
    # labels = paste0(x_sparsity, "\n", data_names)
    # guide = guide_axis(n.dodge = 3)  
  ) +
  #scale_color_brewer(palette = "Set2")+
  scale_color_manual(values = c(
    #"black", "red", "green", "darkgreen", "pink", "blue", "skyblue", "yellow", "purple"
    "black", "purple", "deeppink", "red", "#ff7f0e", "#1f77b4", "#2ca02c", "#8c564b", "#bcbd22"
  )) + 
  labs(
    # title = "AUC vs Sparsity for Different Methods",
    x = "Sparsity",
    y = "AUC Value"
  ) +
  theme_minimal()+
  # geom_smooth(method="lm", se = TRUE, color = NA)+
  # facet_wrap(~ Method)+
  theme(
    # plot.title = element_text(size = 16, hjust = 0.5),  
    
    axis.title.x = element_text(size = 18, color = "black",hjust = 0.5, face = "bold"),  # x-axis label font, face = "bold",
    axis.title.y = element_text(size = 18, color = "black",hjust = 0.5,face = "bold"),  # y-axis label font
    
    
    axis.text.x = element_text(size = 16,  angle = 90, hjust = 0,face = "bold"),  # x-axis ticks font
    axis.text.y = element_text(size = 16, angle = 0, hjust = 0,face = "bold"),   # y-axis ticks font
    
    # legend.title = element_text(size = 16),  # Legend title size
    # legend.text = element_text(size = 14),    # Legend text size
    # legend.position = "bottom"
  )





############# NMF Per with different tests
set.seed(2024)
# HDVD
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/HDVDdata_asso.RData")
Asso_matrix = t(HDVDdata_asso)
a = load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_KStest_1run/NMF_Per_Stat_kstest_HDVDdata.RData")
pred = Stat_CV
b = load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_ttest_1run/NMF_Per_HDVDdata_prediction.RData")
pred = Stat_CV
c = load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_wiltest_1run/HDVD_matrix_wilttest_NMFPer.RData")
pred = HDVD

# LAGCN
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/LAGCNdata_asso.RData")
Asso_matrix = LAGCNdata_asso
a = load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_KStest_1run/NMFPer_LAGCNdata_Stat_kstest_prediction.RData")
b = load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_wiltest_1run/NMF_Per_wiltest_LAGCN_prediction.RData")
pred = Stat_CV

# Fdata
load("C:/Trang/KIProjects/ComprehensionDR/Datasets/RDataFiles/Fdata_asso.RData")
Asso_matrix = Fdata_asso
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_KStest_1run/NMF_Per_Stat_kstest_Fdata.RData")
NMF_Per_Fdata_ttest_prediction
b = load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_wiltest_1run/NMF_Per_Stat_wiltest_Fdata.RData")
pred = Stat_CV
# Cdata
pred_kstest =load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_KStest_1run/NMF_Per_Pr_kstest_Cdata.RData")
pred_ttest =load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_ttest_1run/NMF_Per_Cdata_ttest_prediction.RData")
pred_wilcox =load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_ttest_1run/NMF_Per_Cdata_ttest_prediction.RData")
# LRSSL
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_KStest_1run/NMF_Per_Stat_kstest_LRSSLdata.RData")
NMF_Per_LRSSL_ttest_prediction
# Ydata
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_KStest_1run/NMF_Per_Stat_kstest_Ydata.RData")
NMF_Per_Ydata_ttest_prediction
# oMAt
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_KStest_1run/NMF_Per_Stat_kstest_Omatdata.RData")
NMF_Per_Stat_ttest_Omatdata
# hsdn
load("C:/Trang/KIProjects/ComprehensionDR/NMF_NMFPermutation/Result_KStest_1run/NMF_Per_kstest_hsdn_prediction.RData")
NMF_Per_hsdndata_ttest_prediction




inputObs_matrix=Asso_matrix
prediction_matrix=pred



#sort inputObs_matrix by column using the decreasing order by column of prediction_matrix
res=sort_matrix(prediction_matrix,inputObs_matrix) 

sorted_inputObs_matrix=res$y_sorted
sorted_score_matrix=res$score_sorted
sort_index=res$sort_index

tpr_list = NULL
fpr_list = NULL

recall_list = NULL
precision_list = NULL

accuracy_list = NULL
F1_list = NULL

#now compute performance metrics for top k=cutoff rows (which means top k predicted values for individual diseases or drugs in the columns) against the remaining rows
for (cutoff in 1:nrow(sorted_inputObs_matrix)){
  #for (cutoff in 1:length(sorted_inputObs_matrix)){
  P_matrix = sorted_inputObs_matrix[1:cutoff,] #predicted Positives
  N_matrix = sorted_inputObs_matrix[-c(1:cutoff),] #predicted Negatives
  
  TP = sum(P_matrix > 0)
  FP = sum(P_matrix == 0)
  TN = sum(N_matrix == 0)
  FN = sum(N_matrix > 0)
  tpr = TP / (TP + FN)
  fpr = FP / (FP + TN)
  recall = TP / (TP + FN)
  precision = TP / (TP + FP)
  accuracy = (TN + TP) / (TN + TP + FN + FP)
  F1 = (2 * TP) / (2*TP + FP + FN)
  if ((2*TP + FP + FN)==0)  F1 = 0
  
  F1_list=c(F1_list,F1)
  accuracy_list=c(accuracy_list,accuracy)
  
  tpr_list = c(tpr_list,tpr)
  fpr_list = c(fpr_list,fpr)
  recall_list = c(recall_list, recall)
  precision_list = c(precision_list,precision)
}


### now compute AUC and AUPR
library(caTools)

NMF_AUC=trapz(fpr_list,tpr_list)#AUC

NMF_AUC

NMF_AUPR=trapz(recall_list,precision_list)#AUPR
NMF_AUPR




par(mfrow = c(1, 2))
plot(fpr_list,tpr_list, type="l", main=paste0("AUC=",round(NMF_AUC,3)), xlab="FPR (1-specificity)", ylab = "TPR (sensitivity)", cex=2, cex.lab=1.5, cex.axis=1.5, cex.main=2, lwd=2)
abline(0,1,lty=2)
plot(recall_list,precision_list,type="l",main=paste0("AUPR=",round(NMF_AUPR,3)), xlab = "Recall", ylab="Precision", cex=2, cex.lab=1.5, cex.axis=1.5, cex.main=2, lwd=2)
par(mfrow = c(1,1))








