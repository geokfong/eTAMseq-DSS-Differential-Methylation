## Step 7: Differential Methylation Analysis (R Script)
# Install DSS package if necessary
# BiocManager::install("DSS")
library(DSS)
require(bsseq)
path = file.path(system.file(package="DSS"), "extdata")

# Example: Read in the data files for WT and KO conditions
dat_wt_ut_1 = read.table("WT_UT_1.dss.txt", header=TRUE, sep="\t", stringsAsFactors=F)
dat_wt_ut_2 = read.table("WT_UT_2.dss.txt", header=TRUE, sep="\t", stringsAsFactors=F)
dat_wt_bt2_1 = read.table("WT_BT2_1.dss.txt", header=TRUE, sep="\t", stringsAsFactors=F)
dat_wt_bt2_2 = read.table("WT_BT2_2.dss.txt", header=TRUE, sep="\t", stringsAsFactors=F)
dat_wt_bt4_1 = read.table("WT_BT4_1.dss.txt", header=TRUE, sep="\t", stringsAsFactors=F)
dat_wt_bt4_2 = read.table("WT_BT4_2.dss.txt", header=TRUE, sep="\t", stringsAsFactors=F)
dat_ko_ut_1 = read.table("KO_UT_1.dss.txt", header=TRUE, sep="\t", stringsAsFactors=F)
dat_ko_ut_2 = read.table("KO_UT_2.dss.txt", header=TRUE, sep="\t", stringsAsFactors=F)
dat_ko_bt2_1 = read.table("KO_BT2_1.dss.txt", header=TRUE, sep="\t", stringsAsFactors=F)
dat_ko_bt2_2 = read.table("KO_BT2_2.dss.txt", header=TRUE, sep="\t", stringsAsFactors=F)
dat_ko_bt4_1 = read.table("KO_BT4_1.dss.txt", header=TRUE, sep="\t", stringsAsFactors=F)
dat_ko_bt4_2 = read.table("KO_BT4_2.dss.txt", header=TRUE, sep="\t", stringsAsFactors=F)

# Define your sample groups and their respective labels
comparisons <- list(
  list(data = list(dat_wt_ut_1, dat_wt_ut_2, dat_wt_bt4_1, dat_wt_bt4_2),
       labels = c("wt_ut", "wt_ut_2", "wt_bt4", "wt_bt4_2")),
  
  list(data = list(dat_wt_ut_1, dat_wt_ut_2, dat_ko_ut_1, dat_ko_ut_2),
       labels = c("wt_ut", "wt_ut_2", "ko_ut", "ko_ut_2")),
  
  list(data = list(dat_ko_ut_1, dat_ko_ut_2, dat_ko_bt4_1, dat_ko_bt4_2),
       labels = c("ko_ut", "ko_ut_2", "ko_bt4", "ko_bt4_2")),
  
  list(data = list(dat_wt_bt4_1, dat_wt_bt4_2, dat_ko_bt4_1, dat_ko_bt4_2),
       labels = c("wt_bt4", "wt_bt4_2", "ko_bt4", "ko_bt4_2")),
  
  list(data = list(dat_wt_bt2_1, dat_wt_bt2_2, dat_ko_bt2_1, dat_ko_bt2_2),
       labels = c("wt_bt2", "wt_bt2_2", "ko_bt2", "ko_bt2_2")),
  
  list(data = list(dat_wt_ut_1, dat_wt_ut_2, dat_wt_bt2_1, dat_wt_bt2_2),
       labels = c("wt_ut", "wt_ut_2", "wt_bt2", "wt_bt2_2")),
  
  list(data = list(dat_ko_ut_1, dat_ko_ut_2, dat_ko_bt2_1, dat_ko_bt2_2),
       labels = c("ko_ut", "ko_ut_2", "ko_bt2", "ko_bt2_2"))
)

# Define threshold values for the differential analysis
thresholds <- c(1, 0.1 , 0.05, 0.01)

# Loop through each comparison
for (comparison in comparisons) {
    for (thresh in thresholds) {
        BSobj = makeBSseqData(comparison$data, comparison$labels)
        # Perform DML test
        dmlTest = DMLtest(BSobj, group1=comparison$labels[1:2], group2=comparison$labels[3:4])
        dmls = callDML(dmlTest, p.threshold=thresh)
        write.table(dmls, file=paste0("dmls_", comparison$labels[1], "_vs_", comparison$labels[3], "_thresh", thresh, ".tsv"), sep="\t", row.names=FALSE)
  
        # Call DMR
        dmrs = callDMR(dmlTest, p.threshold=thresh)
        write.table(dmrs, file=paste0("dmrs_", comparison$labels[1], "_vs_", comparison$labels[3], "_thresh", thresh, ".tsv"), sep="\t", row.names=FALSE)
  
        # Prepare the data for BedGraph
        bedgraph_data <- data.frame(chr = dmls$chr, start = dmls$pos, end = dmls$pos + 1, diff = dmls$diff)
        write.table(bedgraph_data, file=paste0("dmls_", comparison$labels[1], "_vs_", comparison$labels[3], "_thresh", thresh, ".bedgraph"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  
        bedgraph_data_dmrs <- data.frame(chr = dmrs$chr, start = dmrs$start, end = dmrs$end, diff = dmrs$diff.Methy)
        write.table(bedgraph_data_dmrs, file=paste0("dmrs_", comparison$labels[1], "_vs_", comparison$labels[3], "_thresh", thresh, ".bedgraph"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    }
}
