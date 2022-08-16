data_dir <- normalizePath("/camp/lab/luscomben/home/users/palk/")
project_dir <- file.path(data_dir, "Projects",
    "0003-collaboration-with-slava/")
setwd(project_dir)
input_dir <- file.path(project_dir, "input_files")
save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
temp_dir <- file.path(project_dir, "temp", "get_geo")
analysis_dir <- file.path(project_dir, "analysis", "final")
Gtex_dir <- file.path("/camp/lab/luscomben/home/shared/projects/tf_splicing/gtex8/")
isoform_expression_level <- file.path(Gtex_dir, "isoforms_expression_per_tissue.tsv")
isoform_tissue_specificity <- file.path(Gtex_dir, "isoforms_tissue_specificity.tsv")
gene_tissue_specificity <- file.path(Gtex_dir, "genes_tissue_specificity.tsv")
gene_expression_file <- file.path(Gtex_dir, "genes_expression_per_tissue.tsv")
nontf_metadata_file <- file.path(Gtex_dir, "ensembl99_tsl1_tsl2_tslNA_mane_nontfs_ID_table_filtered.tsv")
tf_metadata_file <- file.path(Gtex_dir, "ensembl99_tsl1_tsl2_tslNA_mane_tfs_ID_table_filtered.tsv")
tf_family_file <- file.path(Gtex_dir, "tf_family_table.tsv")
domainless_dir <- "/camp/lab/luscomben/home/shared/projects/tf_splicing/domains"
# ==========================================================================
# Imports
# ==========================================================================
library("data.table")
library("stringr")
library("ggplot2")
library("rtracklayer")
library("hexbin")
library("grid")

output_prefix = function(){
    return(format(Sys.Date(),"%Y-%m-%d"))
}
source("analysis-scripts-master/src/plot_utilities.R")
# ==========================================================================
# Analysis
# ==========================================================================


#
# Import transcript to gene mapping
# --------------------------------------------------------------------------

tf_metadata_df <- fread(tf_metadata_file, data.table = FALSE)[, -1]
# notf_metadata_df <- fread(nontf_metadata_file, data.table = FALSE)[, -1]
tf_metadata_df$group <- "DBD+"
# notf_metadata_df$group <- "not_transcription_factor"
# all_tf_metadata_df <- rbind(tf_metadata_df, notf_metadata_df)
# all_tf_metadata_df <- tf_metadata_df
tf_family_df <- fread(tf_family_file, data.table = FALSE)
domainless_df <- fread(file.path(domainless_dir, "isoforms_no_domains_list.tsv"), data.table = FALSE)
domainless_df$group <- "Domainless"
dbdless_df <-  fread(file.path(domainless_dir, "isoforms_no_dbd_list.tsv"), data.table = FALSE)
dbdless_df$group <- "DBD-"

dbdneg_df <- rbind(domainless_df, dbdless_df)

tf_metadata_df$group[tf_metadata_df$ensembl_transcript_id %in% dbdneg_df[,1]] <- "DBD-"
# tf_metadata_df$group[tf_metadata_df$ensembl_transcript_id %in% domainless_df[,1]] <- "Domainless"
# tf_metadata_df$group[tf_metadata_df$ensembl_transcript_id %in% dbdless_df[,1]] <- "DBD-"

#
# Import gene expression and isoform expression matrices
# --------------------------------------------------------------------------

isoform_expr_df <- fread(file = isoform_expression_level, data.table = FALSE)
isoform_expr_matrix <- as.matrix(isoform_expr_df[, -1])
rownames(isoform_expr_matrix) <- isoform_expr_df[, 1]

gene_expression_df <- fread(gene_expression_file, data.table = FALSE)
gene_expression_matrix <- as.matrix(gene_expression_df[, -1])
rownames(gene_expression_matrix) <- gene_expression_df[, 1]


isoform_specificity_df <- fread(isoform_tissue_specificity, data.table = FALSE)
gene_specificity_df <- fread(gene_tissue_specificity, data.table = FALSE)

tf_isoform_specificity_df <- isoform_specificity_df[isoform_specificity_df[,1] %in% tf_metadata_df[,2],]
tf_isoform_specificity_df$group <- tf_metadata_df$group[match(tf_isoform_specificity_df$ensembl_transcript_id, 
	tf_metadata_df$ensembl_transcript_id)]

#
# Plot tissue specificity vs isoform expression
# --------------------------------------------------------------------------
# mean_expr_of_isoforms <- rowMeans(isoform_expr_matrix)
mean_expr_of_isoforms <- vapply(seq_len(nrow(isoform_expr_matrix)), function(x){
	a_vector <- isoform_expr_matrix[x,]
	mean(a_vector[a_vector > quantile(a_vector, 0.80)])
}, 1)
isoform_variance <- vapply(seq_len(nrow(isoform_expr_matrix)), function(x){
	var(isoform_expr_matrix[x,])
},1)


tf_expression_df <- cbind(tf_isoform_specificity_df,
	mean_expr_in_tissues = mean_expr_of_isoforms[match(tf_isoform_specificity_df[,1], 
		rownames(isoform_expr_matrix))],
	expr_variance = isoform_variance[match(tf_isoform_specificity_df[,1], rownames(isoform_expr_matrix))])

tf_expression_df$mean_expr_in_tissues[is.na(tf_expression_df$mean_expr_in_tissues)] <- 0
Plot_df <- tf_expression_df


frequency_table_df <- as.data.frame(table(Plot_df$group))
colnames(frequency_table_df) <- c("group", "frequency")
frequency_table_df$label <- paste("N", 
	frequency_table_df$frequency, sep = "=")
frequency_table_df$expr_yaxis <- max(Plot_df$mean_expr_in_tissues)
frequency_table_df$specificity_yaxis <- min(Plot_df$tissue_specificity)


Plot_df$group <- factor(Plot_df$group, levels = c("DBD+", "DBD-"))

default_theme <- .fetch_barplot_theme()
the_plot <- ggplot(Plot_df, aes(x = group, y = mean_expr_in_tissues)) +
geom_violin(aes(fill = group), colour = "#000000") +
geom_boxplot(fill = "#e5ece9", colour = "#000000", width = 0.1, outlier.shape = NA, outlier.size = NA) +
geom_label(data = frequency_table_df, aes(x = group, y = expr_yaxis, label = label)) +
geom_hline(yintercept = 1, linetype = "dashed") +
scale_fill_brewer("group", palette = "Set1") +
scale_x_discrete("Transcription factor group") +
scale_y_continuous("Isoform expression level (mean(x > Q(x, 0.80)))") +
default_theme + theme(legend.position = "bottom") + 
labs(title = "mean isoform expression in domain groups")
ggsave(file = file.path(analysis_dir, paste(output_prefix(), 
	"isoform_expression_distribution_in_different_tf_domain_groups_quantile_0.80.pdf", sep = "-")), 
the_plot, height = unit(9, "cm"), width = unit(5, "cm"))


# Plot_df <- tf_expression_df[!is.nan(tf_expression_df$mean_expr_in_tissues),]
# Plot_df$group <- factor(Plot_df$group, levels = c("DBD+", "DBD-"))

# default_theme <- .fetch_barplot_theme()
# the_plot <- ggplot(Plot_df, aes(x = group, y = tissue_specificity)) +
# geom_violin(aes(fill = group), colour = "#000000") +
# geom_boxplot(fill = "#e5ece9", colour = "#000000", width = 0.1, outlier.shape = NA, outlier.size = NA) +
# geom_label(data = frequency_table_df, aes(x = group, y = specificity_yaxis, label = label), alpha = 0.6) +
# scale_fill_brewer("group", palette = "Set1") +
# scale_x_discrete("Transcription factor group") +
# scale_y_continuous("Tissue specificity") +
# default_theme + theme(legend.position = "bottom") + 
# labs(title = "Tissue specificity distribution in domain groups")
# ggsave(file = file.path(analysis_dir, paste(output_prefix(), 
# 	"tissue_specificity_distribution_in_different_tf_domain_groups.pdf", sep = "-")), 
# the_plot, height = unit(9, "cm"), width = unit(5, "cm"))


#
# Do pairwise wilcoxon test
# --------------------------------------------------------------------------

unique_groups <- unique(tf_expression_df$group)
combinations <- as.data.frame(cbind(rep(seq_along(unique_groups), 
	times = length(unique_groups)), 
	rep(seq_along(unique_groups), each = length(unique_groups))))
combinations <- combinations[combinations[,2] > combinations[,1],]

tissue_specificity_pval <- do.call(rbind, lapply(seq_len(nrow(combinations)), function(x){
	test_object <- wilcox.test(x = tf_expression_df$tissue_specificity[tf_expression_df$group %in% unique_groups[combinations[x,1]]],
		y = tf_expression_df$tissue_specificity[tf_expression_df$group %in% unique_groups[combinations[x,2]]])
	data.frame(pval = test_object$p.value, 
		group1 = unique_groups[combinations[x,1]],
		group2 = unique_groups[combinations[x,2]])
}))

isoform_expression_pval <- do.call(rbind, lapply(seq_len(nrow(combinations)), function(x){
	test_object <- wilcox.test(x = tf_expression_df$mean_expr_in_tissues[tf_expression_df$group %in% unique_groups[combinations[x,1]]],
		y = tf_expression_df$mean_expr_in_tissues[tf_expression_df$group %in% unique_groups[combinations[x,2]]])
	data.frame(pval = test_object$p.value, 
		group1 = unique_groups[combinations[x,1]],
		group2 = unique_groups[combinations[x,2]])
}))


# #
# # Import isoform specificity and gene specificity
# # --------------------------------------------------------------------------
# isoform_specificity_df <- fread(isoform_tissue_specificity, data.table = FALSE)
# gene_specificity_df <- fread(gene_tissue_specificity, data.table = FALSE)


# isoform_specificity_df$group <- all_tf_metadata_df[match(isoform_specificity_df[,1], all_tf_metadata_df[,2]), 3]
# isoform_split_by_group <- split(isoform_specificity_df, isoform_specificity_df$group)

