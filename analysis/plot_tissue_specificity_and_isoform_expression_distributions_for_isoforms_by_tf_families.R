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
library("ggforce")
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
tf_family_df <- fread(tf_family_file, data.table = FALSE)


tf_metadata_df <- fread(tf_metadata_file, data.table = FALSE)[, -1]

tf_metadata_df$group <- "DBD+"
tf_family_df <- fread(tf_family_file, data.table = FALSE)
domainless_df <- fread(file.path(domainless_dir, "isoforms_no_domains_list.tsv"), data.table = FALSE)
domainless_df$group <- "Domainless"
dbdless_df <-  fread(file.path(domainless_dir, "isoforms_no_dbd_list.tsv"), data.table = FALSE)
dbdless_df$group <- "DBD-"
dbdneg_df <- rbind(domainless_df, dbdless_df)

tf_metadata_df$group[tf_metadata_df$ensembl_transcript_id %in% dbdneg_df[,1]] <- "DBD-"

# tf_metadata_df$group[tf_metadata_df$ensembl_transcript_id %in% domainless_df[,1]] <- "Domainless"
# tf_metadata_df$group[tf_metadata_df$ensembl_transcript_id %in% dbdless_df[,1]] <- "DBD-"

tf_family_df$domain_group <- tf_metadata_df$group[match(tf_family_df$ensembl_transcript_id, tf_metadata_df$ensembl_transcript_id)]
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

tf_isoform_specificity_df <- isoform_specificity_df[isoform_specificity_df[,1] %in% tf_family_df[,5],]
tf_isoform_specificity_df$ensembl_gene_id <- tf_metadata_df$ensembl_gene_id[match(tf_isoform_specificity_df$ensembl_transcript_id, tf_metadata_df$ensembl_transcript_id)]
tf_isoform_specificity_df$family <- tf_family_df$tf_family[match(tf_isoform_specificity_df$ensembl_transcript_id, tf_family_df$ensembl_transcript_id)]
tf_isoform_specificity_df$domain_group <- tf_family_df$domain_group[match(tf_isoform_specificity_df$ensembl_transcript_id, tf_family_df$ensembl_transcript_id)]
#
# Plot tissue specificity vs isoform expression
# --------------------------------------------------------------------------
# mean_expr_of_isoforms <- rowMeans(isoform_expr_matrix)
mean_expr_of_isoforms <- vapply(seq_len(nrow(isoform_expr_matrix)), function(x){
	a_vector <- isoform_expr_matrix[x,]
	a_vector[is.nan(a_vector)] <- 0
	mean(a_vector[a_vector > quantile(a_vector, 0.80)])
}, 1)
isoform_variance <- vapply(seq_len(nrow(isoform_expr_matrix)), function(x){
	var(isoform_expr_matrix[x,])
},1)


tf_expression_df <- cbind(tf_isoform_specificity_df,
	mean_expr_in_tissues = mean_expr_of_isoforms[match(tf_isoform_specificity_df[,1], 
		rownames(isoform_expr_matrix))],
	expr_variance = isoform_variance[match(tf_isoform_specificity_df[,1], rownames(isoform_expr_matrix))])

unique_isoform_family_pair <- unique(tf_expression_df[,c(1,5,4)])

candidate_tf_families_counts <- t(table(unique_isoform_family_pair$family, unique_isoform_family_pair$domain_group) > 10)["DBD-",]
candidate_tf_families <- names(candidate_tf_families_counts[candidate_tf_families_counts])


Plot_df <- tf_expression_df[!is.nan(tf_expression_df$mean_expr_in_tissues) & 
	tf_expression_df$family %in% candidate_tf_families,]



frequency_table_df <- as.data.frame(table(Plot_df$family))
colnames(frequency_table_df) <- c("group", "frequency")
frequency_table_df$label <- paste("N", 
	frequency_table_df$frequency, sep = "=")
frequency_table_df$expr_yaxis <- max(Plot_df$mean_expr_in_tissues)
frequency_table_df$specificity_yaxis <- min(Plot_df$tissue_specificity)


Plot_df$domain_group <- factor(Plot_df$domain_group, levels = c("DBD+", "DBD-"))
Plot_df$family <- factor(Plot_df$family, levels = c("C2H2 ZF", "CxxC ZF", "Forkhead", "Homeodomain", "Nuclear receptor", "bHLH", "bZIP"))

default_theme <- .fetch_barplot_theme()
the_plot <- ggplot(Plot_df, aes(x = family, y = mean_expr_in_tissues)) +
geom_violin(aes(fill = domain_group), colour = "#000000", position = position_dodge(0.7), alpha = 0.6) +
# geom_sina(shape = 21, aes(fill = domain_group),  
# 	alpha = 0.3, position = position_dodge(0.7), size = 0.5) +
geom_boxplot(aes(fill = domain_group), colour = "#000000", 
	width = 0.1, outlier.shape = NA, outlier.size = NA, position = position_dodge(0.7), alpha = 0.8) +
geom_hline(yintercept = 1, linetype = "dashed") +
geom_label(data = frequency_table_df, aes(x = group, y = expr_yaxis, label = label)) +
scale_fill_brewer("DBD group", palette = "Set1") +
# scale_colour_brewer("DBD group", palette = "Set1") +
scale_x_discrete("Transcription factor family") +
scale_y_continuous("Isoform expression level (mean(x > Q(x, 0.80)))") +
# facet_grid(. ~ domain_group, scales = "free") +
default_theme + theme(legend.position = "bottom",
	axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2)) + 
labs(title = "mean isoform expression in domain groups")
ggsave(file = file.path(analysis_dir, paste(output_prefix(), 
	"isoform_expression_distribution_in_different_tf_families_quantile_0.80.pdf", sep = "-")), 
the_plot, height = unit(5, "cm"), width = unit(11, "cm"))


# Plot_df <- tf_expression_df[!is.nan(tf_expression_df$mean_expr_in_tissues),]

default_theme <- .fetch_barplot_theme()
the_plot <- ggplot(Plot_df, aes(x = family, y = tissue_specificity)) +
geom_violin(aes(fill = domain_group), colour = "#000000", position = position_dodge(0.7), alpha = 0.6) +
# geom_sina(shape = 21, aes(fill = domain_group), alpha = 0.3, position = position_dodge(0.7), size = 0.5) +
geom_boxplot(aes(fill = domain_group), colour = "#000000", 
	width = 0.1, outlier.shape = NA, outlier.size = NA, position = position_dodge(0.7), alpha = 0.8) +
geom_label(data = frequency_table_df, aes(x = group, y = specificity_yaxis, label = label), alpha = 0.6) +
scale_fill_brewer("group", palette = "Set1") +
# scale_colour_brewer("DBD group", palette = "Set1") +
scale_x_discrete("Transcription factor family") +
scale_y_continuous("Tissue specificity") +
# facet_grid(. ~ domain_group, scales = "free") +
default_theme + theme(legend.position = "bottom",
	axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2)) + 
labs(title = "Tissue specificity distribution in tf families")
ggsave(file = file.path(analysis_dir, paste(output_prefix(), 
	"tissue_specificity_distribution_in_different_tf_families.pdf", sep = "-")), 
the_plot, height = unit(5, "cm"), width = unit(11, "cm"))






