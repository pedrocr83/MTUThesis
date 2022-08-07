# HELPER FUNCTIONS

phyloseq_aggregation <- function(phyloseq_obj, taxa_rank)
{
    if (tolower(taxa_rank) == 'phylum') {
        phyloseq_obj <- tax_glom(phyloseq_obj, taxrank = "Phylum")
    } else if (tolower(taxa_rank) == 'class') {
        phyloseq_obj <- tax_glom(phyloseq_obj, taxrank = "Class")
    } else if  (tolower(taxa_rank) == 'order') {
        phyloseq_obj <- tax_glom(phyloseq_obj, taxrank = "Order")
    } else if  (tolower(taxa_rank) == 'family') {
        phyloseq_obj <- tax_glom(phyloseq_obj, taxrank = "Family")
    } else if  (tolower(taxa_rank) == 'genus') {
        phyloseq_obj <- tax_glom(phyloseq_obj, taxrank = "Genus")
    }

    return(phyloseq_obj)
}

phyloseq_transformation <- function(phyloseq_obj, transformation)
{
    if (tolower(transformation) == 'clr') {
        phyloseq_obj <- microbiome::transform(phyloseq_obj, 'clr')
    } else if (tolower(transformation) == 'relative') {
        phyloseq_obj <- microbiome::transform(phyloseq_obj, 'compositional')
    } else if  (tolower(transformation) == 'log10') {
        phyloseq_obj <- microbiome::transform(phyloseq_obj, 'log10')
    } else if  (tolower(transformation) == 'counts') {
        phyloseq_obj <- microbiome::transform(phyloseq_obj, 'identity')
    }

    return(phyloseq_obj)
}

phyloseq_results <- function(phyloseq_obj, taxa_rank, transformation, class_variable, phasetxt, betatest)
{

    results_name <- paste(phasetxt, taxa_rank, transformation, sep='_')
    scene <- list(camera = list(eye = list(x=1.5, y=1.5, z = 1.5)))
    mds_path <- file.path(results_path, paste(paste('MDS', results_name, sep='_'), 'png', sep='.'))
    tsne_path <- file.path(results_path, paste(paste('TSNE', results_name, sep='_'), 'png', sep='.'))
    umap_path <- file.path(results_path, paste(paste('UMAP', results_name, sep='_'), 'png', sep='.'))


    mds_obj_path <- file.path(temp_files, paste(paste('MDS', results_name, sep='_'), 'RDS', sep='.'))
    tsne_obj_path <- file.path(temp_files, paste(paste('TSNE', results_name, sep='_'), 'RDS', sep='.'))
    umap_obj_path <- file.path(temp_files, paste(paste('UMAP', results_name, sep='_'), 'RDS', sep='.'))

    phyloseq_obj.data <- sample_data(phyloseq_obj) %>% data.frame() %>% rownames_to_column("sample")
    colnames(phyloseq_obj.data) <- c('sample', 'VAR')

    phyloseq_obj.distance <- phyloseq::distance(phyloseq_obj, 'bray')

    permanova <- adonis2(phyloseq_obj.distance ~ VAR, data = phyloseq_obj.data, permutations = 10000)
    beta <- betadisper(phyloseq_obj.distance, phyloseq_obj.data[['VAR']])
    dispersion <- permutest(beta)
    betatest[phasetxt] <- c(permanova$`Pr(>F)`[1], dispersion$tab$`Pr(>F)`[1])

    # MDS RESULTS
    set.seed(39)
    phyloseq_obj.mds <- vegan::metaMDS(phyloseq_obj.distance, k = 3, try = 200, trymax = 200)
    phyloseq_obj.mds <- phyloseq_obj.mds$points %>% data.frame() %>% rownames_to_column("sample")
    colnames(phyloseq_obj.mds) <- c('sample', 'MDS1', 'MDS2', 'MDS3')
    phyloseq_obj.mds <- phyloseq_obj.mds %>% left_join(phyloseq_obj.data, by = "sample")

    result1 <- plot_ly(phyloseq_obj.mds, x=~MDS1,y=~MDS2,z=~MDS3) %>% add_markers(color=~VAR, colors = c("darkgreen", "red"), size=3)
    result1 <- result1  %>% layout(scene = scene)

    saveRDS(result1, file = mds_obj_path)

    # plotly_IMAGE(result1, format = "png", out_file = mds_path, scale = 2)

    # TSNE RESULTS
    set.seed(39)
    perplexity = floor((nrow(as.matrix(phyloseq_obj.distance)) - 1) / 3)
    phyloseq_obj.tsne <- Rtsne(phyloseq_obj.distance, is_distance = TRUE, dims=3, perplexity=perplexity)
    phyloseq_obj.tsne <- data.frame(phyloseq_obj.tsne$Y) %>% bind_cols(as.vector(sample_data(phyloseq_obj)[[class_variable]]))
    colnames(phyloseq_obj.tsne) <- c('tSNE1', 'tSNE2', 'tSNE3', 'VAR')

    result2 <- plot_ly(phyloseq_obj.tsne, x=~tSNE1,y=~tSNE2,z=~tSNE3) %>% add_markers(color=~VAR, colors = c("darkgreen", "red"), size=3)
    result2 <- result2  %>% layout(scene = scene)

    saveRDS(result2, file = tsne_obj_path)

    # plotly_IMAGE(result2, format = "png", out_file = tsne_path, scale = 2)

    # UMAP RESULTS
    set.seed(39)
    if (nrow(t(otu_table(phyloseq_obj))) < 15){
        nneighbors <- nrow(t(otu_table(phyloseq_obj))) - 1
        umap <- umap(t(otu_table(phyloseq_obj)), n_components=3, n_neighbors=nneighbors, metric='pearson')
    }else{
        umap <- umap(t(otu_table(phyloseq_obj)), n_components=3, metric='pearson', random_state=39)
    }
    phyloseq_obj.umapdata <- umap$layout %>% data.frame() %>% rownames_to_column("sample")
    colnames(phyloseq_obj.umapdata) <- c('sample', 'UMAP1', 'UMAP2', 'UMAP3')
    phyloseq_obj.umapdata <- phyloseq_obj.umapdata %>% left_join(phyloseq_obj.data, by = "sample")

    result3 <- plot_ly(phyloseq_obj.umapdata, x=~UMAP1,y=~UMAP2,z=~UMAP3) %>% add_markers(color=~VAR, colors = c("darkgreen", "red"), size=3)
    result3 <- result3  %>% layout(scene = scene)

    saveRDS(result3, file = umap_obj_path)

    # plotly_IMAGE(result3, format = "png", out_file = umap_path, scale = 2)

    return (betatest)

}

phyloseq_diff_tree <- function(class1_obj, class2_obj, class1_str, class2_str, vars_list)
{

    blank_theme <- theme(plot.background=element_blank(),panel.grid.major = element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank(),panel.background = element_blank())
    cccolors <- c("lightblue", "#e0e0e0", "#FFCCCB")

    class1 <- tax_glom(class1_obj, taxrank = "Genus")
    tax_table(class1) <- tax_table(class1) %>% gsub("g__", "", .) %>% gsub("k__", "", .) %>% gsub("p__", "", .) %>% gsub("c__", "", .) %>% gsub("o__", "", .) %>% gsub("f__", "", .)
    class1.metacoder <- parse_phyloseq(class1)
    class1.metacoder$data$otu_table <- zero_low_counts(class1.metacoder, data = "otu_table", min_count = 5)
    no_reads <- rowSums(class1.metacoder$data$otu_table[, class1.metacoder$data$sample_data$sample_id]) == 0
    class1.metacoder <- filter_obs(class1.metacoder, data = "otu_table", ! no_reads, drop_taxa = TRUE)
    class1.metacoder$data$tax_abund <- calc_taxon_abund(class1.metacoder, "otu_table", cols = class1.metacoder$data$sample_data$sample_id)

    class2 <- tax_glom(class2_obj, taxrank = "Genus")
    tax_table(class2) <- tax_table(class2) %>% gsub("g__", "", .) %>% gsub("k__", "", .) %>% gsub("p__", "", .) %>% gsub("c__", "", .) %>% gsub("o__", "", .) %>% gsub("f__", "", .)
    class2.metacoder <- parse_phyloseq(class2)
    class2.metacoder$data$otu_table <- zero_low_counts(class2.metacoder, data = "otu_table", min_count = 5)
    no_reads <- rowSums(class2.metacoder$data$otu_table[, class2.metacoder$data$sample_data$sample_id]) == 0
    class2.metacoder <- filter_obs(class2.metacoder, data = "otu_table", ! no_reads, drop_taxa = TRUE)
    class2.metacoder$data$tax_abund <- calc_taxon_abund(class2.metacoder, "otu_table", cols = class2.metacoder$data$sample_data$sample_id)

    for (var in vars_list)
        {
            class1.metacoder$data$diff_table <- compare_groups(class1.metacoder, data = "tax_abund",
                                                               cols = class1.metacoder$data$sample_data$sample_id,
                                                               groups = class1.metacoder$data$sample_data[[var]],
                                                               combinations = list(c("notanomaly", "anomaly")))
            class2.metacoder$data$diff_table <- compare_groups(class2.metacoder, data = "tax_abund",
                                                               cols = class2.metacoder$data$sample_data$sample_id,
                                                               groups = class2.metacoder$data$sample_data[[var]],
                                                               combinations = list(c("notanomaly", "anomaly")))

            class1_tree_path <- file.path(results_path, paste(paste(class1_str, var, sep='_'), 'png', sep='.'))
            class2_tree_path <- file.path(results_path, paste(paste(class2_str, var, sep='_'), 'png', sep='.'))
            set.seed(39)

            heat_tree(class1.metacoder, node_label = taxon_names,
                      node_size = n_obs, node_color = mean_diff,
                      node_color_interval = c(-20, 20), node_color_range = cccolors,
                      make_node_legend = FALSE, make_edge_legend = FALSE,
                      layout = "davidson-harel", initial_layout = "reingold-tilford",
                      output_file = class1_tree_path, background= "white")

            heat_tree(class2.metacoder, node_label = taxon_names,
                      node_size = n_obs, node_color = mean_diff,
                      node_color_interval = c(-20, 20), node_color_range = cccolors,
                      make_node_legend = FALSE, make_edge_legend = FALSE,
                      layout = "davidson-harel", initial_layout = "reingold-tilford",
                      output_file = class2_tree_path, background= "white")
        }
}

# IMPORT LIBRARIES

suppressMessages(library(phyloseq))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(microbiome))
suppressMessages(library(metacoder))
suppressMessages(library(GUniFrac))
suppressMessages(library(microViz))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(getopt))
suppressMessages(library(tools))
suppressMessages(library(glue))
suppressMessages(library(vegan))
suppressMessages(library(Rtsne))
suppressMessages(library(scater))
suppressMessages(library(umap))
suppressMessages(library(plotly))
suppressMessages(library(plotly.microbiome))

Sys.setenv("plotly_username" = "pedrocr83")
Sys.setenv("plotly_api_key" = "j3lrXQabnJ8X27WCZxrt")

# SET THEME
theme_set(theme_bw())

# 0. Get CLI parameters:

spec <- matrix(c(
  'study_name'           , 's', 1, "character", "STUDY NAME  (required)",
  'biom_file'            , 'b', 1, "character", "BIOM or OTU FILE PATH  (required)",
  'mapping_file'         , 'm', 1, "character", "STUDY META MAPPING FILE PATH (with headers) (required)",
  'taxonomy_file'        , 't', 1, "character", "GREENGENES TAXONOMY FILE PATH (required)",
  'tree_file'            , 'e', 1, "character", "GREENGENES TREE FILE PATH (required)",
  'class_variable'       , 'v', 1, "character", "STUDY META MAPPING COLUMN USED BY CLASSES (required)",
  'classes'              , 'c', 1, "character", "STUDY TWO CLASSES USED (SEPARATED BY ',') (required)",
  'taxa_rank'            , 'r', 1, "character", "TAXA RANK AGGREGATION (phylum, genus)",
  'transformation'       , 'a', 1, "character", "OTU ABUNDANCE TRANSFORMATION (relative, clr , log10)",
  'distance_measure'     , 'u', 1, "character", "DISTANCE TO USE FOR ML ANALYSIS (gunifrac, unifrac, wunifrac",
  'phase1'               , 'p', 1, "logical", "TRUE/FALSE",
  'phase2'               , 'q', 1, "logical", "TRUE/FALSE",
  'help'                 , 'h', 0, "logical",   "this help"
),ncol=5,byrow=T)

opt = getopt(spec);
options(warn=-1)

if (!is.null(opt$help) || is.null(opt$study_name) || is.null(opt$biom_file) || is.null(opt$mapping_file)
|| is.null(opt$taxonomy_file) || is.null(opt$tree_file) || is.null(opt$class_variable) || is.null(opt$classes)
|| is.null(opt$taxa_rank) || is.null(opt$transformation) || is.null(opt$distance_measure)
|| is.null(opt$phase1) || is.null(opt$phase2)) {
  cat(paste(getopt(spec, usage=T),"\n"));
  q();
}

# IMPORT AND PREPROCESS DATASET

cat('\nIMPORT AND PREPROCESS DATASET\n')

study_name <- opt$study_name
study_biom <- opt$biom_file
mapping_file <- opt$mapping_file
taxonomy_file <- opt$taxonomy_file
tree_file <- opt$tree_file
class_variable <- opt$class_variable
classes <- c(unlist(strsplit(opt$classes, ",")))
class1 <- classes[1]
class2 <- classes[2]

phase1 <- opt$phase1
phase2 <- opt$phase2

taxa_rank <- opt$taxa_rank
transformation <- opt$transformation
distance_measure <- opt$distance_measure
distance_measure_results <- opt$distance_measure_results

run_name <- paste(taxa_rank, transformation, distance_measure, sep='_')
run_name_path <- file.path('results', study_name, run_name)
results_path <- file.path('results', study_name, run_name, 'results')
temp_files <- file.path('results', study_name, run_name, 'tempfiles')

phase1.class1.otus.path <- file.path(temp_files, paste(study_name, taxa_rank, transformation, distance_measure, 'class1_otus.csv', sep='_'))
phase1.class2.otus.path <- file.path(temp_files, paste(study_name, taxa_rank, transformation, distance_measure, 'class2_otus.csv', sep='_'))
phase1.class1.distances.path <- file.path(temp_files, paste(study_name, taxa_rank, transformation, distance_measure, 'class1_distances.csv', sep='_'))
phase1.class2.distances.path <- file.path(temp_files, paste(study_name, taxa_rank, transformation, distance_measure, 'class2_distances.csv', sep='_'))

phase1.biom_obj.path <- file.path(temp_files, paste(study_name, taxa_rank, transformation, distance_measure, 'biom.RDS', sep='_'))
phase1.class1.biom_obj.path <- file.path(temp_files, paste(study_name, taxa_rank, transformation, distance_measure, 'class1_biom.RDS', sep='_'))
phase1.class2.biom_obj.path <- file.path(temp_files, paste(study_name, taxa_rank, transformation, distance_measure, 'class2_biom.RDS', sep='_'))

phase2.class1.anomalyresults <- file.path(results_path, 'anomaly_class1.csv')
phase2.class2.anomalyresults <- file.path(results_path, 'anomaly_class2.csv')

betatest.results.path <- file.path(results_path, 'betatest.csv')

working_directory <- getwd()
glue("Current working dir: {working_directory}")

glue("Taxonomy file: {taxonomy_file}")
glue("Tree file: {tree_file}")

glue("FILES: {study_biom} (biom/otu); {mapping_file} (meta file)")
glue("VARIABLES: {class_variable} (class variable); {classes} (class); {opt$taxa_rank} (rank)")

if (phase1 == TRUE) {

    dir.create(run_name_path)
    dir.create(results_path)
    dir.create(temp_files)

    if (file_ext(study_biom) == 'biom') {
    study.otu <- import_biom(study_biom)
    } else {
    study.otu <- read.csv(study_biom, header=T, sep='\t', row.names=1)
    }

    study.mapping <- read.csv(mapping_file, header=T, sep='\t')
    study.taxonomy <- read.csv(taxonomy_file, header=F, sep='\t')
    study.tree <- read_tree(tree_file)

    betatest <- data.frame(matrix(ncol = 0, nrow = 2))
    row.names(betatest) <- c('permanova', 'betapermutest')

    row.names(study.taxonomy) <- study.taxonomy$V1
    study.taxonomy <- study.taxonomy %>% separate(V2, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=';', remove=T)
    study.taxonomy <- subset(study.taxonomy, select = -c(V1))

    row.names(study.mapping) <- study.mapping$X.SampleID
    study.samples <- row.names(study.mapping[study.mapping[[class_variable]] %in% classes, ])

    study.samples.class1 <- row.names(study.mapping[study.mapping[[class_variable]] == class1, ])
    study.samples.class2 <- row.names(study.mapping[study.mapping[[class_variable]] == class2, ])

    study.mapping <- subset(study.mapping, select = -c(X.SampleID))
    study <- phyloseq(otu_table(as.matrix(study.otu), taxa_are_rows=T), sample_data(study.mapping), tax_table(as.matrix(study.taxonomy)), study.tree)
    study <- subset_samples(study, (row.names(sample_data(study)) %in% study.samples))

    # REMOVE TAXA WITH ==0 SUM ACROSS SAMPLES AND CLASSES
    study <- prune_taxa(taxa_sums(study) > 0, study)

    study.class1 <- subset_samples(study, (row.names(sample_data(study)) %in% study.samples.class1))
    study.class2 <- subset_samples(study, (row.names(sample_data(study)) %in% study.samples.class2))

    study.class1 <- prune_taxa(taxa_sums(study.class1) > 0, study.class1)
    study.class2 <- prune_taxa(taxa_sums(study.class2) > 0, study.class2)

    # SAVE UNIFRAC DISTANCES AS CSV FOR ML ANALYSIS
    cat('\nSAVE DISTANCES AS CSV FOR ML ANALYSIS\n')
    # CLASS 1
    study.class1.distances <- dist_get(dist_calc(study.class1, dist=distance_measure))
    study.class1.distances <- melt(as.matrix(study.class1.distances), varnames = c("row", "col"))
    study.class1.distances <- acast(study.class1.distances, row ~ col)
    write.csv(study.class1.distances, file = phase1.class1.distances.path)

    # CLASS 2
    study.class2.distances <- dist_get(dist_calc(study.class2, dist=distance_measure))
    study.class2.distances <- melt(as.matrix(study.class2.distances), varnames = c("row", "col"))
    study.class2.distances <- acast(study.class2.distances, row ~ col)
    write.csv(study.class2.distances, file = phase1.class2.distances.path)

    saveRDS(study, file = phase1.biom_obj.path)
    saveRDS(study.class1, file = phase1.class1.biom_obj.path)
    saveRDS(study.class2, file = phase1.class2.biom_obj.path)

    # PHASE 1 PROCESSING AND RESULTS

    cat('\nPHASE 1 PROCESSING AND RESULTS\n')

    cat('\nAGGREGATION\n')
    study.phase1 <- phyloseq_aggregation(study, taxa_rank)
    cat('\nTRANSFORMATION\n')
    study.phase1 <- phyloseq_transformation(study.phase1, transformation)

    betatest <- phyloseq_results(study.phase1, taxa_rank, transformation, class_variable, 'phase1', betatest)

    # SAVE OTU VALUES AS CSV FOR ML ANALYSIS
    cat('\nSAVE OTU VALUES AS CSV FOR ML ANALYSIS\n')

    study.phase1.class1 <- subset_samples(study.phase1, (row.names(sample_data(study)) %in% study.samples.class1))
    study.phase1.class2 <- subset_samples(study.phase1, (row.names(sample_data(study)) %in% study.samples.class2))

    # CLASS 1
    study.class1.otus <- as(otu_table(study.phase1.class1), "matrix")
    if(taxa_are_rows(study.class1)){study.class1.otus <- t(study.class1.otus)}
    study.class1.otus <- as.data.frame(study.class1.otus)
    write.csv(study.class1.otus, file = phase1.class1.otus.path)

    # CLASS 2
    study.class2.otus <- as(otu_table(study.phase1.class2), "matrix")
    if(taxa_are_rows(study.class2)){study.class2.otus <- t(study.class2.otus)}
    study.class2.otus <- as.data.frame(study.class2.otus)
    write.csv(study.class2.otus, file = phase1.class2.otus.path)

    write.csv(betatest, file = betatest.results.path)

    cat('\nFINISHED PHASE 1\n')

} else if (phase2 == TRUE) {

    # PHASE 2 PROCESSING AND RESULTS

    study.class1 <- readRDS(phase1.class1.biom_obj.path)
    study.class2 <- readRDS(phase1.class2.biom_obj.path)

    cat('\nPHASE 2 PROCESSING AND RESULTS\n')

    class1.anomaly.mapping <- read.csv(phase2.class1.anomalyresults, header=T, sep=',')
    class2.anomaly.mapping <- read.csv(phase2.class2.anomalyresults, header=T, sep=',')

    row.names(class1.anomaly.mapping) <- class1.anomaly.mapping$id
    class1.anomaly.mapping <- subset(class1.anomaly.mapping, select = -c(id))
    row.names(class2.anomaly.mapping) <- class2.anomaly.mapping$id
    class2.anomaly.mapping <- subset(class2.anomaly.mapping, select = -c(id))

    # ADD ALL CLASSIFICATION COLUMNS TO SAMPLE_DATA CLASS 1 AND CLASS 2 and CALCULATE TREES

    # CLASS 1
    study.class1.mapping <- merge(as.data.frame(sample_data(study.class1)), class1.anomaly.mapping, by=0, all=TRUE)
    row.names(study.class1.mapping) <- study.class1.mapping$Row.names
    study.class1.mapping <- subset(study.class1.mapping, select = -c(Row.names))
    sample_data(study.class1) <- sample_data(study.class1.mapping)

    # CLASS 2
    study.class2.mapping <- merge(as.data.frame(sample_data(study.class2)), class2.anomaly.mapping, by=0, all=TRUE)
    row.names(study.class2.mapping) <- study.class2.mapping$Row.names
    study.class2.mapping <- subset(study.class2.mapping, select = -c(Row.names))
    sample_data(study.class2) <- sample_data(study.class2.mapping)

    # DIFF TREE
    tree_graph_cols <- colnames(class1.anomaly.mapping)
    phyloseq_diff_tree(study.class1, study.class2, class1, class2, tree_graph_cols)

    col_df <- data.frame(class1=colnames(class1.anomaly.mapping), class2=colnames(class2.anomaly.mapping))
    betatest <- read.csv(betatest.results.path, header=T, sep=',')

    for(i in 1:nrow(col_df)) {

        file_name <- col_df[i,"class1"]

        cat("Processing" , file_name ,"\n")

        study <- readRDS(phase1.biom_obj.path)
        class1.anomaly.anomalies <- row.names(class1.anomaly.mapping[class1.anomaly.mapping[[col_df[i,'class1']]] == 'anomaly', ])
        class2.anomaly.anomalies <- row.names(class2.anomaly.mapping[class2.anomaly.mapping[[col_df[i,'class2']]] == 'anomaly', ])

        anomalies <- c(class1.anomaly.anomalies, class2.anomaly.anomalies)
        study <- subset_samples(study, !(row.names(sample_data(study)) %in% anomalies))
        study <- prune_taxa(taxa_sums(study) > 0, study)

        cat('\nAGGREGATION\n')
        study.phase1 <- phyloseq_aggregation(study, taxa_rank)
        cat('\nTRANSFORMATION\n')
        study.phase1 <- phyloseq_transformation(study.phase1, transformation)

        betatest <- phyloseq_results(study.phase1, taxa_rank, transformation, class_variable, paste(file_name, 'phase2', sep='_'), betatest)

    }

    write.csv(betatest, file = betatest.results.path)

} else {

    print("PLEASE CHOOSE PHASE 1 OR PHASE 2!")

}