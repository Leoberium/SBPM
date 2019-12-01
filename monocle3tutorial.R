library(monocle3)
library(dplyr)

#### Getting started with Monocle 3 ####

# The cell_data_set class
expression_matrix <- readRDS(url('http://staff.washington.edu/hpliner/data/cao_l2_expression.rds'))
cell_metadata <- readRDS(url('http://staff.washington.edu/hpliner/data/cao_l2_colData.rds'))
gene_annotation <- readRDS(url('http://staff.washington.edu/hpliner/data/cao_l2_rowData.rds'))

# Make the CDS object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
colData(cds)

#### Clustering and classifying your cells ####

# Pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)

# Reduce dimensionality and visualize the cells
cds <- reduce_dimension(cds, reduction_method = 'UMAP', preprocess_method = 'PCA', cores = 8)
cds <- reduce_dimension(cds, reduction_method = 'tSNE', preprocess_method = 'PCA', cores = 8)
plot_cells(cds, color_cells_by = 'cao_cell_type', group_label_size = 3) # UMAP
plot_cells(cds, genes = c('cpna-2', 'egl-21', 'ram-2', 'inos-1')) # by gene
plot_cells(cds, color_cells_by = 'cao_cell_type', group_label_size = 3, reduction_method = 'tSNE') # tSNE

# Check for and remove batch effects
plot_cells(cds, color_cells_by = 'plate', label_cell_groups = FALSE)
cds <- align_cds(cds, num_dim = 100, alignment_group = 'plate')
cds <- reduce_dimension(cds)
plot_cells(cds, color_cells_by = 'plate', label_cell_groups = FALSE)

# Group cells into clusters
cds <- cluster_cells(cds, resolution=1e-5, reduction_method='UMAP')
plot_cells(cds, group_label_size = 3)
plot_cells(cds, color_cells_by = 'partition', group_cells_by = 'partition')
plot_cells(cds, color_cells_by = 'cao_cell_type')
plot_cells(cds, color_cells_by = 'cao_cell_type', label_groups_by_cluster = FALSE)

# Find marker genes expressed by each cluster
marker_test_res <- top_markers(cds, group_cells_by = 'partition', 
                               reference_cells = 1000, cores = 8)
top_specific_markers <- marker_test_res %>%
                            filter(fraction_expressing >= 0.10) %>%
                            group_by(cell_group) %>%
                            top_n(1, pseudo_R2)
top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by = 'partition',
                    ordering_type = 'maximal_on_diag',
                    max.size = 3)

top_specific_markers <- marker_test_res %>%
                            filter(fraction_expressing >= 0.10) %>%
                            group_by(cell_group) %>%
                            top_n(3, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by = 'partition',
                    ordering_type = 'cluster_row_col',
                    max.size=3)

# Annotate your cells according to type
colData(cds)$assigned_cell_type <- as.character(partitions(cds))
colData(cds)$assigned_cell_type = dplyr::recode(colData(cds)$assigned_cell_type,
                                                '1' = 'Germline',
                                                '2' = 'Body wall muscle',
                                                '3' = 'Unclassified neurons',
                                                '4' = 'Vulval precursors',
                                                '5' = 'Failed QC',
                                                '6' = 'Seam cells',
                                                '7' = 'Pharyngeal epithelia',
                                                '8' = 'Coelomocytes',
                                                '9' = 'Am/PH sheath cells',
                                                '10' = 'Failed QC',
                                                '11' = 'Touch receptor neurons',
                                                '12' = 'Intestinal/rectal muscle',
                                                '13' = 'Pharyngeal neurons',
                                                '14' = 'NA',
                                                '15' = 'flp-1(+) interneurons',
                                                '16' = 'Canal associated neurons',
                                                '17' = 'Ciliated sensory neurons',
                                                '18' = 'Other interneurons',
                                                '19' = 'Pharyngeal gland',
                                                '20' = 'Failed QC',
                                                '21' = 'Ciliated sensory neurons',
                                                '22' = 'Oxygen sensory neurons',
                                                '23' = 'Ciliated sensory neurons',
                                                '24' = 'Ciliated sensory neurons',
                                                '25' = 'Ciliated sensory neurons',
                                                '26' = 'Ciliated sensory neurons',
                                                '27' = 'Oxygen sensory neurons',
                                                '28' = 'Ciliated sensory neurons',
                                                '29' = 'Unclassified neurons',
                                                '30' = 'Socket cells',
                                                '31' = 'Failed QC',
                                                '32' = 'Pharyngeal gland',
                                                '33' = 'Ciliated sensory neurons',
                                                '34' = 'Ciliated sensory neurons',
                                                '35' = 'Ciliated sensory neurons',
                                                '36' = 'Failed QC',
                                                '37' = 'Ciliated sensory neurons',
                                                '38' = 'Pharyngeal muscle')
plot_cells(cds, group_cells_by = 'partition', color_cells_by = 'assigned_cell_type')
cds_subset <- choose_cells(cds)

#### Constructing single-cell trajectories ####

expression_matrix <- readRDS(url('http://staff.washington.edu/hpliner/data/packer_embryo_expression.rds'))
cell_metadata <- readRDS(url('http://staff.washington.edu/hpliner/data/packer_embryo_colData.rds'))
gene_annotation <- readRDS(url('http://staff.washington.edu/hpliner/data/packer_embryo_rowData.rds'))

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

# Pre-process the data

cds <- preprocess_cds(cds = cds, num_dim = 50, method = 'PCA')
cds <- align_cds(cds = cds, alignment_group = 'batch', 
                 residual_model_formula_str = '~ bg.300.loading + bg.400.loading +
                 bg.500.1.loading + bg.500.2.loading + bg.r17.loading +
                 bg.b01.loading + bg.b02.loading')

# Reduce dimensionality and visualize the results
cds <- reduce_dimension(cds = cds, reduction_method = 'UMAP', preprocess_method = 'PCA')
plot_cells(cds = cds, label_groups_by_cluster = FALSE, color_cells_by = 'cell.type')
ciliated_genes <- c('che-1',
                    'hlh-17',
                    'nhr-6',
                    'dmd-6',
                    'ceh-36',
                    'ham-1')
plot_cells(cds = cds, genes = ciliated_genes, label_cell_groups = FALSE,
           show_trajectory_graph = FALSE)

# Cluster your cells
cds <- cluster_cells(cds = cds, reduction_method = 'UMAP')
plot_cells(cds = cds, color_cells_by = 'partition', group_label_size = 4)

# Learn the trajectory graph
cds <- learn_graph(cds = cds)
plot_cells(cds = cds,
           color_cells_by = 'cell.type',
           label_groups_by_cluster = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           group_label_size = 4)

# Order the cells in pseudotime
plot_cells(cds = cds,
           color_cells_by = 'embryo.time.bin',
           label_cell_groups = FALSE,
           label_leaves = TRUE,
           label_branch_points = TRUE,
           graph_label_size = 4)
cds <- order_cells(cds = cds)
plot_cells(cds = cds,
           color_cells_by = 'pseudotime',
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 4)

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin='130-170') {
  cell_ids <- which(colData(cds)[, 'embryo.time.bin'] == time_bin)
  closest_vertex <- cds@principal_graph_aux[['UMAP']]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[['UMAP']])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  return(root_pr_nodes)
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
plot_cells(cds = cds,
           color_cells_by = 'pseudotime',
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 4)

# Working with 3D trajectories
cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))
cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="partition")
