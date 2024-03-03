import numpy as np
import os
import pandas as pd
import scanpy as sc


sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

absolute_path_to_root = "C:/Users/Lenovo/Desktop/Dr_Rueda_Research/scanpy_tutorial/"

results_file_path = absolute_path_to_root + "Results/write"
results_file = results_file_path + "pbmc3k.h5ad"  # the file that will store the analysis results

data_file_path = absolute_path_to_root + "data/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19"
adata = sc.read_10x_mtx(data_file_path,  # the directory with the `.mtx` file
                        var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
                        cache=True,  # write a cache file for faster subsequent reading
                        )

adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`

# Show those genes that yield the highest fraction of counts in each single cell, across all cells.
# sc.pl.highest_expr_genes(adata, n_top=20)

print(adata)
# Filter cells that have less than 200 genes expressed (likely dead)
sc.pp.filter_cells(adata, min_genes=200)
# Filter gene columns that have less than 3 cells with expression (not informative)
sc.pp.filter_genes(adata, min_cells=3)
print(adata)

# QC Check by mitochonrdial genes to total gene ratio
# Higher ratio means cell likely broken and leaked RNA
# while mitochondria is preserved.
# Lower ratio means likely cell clumps in single read
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)

# For plotting the number of genes expressed in the count matrix
# the total counts per cell
# the percentage of counts in mitochondrial genes
# sc.pl.violin(
#     adata,
#     ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
#     jitter=0.4,
#     multi_panel=True,
# )

# Remove cells that have too many mitochondrial genes expressed or
# too many total counts

# For looking at scatter plot to see where our cutoff should be
# sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt")
# sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts")


