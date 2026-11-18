# Small synthetic metadata: 5 clusters x 40 cells; cluster 5 is a low-count
# (debris-like) cluster.
toy_meta <- function() {
  set.seed(1)
  m <- data.frame(cl = rep(1:5, each = 40),
                  nCount_RNA = rpois(200, 3000),
                  percent.mt = runif(200, 0, 5))
  m$nCount_RNA[m$cl == 5] <- 200
  m
}

toy_metrics <- function() {
  list(nCount_RNA = list(direction = "low_bad"),
       percent.mt = list(direction = "high_bad", scale = "raw"))
}

toy_seurat <- function() {
  set.seed(1)
  counts <- matrix(rpois(50 * 60, 5), nrow = 50,
                   dimnames = list(paste0("g", 1:50), paste0("c", 1:60)))
  obj <- Seurat::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  obj$cl <- rep(c("1", "2", "10"), each = 20)
  obj
}

# pbmc_small ships with SeuratObject (80 cells, 230 genes). It has no UMAP, so
# reuse the t-SNE embedding under the name "umap" for functions that need it.
pbmc_toy <- function() {
  obj <- SeuratObject::pbmc_small
  emb <- SeuratObject::Embeddings(obj, "tsne")
  colnames(emb) <- c("UMAP_1", "UMAP_2")
  obj[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = emb, key = "UMAP_", assay = "RNA")
  obj
}

pbmc_genes <- function(n = 4) head(rownames(SeuratObject::pbmc_small), n)
