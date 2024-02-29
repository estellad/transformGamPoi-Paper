alpha = 0.05

logp1_fnc <- function(UMI, sf, alpha){
  transformGamPoi::shifted_log_transform(UMI, pseudo_count = 1, size_factors = sf)
}

logp1_size_normed_fnc <- function(UMI, sf, alpha){
  x <- transformGamPoi::shifted_log_transform(UMI, pseudo_count = 1, size_factors = sf)
  # This differs from the code in Suppl. Fig. 12 of the
  # _Depth normalization for single-cell genomics count data_ paper
  # because (I think) they confused the row and colSums function
  # But as we want to normalize the depth per cell, we need to works with the
  # *colSums*.
  pf <- MatrixGenerics::colSums2(x)
  pf <- pf / mean(pf)
  sweep(x, MARGIN = 2, STATS = pf, FUN = "/")
}

logp_cpm_fnc  <- function(UMI, sf, alpha){
  colsums <- MatrixGenerics::colSums2(UMI)
  log1p(t(t(UMI) / colsums * 1e6))
}

logp_alpha_fnc  <- function(UMI, sf, alpha){
  transformGamPoi::shifted_log_transform(UMI, overdispersion = alpha, size_factors = sf, on_disk = FALSE)
}

acosh_fnc  <- function(UMI, sf, alpha){
  transformGamPoi::acosh_transform(UMI, overdispersion = alpha, size_factors = sf, on_disk = FALSE)
}

pearson_fnc  <- function(UMI, sf, alpha){
  transformGamPoi::residual_transform(UMI, overdispersion = alpha, residual_type = "pearson", size_factor = sf, on_disk = FALSE)
}

pearson_clip_fnc  <- function(UMI, sf, alpha){
  transformGamPoi::residual_transform(UMI, overdispersion = alpha, residual_type = "pearson", size_factor = sf, clipping = TRUE, on_disk = FALSE)
}

pearson_clip_zscore_fnc <- function(UMI, sf, alpha){
  res <- pearson_clip_fnc(UMI, sf, alpha)
  z_score_rows(res)
}

pearson_clip_hvg_zscore_fnc <- function(UMI, sf, alpha){
  res <- pearson_clip_fnc(UMI, sf, alpha)
  z_score_rows(select_hvg(res))
}

pearson_clip_hvg_fnc <- function(UMI, sf, alpha){
  res <- pearson_clip_fnc(UMI, sf, alpha)
  select_hvg(res)
}

pearson_analytic_fnc  <- function(UMI, sf, alpha){
  transformGamPoi::residual_transform(UMI, overdispersion = alpha, residual_type = "analytic_pearson", size_factor = sf, clipping = TRUE, on_disk = FALSE)
}

sctransform_fnc <- function(UMI, sf, alpha){
  UMI <- as(UMI, "dgCMatrix")
  colnames(UMI) <- paste0("Cell_", seq_len(ncol(UMI)))
  rownames(UMI) <- paste0("Gene_", seq_len(nrow(UMI)))
  res <- sctransform::vst(UMI, vst.flavor = "v2")$y
  stopifnot(ncol(res) == ncol(UMI))
  res
}

rand_quantile_fnc  <- function(UMI, sf, alpha){
  rand_quantile <- transformGamPoi::residual_transform(UMI, overdispersion = alpha, residual_type = "randomized_quantile", size_factor = sf, on_disk = FALSE)
  rand_quantile[matrixStats::rowAnyMissings(rand_quantile), ] <- 0
  rand_quantile
}

all_transformations <- list(logp1 = logp1_fnc,
                            acosh = acosh_fnc,
                            logp_alpha = logp_alpha_fnc,
                            logp_cpm = logp_cpm_fnc,

                            logp1_size_normed = logp1_size_normed_fnc,
                            # logp1_hvg = logp1_hvg_fnc,
                            # logp1_zscore = logp1_zscore_fnc,
                            # logp1_hvg_zscore = logp1_hvg_zscore_fnc,


                            pearson = pearson_fnc,
                            sctransform = sctransform_fnc,
                            pearson_analytic = pearson_analytic_fnc,
                            rand_quantile = rand_quantile_fnc,
                            pearson_clip = pearson_clip_fnc,

                            # pearson_clip_hvg = pearson_clip_hvg_fnc,
                            # pearson_clip_zscore = pearson_clip_zscore_fnc,
                            # pearson_clip_hvg_zscore = pearson_clip_hvg_zscore_fnc,


                            # sanity_map = sanity_map_fnc,
                            # sanity_dists = sanity_dists_fnc,
                            # dino = dino_fnc,
                            # normalisr_norm = normalisr_norm_fnc,


                            # glmpca = glmpca_fnc,
                            # newwave = newwave_fnc,


                            raw_counts = raw_counts_fnc,
                            scaled_raw_counts = scaled_raw_counts_fnc)
