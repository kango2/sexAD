inputs <- commandArgs(trailingOnly = TRUE)
workdir <- inputs[1]
metadata_path <- paste0(workdir, "/metadata_input.csv")
metadata_vars <- character(0)
metadata_plot_metrics <- character(0)
plot_metrics <- character(0)
plot_bar <- character(0)
plot_heatmap <- character(0)
plot_scatter <- character(0)
if (length(inputs) >= 2 && nzchar(inputs[2])) metadata_vars <- strsplit(inputs[2], " ")[[1]]
if (length(inputs) >= 3 && nzchar(inputs[3])) metadata_plot_metrics <- strsplit(inputs[3], " ")[[1]]
if (length(inputs) >= 4 && nzchar(inputs[4])) plot_metrics <- strsplit(inputs[4], " ")[[1]]
if (length(inputs) >= 5 && nzchar(inputs[5])) plot_bar <- strsplit(inputs[5], " ")[[1]]
if (length(inputs) >= 6 && nzchar(inputs[6])) plot_heatmap <- strsplit(inputs[6], " ")[[1]]
if (length(inputs) >= 7 && nzchar(inputs[7])) plot_scatter <- strsplit(inputs[7], "\\|", fixed = FALSE)[[1]]

qc_path <- paste0(workdir, "/02_Metadata/combined_qc.tsv")
if (!file.exists(qc_path)) {
  stop(paste("QC table not found:", qc_path))
}

qc <- read.delim(
  qc_path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"Sample" %in% colnames(qc)) {
  stop("combined_qc.tsv must contain a Sample column")
}

meta <- NULL
if (file.exists(metadata_path)) {
  meta <- read.csv(metadata_path, header = TRUE, stringsAsFactors = FALSE)
  if (!"file" %in% colnames(meta)) {
    stop("metadata_input.csv must contain a file column")
  }
  meta <- meta[!duplicated(meta$file), , drop = FALSE]
  qc <- merge(qc, meta, by.x = "Sample", by.y = "file", all.x = TRUE)
} else {
  warning(paste("Metadata file not found, skipping metadata-aware plots:", metadata_path))
}

plot_dir <- paste0(workdir, "/03_Plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

numeric_cols <- setdiff(names(qc), "Sample")
numeric_cols <- numeric_cols[sapply(qc[numeric_cols], function(x) suppressWarnings(!all(is.na(as.numeric(x)))))]

if (length(numeric_cols) == 0) {
  stop("No numeric columns available for plotting")
}

for (nm in numeric_cols) {
  qc[[nm]] <- suppressWarnings(as.numeric(qc[[nm]]))
}

if (length(plot_metrics) > 0) {
  missing_metrics <- setdiff(plot_metrics, names(qc))
  if (length(missing_metrics) > 0) {
    warning(paste("Skipping missing plot metrics:", paste(missing_metrics, collapse = ", ")))
  }
  plot_metrics <- intersect(plot_metrics, numeric_cols)
} else {
  plot_metrics <- head(numeric_cols, 4)
  missing_metrics <- character(0)
}

plot_bar <- intersect(unique(c(plot_bar, plot_metrics)), numeric_cols)
plot_heatmap <- intersect(unique(c(plot_heatmap, plot_metrics)), numeric_cols)

scatter_specs <- list()
if (length(plot_scatter) > 0) {
  for (spec in plot_scatter) {
    parts <- strsplit(spec, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 3) {
      warning(paste("Skipping malformed scatter spec:", spec))
      next
    }
    scatter_specs[[length(scatter_specs) + 1]] <- list(
      x = parts[1],
      y = parts[2],
      label = parts[3]
    )
  }
}

run_summary <- c(
  paste0("QC plot summary for: ", workdir),
  paste0("Samples: ", nrow(qc)),
  paste0("Metadata variables requested: ", paste(metadata_vars, collapse = ", ")),
  paste0("Metadata plot metrics: ", paste(metadata_plot_metrics, collapse = ", ")),
  paste0("Numeric metrics available: ", length(numeric_cols)),
  paste0("Bar plots: ", length(plot_bar)),
  paste0("Heatmap metrics: ", length(plot_heatmap)),
  paste0("Scatter panels: ", length(scatter_specs)),
  paste0("Plotted metrics: ", paste(plot_metrics, collapse = ", "))
)
if (length(missing_metrics) > 0) {
  run_summary <- c(run_summary, paste0("Missing plot metrics: ", paste(missing_metrics, collapse = ", ")))
}
writeLines(run_summary, paste0(plot_dir, "/run_summary.txt"))

summary_pdf <- paste0(plot_dir, "/qc_summary_plots.pdf")
pdf(summary_pdf, width = 11, height = 8.5)
op <- par(no.readonly = TRUE)
on.exit({
  par(op)
  dev.off()
}, add = FALSE)

layout_n <- max(1, min(4, length(plot_bar)))
par(mfrow = c(ceiling(layout_n / 2), min(2, layout_n)), mar = c(8, 4, 3, 1))

for (nm in plot_bar) {
  values <- qc[[nm]]
  bar_cols <- ifelse(is.na(values), "grey80", "steelblue")
  names_arg <- qc$Sample
  if (length(values) > 40) {
    ord <- order(values, decreasing = TRUE, na.last = NA)
    values <- values[ord]
    names_arg <- names_arg[ord]
    bar_cols <- bar_cols[ord]
  }
  barplot(
    values,
    names.arg = names_arg,
    las = 2,
    col = bar_cols,
    main = nm,
    ylab = nm,
    cex.names = 0.6
  )
}

make_metric_plot <- function(metric, file_base) {
  png(paste0(plot_dir, "/", file_base, ".png"), width = 1800, height = 900, res = 150)
  ord <- order(qc[[metric]], decreasing = TRUE, na.last = NA)
  values <- qc[[metric]][ord]
  samples <- qc$Sample[ord]
  barplot(
    values,
    names.arg = samples,
    las = 2,
    col = "steelblue4",
    main = metric,
    ylab = metric,
    cex.names = 0.5,
    cex.axis = 0.8
  )
  dev.off()
}

for (nm in plot_bar) {
  safe_name <- gsub("[^A-Za-z0-9_]+", "_", nm)
  make_metric_plot(nm, paste0("metric_", safe_name))
}

if (length(plot_heatmap) >= 2) {
  heat_cols <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
  heat_mat <- as.matrix(qc[plot_heatmap])
  rownames(heat_mat) <- qc$Sample
  heat_mat <- scale(heat_mat)
  heat_mat[is.na(heat_mat)] <- 0
  png(paste0(plot_dir, "/qc_metric_heatmap.png"), width = 1600, height = 1200, res = 160)
  image(
    1:ncol(heat_mat),
    1:nrow(heat_mat),
    t(heat_mat[nrow(heat_mat):1, ]),
    axes = FALSE,
    col = heat_cols,
    main = "QC Metric Heatmap"
  )
  axis(1, at = 1:ncol(heat_mat), labels = colnames(heat_mat), las = 2, cex.axis = 0.6)
  axis(2, at = 1:nrow(heat_mat), labels = rev(rownames(heat_mat)), las = 2, cex.axis = 0.4)
  box()
  dev.off()
}

plot_group_boxplot <- function(metric, group_col, file_base) {
  if (!group_col %in% names(qc)) {
    return(invisible(NULL))
  }
  groups <- qc[[group_col]]
  values <- qc[[metric]]
  keep <- !is.na(groups) & !is.na(values)
  if (sum(keep) < 2 || length(unique(groups[keep])) < 2) {
    return(invisible(NULL))
  }
  png(paste0(plot_dir, "/", file_base, ".png"), width = 1600, height = 1200, res = 160)
  boxplot(
    values[keep] ~ as.factor(groups[keep]),
    las = 2,
    col = "steelblue2",
    main = paste(metric, "by", group_col),
    xlab = group_col,
    ylab = metric
  )
  dev.off()
}

if (length(metadata_vars) > 0) {
  metrics_for_group_plots <- intersect(metadata_plot_metrics, numeric_cols)
  for (meta_var in metadata_vars) {
    if (!meta_var %in% names(qc)) {
      warning(paste("Skipping missing metadata variable:", meta_var))
      next
    }
    if (length(metrics_for_group_plots) == 0) {
      warning(paste("No metadata plot metrics available for", meta_var))
      next
    }
    for (nm in metrics_for_group_plots) {
      safe_name <- gsub("[^A-Za-z0-9_]+", "_", nm)
      plot_group_boxplot(nm, meta_var, paste0(meta_var, "_boxplot_", safe_name))
    }
  }
}

if (length(scatter_specs) > 0) {
  for (spec in scatter_specs) {
    x <- spec$x
    y <- spec$y
    label <- spec$label
    if (!all(c(x, y) %in% names(qc))) {
      warning(paste("Skipping scatter plot for missing columns:", x, y))
      next
    }
    xvals <- qc[[x]]
    yvals <- qc[[y]]
    keep <- is.finite(xvals) & is.finite(yvals)
    png(paste0(plot_dir, "/scatter_", gsub("[^A-Za-z0-9_]+", "_", label), ".png"), width = 1400, height = 1200, res = 160)
    plot(
      xvals,
      yvals,
      xlab = x,
      ylab = y,
      main = label,
      pch = 16,
      col = "steelblue4"
    )
    if (sum(keep) >= 5) {
      fit <- lm(yvals[keep] ~ xvals[keep])
      if (all(is.finite(coef(fit)))) {
        abline(fit, col = "firebrick3", lwd = 2)
      }
    } else {
      warning(paste("Skipping scatter trend line for", label, "- fewer than 5 finite points"))
    }
    dev.off()
  }
}

message(paste("Wrote plots to", plot_dir))
