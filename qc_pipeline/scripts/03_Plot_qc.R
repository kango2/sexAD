inputs <- commandArgs(trailingOnly = TRUE)
workdir <- inputs[1]
plot_metrics <- character(0)
plot_bar <- character(0)
plot_heatmap <- character(0)
plot_scatter <- character(0)
if (length(inputs) >= 2 && nzchar(inputs[2])) plot_metrics <- strsplit(inputs[2], " ")[[1]]
if (length(inputs) >= 3 && nzchar(inputs[3])) plot_bar <- strsplit(inputs[3], " ")[[1]]
if (length(inputs) >= 4 && nzchar(inputs[4])) plot_heatmap <- strsplit(inputs[4], " ")[[1]]
if (length(inputs) >= 5 && nzchar(inputs[5])) plot_scatter <- strsplit(inputs[5], "\\|", fixed = FALSE)[[1]]

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

plot_dir <- paste0(workdir, "/04_Plots")
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
    if (sum(keep) < 2) {
      warning(paste("Skipping scatter trend line for", label, "- not enough finite points"))
      next
    }
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
    fit <- lm(yvals[keep] ~ xvals[keep])
    if (all(is.finite(coef(fit)))) {
      abline(fit, col = "firebrick3", lwd = 2)
    }
    dev.off()
  }
}

message(paste("Wrote plots to", plot_dir))
