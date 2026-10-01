## =============================================================================
##  Utility functions used 


## ---- Small utilities -----------------------------------------------------

.have <- function(p) requireNamespace(p, quietly = TRUE)


`%||%` <- function(a, b) if (is.null(a)) b else a

## Target model size.  
target_size <- function(n) floor(n / log(n))

## Split cluster ids into two halves (proxy estimation vs screening).
split_clusters <- function(ID, frac = 0.5) {
  u <- unique(ID); k <- max(1L, floor(frac * length(u)))
  s1 <- sample(u, k)
  list(prox = which(ID %in% s1), screen = which(!(ID %in% s1)))
}


## --- signal magnitude helper -------------------------------------------------
##  B = 1 and 0.5 are constant-order signals; "n^-0.2" and "n^-0.35" put the
##  design inside the minimal-signal regime of assumption (B3), which is where
##  the theory (and the validity of the I0-P proxy) actually lives.
signal_B <- function(spec, n) {
  if (is.numeric(spec)) return(spec)
  switch(spec,
         "1"      = 1,
         "0.5"    = 0.5,
         "n^-0.2" = n^(-0.2),
         "n^-0.35" = n^(-0.35),
         stop("unknown signal spec: ", spec))
}

## ==============================================================
##  Here the target column is DERIVED from the true active set, 
##  so the intended contrast holds in every signal scenario:
##      target = "active"   -> the first active screened column
##      target = "inactive" -> the first column that is not active
##  Every call returns, in $meta, which column was hit and whether it was
##  active, so the results tables can state this per configuration rather than
##  relying on prose.  
## =============================================================================

pick_col <- function(target, active, p) {
  switch(target,
         "active"   = active[1],
         "inactive" = setdiff(1:p, active)[1],
         "col2"     = 2L,
         "col6"     = 6L,
         stop("unknown target: ", target))
}


## -----------------------------------------------------------------------------
## Aggregation across replications, in the format used for the paper's tables:
## median TPR with EmpSSP in parentheses, plus the new columns.
## -----------------------------------------------------------------------------
summarise_runs <- function(df, by = c("method")) {
  for (col in by) df[[col]] <- factor(df[[col]], levels = unique(df[[col]]))
  spl <- split(df, df[by], drop = TRUE)
  do.call(rbind, lapply(levels(df[[by[1]]]), function(k) {
    x <- spl[[k]]
    data.frame(
      key       = k,
      n_rep     = nrow(x),
      TPR_med   = stats::median(x$TPR,   na.rm = TRUE),
      EmpSSP    = mean(x$EmpSSP,         na.rm = TRUE),
      MinMS_med = stats::median(x$MinMS, na.rm = TRUE),
      MinMS_mean= mean(x$MinMS,          na.rm = TRUE),
      FP_med    = stats::median(x$FP,    na.rm = TRUE),
      FDP_med   = stats::median(x$FDP,   na.rm = TRUE),
      rank_med  = stats::median(x$rank_med, na.rm = TRUE),
      rank_q90  = stats::median(x$rank_q90, na.rm = TRUE),
      time_mean = mean(x$time,           na.rm = TRUE),
      errSigma_med = stats::median(x$err_Sigma, na.rm = TRUE),
      errV_med     = stats::median(x$err_V,     na.rm = TRUE),
      nonconv_mean = mean(x$n_nonconv,   na.rm = TRUE),
      stringsAsFactors = FALSE)
  }))
}




################################################
##  boxplot_enhanced For personalized boxplot figure generations 

boxplot_enhanced <- function(data, group_var = NULL, facet_var = NULL,
                             show = c("box", "violin", "outliers",
                                      "summary_point", "summary_text", "jitter"),
                             summary_fun = mean,
                             ref_lines   = NULL,
                             grid_lines = TRUE,
                             horiz = FALSE,
                             data_lim=NULL,
                             legend_pos   = c("topleft","topright","bottomleft","bottomright",
                                              "top","bottom","left","right","center","none"),
                             box_args     = list(border = "gray40", lwd = 1),
                             violin_args  = list(border = NA),
                             summary_args = list(pch = 19, col = "black"),
                             outlier_args = list(pch = 1, col = "gray40"),
                             jitter_args  = list(pch = 16, col = rgb(0.3,0.3,0.3,0.4), cex = 0.6),
                             text_args    = list(cex = 1, col = "black", pos = 3),
                             group_palette = NULL,
                             plot_margins = c(4,5,2,2), outer_margins = c(0,0,0,0),
                             ... ) {

  merge_args <- function(defaults, user) { defaults[names(user)] <- NULL; c(defaults, user) }

  valid_show <- c("box", "violin", "summary_point", "summary_text", "outliers", "jitter")
  if(!all(show %in% valid_show)) stop("Invalid options in 'show' argument")

  legend_pos <- match.arg(legend_pos)

  if (is.matrix(data)) data <- as.data.frame(data)
  stopifnot(is.data.frame(data))
  num_cols <- names(data)[sapply(data, is.numeric)]
  if (!length(num_cols)) stop("No numeric columns found")

  df <- data.frame(
    Value = unlist(data[num_cols], use.names = FALSE),
    Method = rep(num_cols, each = nrow(data))
  )
  df$Group <- if (!is.null(group_var)) rep(data[[group_var]], times = length(num_cols)) else "All"
  df$Facet <- if (!is.null(facet_var)) rep(data[[facet_var]], times = length(num_cols)) else NULL

  methods <- unique(df$Method)
  n_methods <- length(methods)

  groups  <- unique(df$Group)
  n_groups  <- length(groups)

  facets  <- if (!is.null(df$Facet)) unique(df$Facet) else NULL
  n_facets  <- if (!is.null(facets)) length(facets) else 1

  ## ---- assign colors for groups ----
  if (!is.null(group_palette)) {
    if(length(group_palette) < n_groups) stop("group_palette has fewer colors than groups")
    group_cols <- setNames(group_palette[1:n_groups], groups)
  } else {
    default_palette <- c("#A6CEE3","#B2DF8A","#FB9A99","#FDBF6F","#CAB2D6","#FFFF99","#1F78B4","#33A02C")
    group_cols <- setNames(default_palette[1:n_groups], groups)
  }
  legend_cols <- sapply(group_cols, function(col) {
    rgb_val <- grDevices::col2rgb(col)/255
    rgb(rgb_val[1], rgb_val[2], rgb_val[3], alpha = 0.9)
  })

  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mfrow = c(1, n_facets), mar = plot_margins, oma = outer_margins)

  max_dodge <- 0.8
  dodge_width <- max_dodge / max(1, n_groups)
  box_width <- 0.7 * dodge_width
  violin_width <- 0.85 * dodge_width

  min_width <- 0.1
  box_width <- max(box_width, min_width)
  violin_width <- max(violin_width, min_width)

  for (f in seq_len(n_facets)) {
    sub_df <- if (!is.null(facets)) df[df$Facet == facets[f], , drop = FALSE] else df
    main_main <- if (!is.null(facets)) facets[f] else ""

    data_range <- if(is.null(data_lim)) range(sub_df$Value, na.rm = TRUE) else data_lim
    x_range <- if(!horiz) c(0.9, n_methods+0.1) else data_range
    y_range <- if(!horiz) data_range else c(0.9, n_methods+0.1)

    plot(1, type="n", xlim=x_range, ylim=y_range, main = main_main,
         xaxt="n", yaxt="n", ...)
    box()

    if(grid_lines){
      if(!horiz) abline(h = pretty(y_range), col="gray80", lty=2)
      else abline(v = pretty(x_range), col="gray80", lty=2)
    }

    if(!horiz){
      # axis(1, at = 1:n_methods, labels = methods, las=1, ...)
      axis(1, at = 1:n_methods, labels = FALSE, las=1, ...)
      text(
        x = 1:n_methods, 
        y = par("usr")[3] - 5,      # Position slightly below the bottom of the plot
        labels = methods, las =1,
        srt = 45,                      # Angle of rotation (e.g., 45 degrees)
        adj = 1,                       # Right-justify text so the end aligns with the tick
        xpd = T,                    # Allow text to be drawn outside the plot region
        ...                      # Adjust font size if needed
      )
      axis(2, ...)
    } else {
      axis(2, at = 1:n_methods, labels = methods, las=1, ...)
      axis(1, ...)
    }

    if(!is.null(ref_lines)){
      if(!horiz) abline(h=ref_lines, lty=2, col="gray50") else abline(v=ref_lines, lty=2, col="gray50")
    }

    for (i in seq_along(methods)) {
      for (g in seq_along(groups)) {

        vals <- sub_df$Value[sub_df$Method==methods[i] & sub_df$Group==groups[g]]
        vals <- vals[!is.na(vals)]
        if(length(vals)<1) next

        x_pos <- i + (g-(n_groups+1)/2)*dodge_width
        grp_col <- group_cols[groups[g]]

        if (length(vals) >= 1) {
          qs <- quantile(vals, c(0.25, 0.5, 0.75), na.rm = TRUE)
          iqr <- qs[3] - qs[1]
          lo  <- qs[1] - 1.5 * iqr
          hi  <- qs[3] + 1.5 * iqr
          in_range <- vals[vals >= lo & vals <= hi]
          w <- if (length(in_range))
            range(in_range) else range(vals, na.rm = TRUE)
        }

        if ("violin" %in% show && length(vals) > 1) {
          dens <- density(vals, from = min(vals), to = max(vals), cut = 0)
          dens$y <- dens$y / max(dens$y) * violin_width / 2
          col_violin <- rgb(t(col2rgb(grp_col)/255), alpha = 0.6)
          violin_args_final <- merge_args(violin_args, list(col = col_violin))

          if(!horiz){
            x_poly <- c(x_pos + dens$y, rev(x_pos - dens$y))
            y_poly <- c(dens$x, rev(dens$x))
            do.call(polygon, c(list(x = x_poly, y = y_poly), violin_args_final))
          } else {
            x_poly <- c(dens$x, rev(dens$x))
            y_poly <- c(x_pos + dens$y, rev(x_pos - dens$y))
            do.call(polygon, c(list(x = x_poly, y = y_poly), violin_args_final))
          }
        }

        if ("box" %in% show && length(vals)>1) {
          col_box <- rgb(t(col2rgb(grp_col)/255), alpha=1)
          box_args_final <- merge_args(box_args, list(col=col_box))
          if(!horiz){
            do.call(rect, c(list(xleft=x_pos-box_width/2, ybottom=qs[1],
                                 xright=x_pos+box_width/2, ytop=qs[3]), box_args_final))
            segments(x_pos, w[1], x_pos, qs[1])
            segments(x_pos, qs[3], x_pos, w[2])
            segments(x_pos - box_width/2, qs[2], x_pos + box_width/2, qs[2], lwd = 2)
          } else {
            do.call(rect, c(list(xleft=qs[1], ybottom=x_pos-box_width/2,
                                 xright=qs[3], ytop=x_pos+box_width/2), box_args_final))
            segments(w[1], x_pos, qs[1], x_pos)
            segments(qs[3], x_pos, w[2], x_pos)
            segments(qs[2], x_pos - box_width/2, qs[2], x_pos + box_width/2, lwd = 2)
          }
        }

        if("summary_point" %in% show){
          val <- summary_fun(vals, na.rm = TRUE)
          summary_args_final <- merge_args(summary_args, list())
          if(!horiz) do.call(points, c(list(x=x_pos, y=val), summary_args_final))
          else do.call(points, c(list(x=val, y=x_pos), summary_args_final))
          if("summary_text" %in% show){
            text_args_final <- merge_args(text_args, list())
            lab <- round(val, 2)
            if(!horiz) do.call(text, c(list(x = x_pos, y = val, labels = lab), text_args_final))
            else do.call(text, c(list(x = val, y = x_pos, labels = lab), text_args_final))
          }
        }

        if("outliers" %in% show){
          out <- vals[vals<lo | vals>hi]
          if(length(out)){
            outlier_args_final <- merge_args(outlier_args, list())
            if(!horiz) do.call(points, c(list(x=rep(x_pos,length(out)), y=out), outlier_args_final))
            else do.call(points, c(list(x=out, y=rep(x_pos,length(out))), outlier_args_final))
          }
        }

        if("jitter" %in% show){
          jitter_args_final <- merge_args(jitter_args, list())
          if(!horiz) do.call(points,
                             c(list(x = jitter(rep(x_pos,length(vals)), box_width/4),
                                    y = vals), jitter_args_final))
          else do.call(points, c(list(x = vals,
                                      y = jitter(rep(x_pos,length(vals)),box_width/4)),
                                 jitter_args_final))
        }
      }
    }

    if(n_groups>1 && legend_pos != "none")
      legend(legend_pos, legend = groups, fill = legend_cols)
  }
}




