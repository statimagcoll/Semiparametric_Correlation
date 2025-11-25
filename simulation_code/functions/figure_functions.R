### Nicer and more organized plotting functions ###
library(dplyr)
library(ggplot2)
library(latex2exp)
library(tidyr)
library(patchwork)
library(grid)
library(colorspace)


# Fig 1 panels
# A - parameters/bias by model sample size
# B - Coverage of typical Fisher CI for the parameters

fig_1_bias = function(r, dt, dtr2, r0, r1, r2, ns){
  #browser()
  #dt = get(paste0("allResults_",r))
  
  pdt = dt %>% # confidence interval bands for rho_1
    group_by(n) %>%
    summarize(plower=mean(lCI_p_inst),
              pupper=mean(uCI_p_inst),
              p = mean(p_inst))
  
  #dtr2 = get(paste0("allResults_rho_2_",r))
  pdt_r2 =  dtr2 %>% # confidence interval bands for rho_2
    group_by(n) %>%
    summarize(plower=mean(lCI_p_hold),
              pupper=mean(uCI_p_hold),
              p = mean(p_hold))
  
  ggplot(data=data.frame(ns, y=r1), mapping = aes(x=ns, y=y))+
 
    geom_ribbon(data=pdt_r2, mapping=aes(x=n, ymin=plower, ymax=pupper, y=NULL, fill="8 CIr2"), alpha =0.2) + # confidence interval for pearson (rho_2)
    
    geom_ribbon(data=pdt, mapping=aes(x=n, ymin=plower, ymax=pupper, y=NULL, fill="7 CIr1"), alpha =0.2) + # confidence interval for pearson(rho_1)
  
    
     # MAPA parameter
    geom_hline(aes(yintercept = r0, color="2 rho_0"), size = 1.5) + # the correlation from OLS in the population
    # MPA parameter
    geom_hline(aes(yintercept = r/100, color="1 rho_true"), size = 1.5) +  # rho_true 
    
    # Pearson estimator MPA/MAPA
    geom_line(data = pdt,mapping= aes(x = n, y=p), color="#D55E00", lty =1, size = 1.5) +
    geom_point(data = pdt,mapping= aes(x = n, y=p, color="5 mean_p_inst"), size = 3) +
    
    # TPA parameter
    geom_hline(aes(yintercept = r2[length(r2)], color="4 rho_2"), size = 1.5) +
    # TPA estimates
    stat_summary(data = dtr2,mapping= aes(x = n, y=p_hold), color="#009E73",fun = "mean", size=1.5, geom="line", lty = 1) + 
    geom_point(data = pdt_r2,mapping= aes(x = n, y=p, color="6 mean_p_model"), size = 4, pch = 18) +
    
    # axis labels
    scale_x_continuous(name = "Model Sample Size") +
    scale_y_continuous(name="Correlation") +
    
    # legend
    scale_color_manual(name = NULL, 
                       values = c("5 mean_p_inst" = "#D55E00", 
                                  "2 rho_0" = "black", 
                                  "1 rho_true" = "#CC79A7",
                                  "4 rho_2"="#F0E442",
                                  "6 mean_p_model"="#009E73"), 
                       labels = c("5 mean_p_inst" = TeX("$\\hat{\\rho}_{P, target}$"), 
                                  "2 rho_0" = TeX("$\\rho_{target}$"), 
                                  "1 rho_true" = TeX("$\\rho_{true}$"),
                                  "4 rho_2" = TeX("$\\rho_{model}$"),
                                  "6 mean_p_model" = TeX("$\\hat{\\rho}_{P, model}$"))) +
    scale_fill_manual(name=NULL,
                      values=c("7 CIr1"="#D55E00", "8 CIr2"="#009E73"),
                      labels=c("7 CIr1"=TeX("Mean 95% CI Bounds"),
                               "8 CIr2"=TeX("Mean 95% CI Bounds"))) +
    
    # theming
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          legend.box = "horizontal",legend.spacing = unit(0,'pt'),
          legend.margin = margin(t=0,b=0,unit='pt'), legend.position = "bottom") +
    
    ggtitle("Parameters and Pearson Estimate")
  
  
}


fig_1_coverage = function(r, dt, dtr2, r0, r1, r2, nsim = 1000){
  #browser()
  cptrue = get_cp(dt, "p", r/100) %>% filter(est=="p" & ci=="inst") %>% dplyr::select(n, coverage_probability)
  cp0 = get_cp(dt, "p", r0)  %>% filter(est=="p" & ci=="inst") %>% dplyr::select(coverage_probability)
  cp1 = get_cp(dt, "p", r1)  %>% filter(est=="p" & ci=="inst") %>% dplyr::select(coverage_probability)
  cp2 = get_cp(dtr2, "p", r2)  %>% filter(est=="p" & ci=="hold") %>% dplyr::select(coverage_probability)
  
  
  dt = cbind(cptrue, cp0, cp1, cp2)
  colnames(dt) = c("n", "rtrue", "r0", "r1", "r2")
  
  
  ggplot(data=dt, mapping = aes(x=n, y=r2))+
    # expected range given nsim
    geom_ribbon(mapping=aes(ymin = qbinom(0.025,nsim, 0.95)/nsim, ymax=qbinom(0.975,nsim, 0.95)/nsim), fill = "grey92") +
    # nominal coverage
    geom_hline(aes(yintercept = 0.95), col="black", lty=3, size = 1.5) +
    
    # rho_2 (TPA) coverage
    geom_point(aes(color="rho_2"), pch=18, size = 4) +
    geom_line(aes(color="rho_2"), size = 1.5) +
    
    
    # rho_0 (MAPA) coverage
    geom_point(aes(x=n, y=r0, color="rho_0"), size = 3) +
    geom_line(aes(x=n, y=r0, color="rho_0"), size = 1.5) +
    
    # axis labels
    scale_x_continuous(name = "Model Sample Size") +
    scale_y_continuous(name="Coverage Probability") +
    
    # legend
    scale_color_manual(name = "Parameter", 
                       values = c("rho_1" = "#D55E00",
                                  "rho_0" = "#D55E00", 
                                  "1rho_true" = "#CC79A7",
                                  "rho_2"="#009E73"), 
                       labels = c("rho_1" = TeX("$\\rho_1$"), 
                                  "rho_0" = TeX("$\\rho_0$"), 
                                  "1rho_true" = TeX("$\\rho_{\\xi}$"),
                                  "rho_2" = TeX("$\\rho_1$"))) +
    
    
    # theming
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), legend.position = "bottom") +
    
    ggtitle("Coverage Probability of Pearson CI")
}


# Function to extract the legend from a ggplot
get_legend <- function(p, legend_text_size = 24) {
  p <- p + guides(color = guide_legend(direction = "horizontal")) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = legend_text_size),
          legend.title = element_text(size = legend_text_size))
  tmp <- ggplotGrob(p)
  legend <- tmp$grobs[which(sapply(tmp$grobs, function(x) x$name) == "guide-box")]
  if (length(legend) > 0) return(legend[[1]]) else return(NULL)
}


## Compare bias from multiple simulation settings (e.g., different error dists)
fig_compare_bias_mult = function(dt_list, dt_labels, r0, inst = TRUE, ylims = NA, type = "pq", 
                                 title = "Bias", source_name = "source", ns = NULL, breaks_n = c(100, 5000, 10000, 20000), ...) {
  # Assign labels to names
  names(dt_list) = dt_labels
  
  # Combine datasets
  plot_data = bind_rows(lapply(names(dt_list), function(name) {
    dt = dt_list[[name]]
    if (inst) {
      dt$if_est = dt[[paste0(type, "_inst")]]
      dt$p_val  = dt$p_inst
    } else {
      dt$if_est = dt[[paste0(type, "_hold")]]
      dt$p_val  = dt$p_hold
    }
    dt %>%
      group_by(n) %>%
      summarize(
        p = mean(p_val),
        t = mean(if_est),
        .groups = "drop"
      ) %>%
      mutate(source = name)
  }))
  
  # Filter by sample sizes
  if (!is.null(ns)) {
    plot_data = plot_data %>% filter(n %in% ns)
  }
  
  # Convert to long format
  plot_data_long = tidyr::pivot_longer(
    plot_data,
    cols = c("p", "t"),
    names_to = "estimator",
    values_to = "value"
  )
  
  # Dynamically assign shapes/linetypes for all sources
  sources <- unique(plot_data_long$source)
  shape_vals <- seq(0, length.out = length(sources)) %% 6 + 1  # Use shapes 1–6 cyclically
  names(shape_vals) <- sources
  
  linetype_vals <- rep_len(1:3, length(sources))
  names(linetype_vals) <- sources
  
  # Convert source to factor in desired order
  plot_data_long$source <- factor(plot_data_long$source, levels = rev(dt_labels))
  
  # Base plot
  p = ggplot(plot_data_long, aes(x = n, y = value, color = estimator,
                                 shape = source, linetype = source)) +
    geom_hline(aes(yintercept = r0, color = "truth"), size = 1.5) +
    geom_point(size = 3) +
    geom_line(size = 1.2) +
    scale_x_continuous(name = "Model Sample Size", breaks = breaks_n) +
    scale_y_continuous(name = "Correlation") +
    scale_color_manual(
      name = NULL,
      values = c(
        "truth" = "black",
        "p" = "#D55E00",
        "t" = "#56B4E9"
      ),
      labels = c(
        "truth" = TeX("$\\rho_{target}$"),
        "p" = TeX("$\\hat{\\rho}_{P}$"),
        "t" = TeX("$\\hat{\\rho}_{OS}$")
      )
    ) +
    scale_shape_manual(
      name = source_name,
      values = shape_vals
    ) +
    scale_linetype_manual(
      name = source_name,
      values = linetype_vals
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.title = element_text(size = 22),
      legend.text = element_text(size = 22),
      legend.key.size = unit(40, "pt"),
      legend.position = "bottom"
    ) +
    ggtitle(title)
  
  if (!is.na(ylims[1])) {
    p = p + ylim(ylims)
  }
  
  p = p +
    guides(
      shape = guide_legend(title = source_name, reverse = TRUE),
      linetype = guide_legend(title = source_name, reverse = TRUE)
    )
  
  return(p)
}

## Compare coverage from multiple simulation settings (e.g., different error dists)
fig_compare_coverage_mult = function(dt_list, dt_labels, r0, cf = "inst", inst = TRUE,
                                     type = "pq", nsim = 500, title = "Coverage",
                                     source_name = "source", ns = NULL,
                                     breaks_n = c(100, 5000, 10000, 20000), ...) {
  
  names(dt_list) = dt_labels
  ci_type = ifelse(inst, "inst", "hold")
  
  # Extract all sample sizes
  all_ns = unique(unlist(lapply(dt_list, function(dt) unique(dt$n))))
  if (!is.null(ns)) {
    all_ns = all_ns[all_ns %in% ns]
  }
  
  # CI ribbon bounds for each sample size
  ribbon_df = data.frame(
    n = all_ns,
    ymin = qbinom(0.025, nsim, 0.95) / nsim,
    ymax = qbinom(0.975, nsim, 0.95) / nsim
  )
  
  # Collect coverage data
  coverage_data = bind_rows(lapply(names(dt_list), function(name) {
    dt = dt_list[[name]]
    
    cp_os = get_cp(dt, c(type), r0) %>%
      filter(est == type, ci == cf) %>%
      transmute(n, coverage = coverage_probability, estimator = "t", source = name)
    
    cp_p = get_cp(dt, c("p"), r0) %>%
      filter(est == "p", ci == ci_type) %>%
      transmute(n, coverage = coverage_probability, estimator = "p", source = name)
    
    bind_rows(cp_p, cp_os)
  }))
  
  # Filter coverage data if ns provided
  if (!is.null(ns)) {
    coverage_data = coverage_data %>% filter(n %in% ns)
  }
  
  # Dynamically generate shape/linetype values for sources
  sources <- unique(coverage_data$source)
  shape_vals <- seq(0, length.out = length(sources)) %% 6 + 1
  names(shape_vals) <- sources
  
  linetype_vals <- rep_len(1:3, length(sources))
  names(linetype_vals) <- sources
  
  # Final plot
  ggplot(coverage_data, aes(x = n, y = coverage,
                            color = estimator, shape = source, linetype = source)) +
    # Binomial coverage ribbon
    geom_ribbon(data = ribbon_df, aes(x = n, ymin = ymin, ymax = ymax),
                fill = "grey92", inherit.aes = FALSE) +
    
    # Nominal 0.95 coverage line
    geom_hline(yintercept = 0.95, color = "black", linetype = 3, size = 1.2) +
    
    # Points and lines
    geom_point(size = 3) +
    geom_line(size = 1.2) +
    
    # Axes
    scale_x_continuous(name = "Model Sample Size", breaks = breaks_n) +
    scale_y_continuous(name = "Coverage Probability", limits = c(0, 1)) +
    
    # Manual color/shape mapping
    scale_color_manual(name = NULL,
                       values = c("p" = "#D55E00", "t" = "#56B4E9"),
                       labels = c("p" = TeX("$\\hat{\\rho}_{P}$"),
                                  "t" = TeX("$\\hat{\\rho}_{OS}$"))) +
    scale_shape_manual(name = source_name, values = shape_vals) +
    scale_linetype_manual(name = source_name, values = linetype_vals) +
    
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 14),
          legend.key.size = unit(20, "pt"),
          legend.position = "bottom") +
    ggtitle(title)
}

## Compare CI width from multiple simulation settings (e.g., different error dists)
fig_compare_width_mult = function(dt_list, dt_labels, cf = "inst", inst = TRUE,
                                         type = "pq", title = "Average CI Width",
                                         source_name = "source", ns = NULL,
                                         breaks_n = c(100, 5000, 10000, 20000), ...) {
  
  names(dt_list) = dt_labels
  ci_type = ifelse(inst, "inst", "hold")
  
  # Construct CI column names
  lci_os_col = paste0("lCI_", type, "_", cf)
  uci_os_col = paste0("uCI_", type, "_", cf)
  lci_p_col  = paste0("lCI_p_", ci_type)
  uci_p_col  = paste0("uCI_p_", ci_type)
  
  # Collect width data
  width_data = bind_rows(lapply(names(dt_list), function(name) {
    dt = dt_list[[name]]
    
    # One-step width
    width_os = dt %>%
      mutate(width = !!sym(uci_os_col) - !!sym(lci_os_col)) %>%
      group_by(n) %>%
      summarize(width = mean(width, na.rm = TRUE), .groups = "drop") %>%
      mutate(width = width,
             estimator = "t", source = name)
    
    # Plug-in width
    width_p = dt %>%
      mutate(width = !!sym(uci_p_col) - !!sym(lci_p_col)) %>%
      group_by(n) %>%
      summarize(width = mean(width, na.rm = TRUE), .groups = "drop") %>%
      mutate(width = width,
             estimator = "p", source = name)
    
    bind_rows(width_p, width_os)
  }))
  
  # Filter by specified sample sizes
  if (!is.null(ns)) {
    width_data = width_data %>% filter(n %in% ns)
  }
  
  # Dynamically generate shape/linetype values for sources
  sources <- unique(width_data$source)
  shape_vals <- seq(0, length.out = length(sources)) %% 6 + 1
  names(shape_vals) <- sources
  
  linetype_vals <- rep_len(1:3, length(sources))
  names(linetype_vals) <- sources
  
  # Plot
  ggplot(width_data, aes(x = n, y = width,
                         color = estimator, shape = source, linetype = source)) +
    geom_point(size = 3) +
    geom_line(size = 1.2) +
    
    scale_x_log10(name = "Model Sample Size", breaks = breaks_n) +
    scale_y_continuous(name = TeX("Width")) +
    
    scale_color_manual(name = NULL,
                       values = c("p" = "#D55E00", "t" = "#56B4E9"),
                       labels = c("p" = TeX("$\\hat{\\rho}_{P}$"),
                                  "t" = TeX("$\\hat{\\rho}_{OS}$"))) +
    scale_shape_manual(name = source_name, values = shape_vals) +
    scale_linetype_manual(name = source_name, values = linetype_vals) +
    
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.text = element_text(size = 14),
          legend.key.size = unit(20, "pt"),
          legend.position = "bottom") +
    ggtitle(title)
}



# 3x4 panel of results (bias, coverage, width)
plot_panel_mult_full <- function(dt_suffixes, dt_labels, inst = TRUE, type = "pq", cf = "inst", nsim = 500, ylims = c(0.3, 1), ...) {
  #browser()
  # Create list of all data.frames
  dt_list = list(r10 = list(), r25 = list(), r50 = list(), r75 = list())
  i = 1
  for (s in dt_suffixes){
    dt_list$r10[[i]] = get(paste0("allResults_10_",s))
    dt_list$r25[[i]] = get(paste0("allResults_25_",s))
    dt_list$r50[[i]] = get(paste0("allResults_50_",s))
    dt_list$r75[[i]] = get(paste0("allResults_75_",s))
    i = i + 1
  }
  names(dt_list$r10) = names(dt_list$r25) = names(dt_list$r50) = names(dt_list$r75) = dt_labels
  
  titles <- list(
    r10 = "MAPA = 0.1",
    r25 = "MAPA = 0.25",
    r50 = "MAPA = 0.5",
    r75 = "MAPA = 0.75"
  )

  # Bias plots (row 1)
  pb <- list(
    fig_compare_bias_mult(dt_list = dt_list$r10, dt_labels = dt_labels, r0 = 0.1, inst = inst, type = type, title = titles$r10, ...),
    fig_compare_bias_mult(dt_list = dt_list$r25, dt_labels = dt_labels, r0 = 0.25, inst = inst, type = type, title = titles$r25, ...),
    fig_compare_bias_mult(dt_list = dt_list$r50, dt_labels = dt_labels, r0 = 0.5, inst = inst, type = type, title = titles$r50, ...),
    fig_compare_bias_mult(dt_list = dt_list$r75, dt_labels = dt_labels, r0 = 0.75, inst = inst, type = type, title = titles$r75, ...)
  )
  
  # Coverage plots (row 2)
  pc <- list(
    fig_compare_coverage_mult(dt_list = dt_list$r10, dt_labels = dt_labels, r0 = 0.1, inst = inst, type = type, title = NULL, nsim = nsim, cf = cf),
    fig_compare_coverage_mult(dt_list = dt_list$r25, dt_labels = dt_labels, r0 = 0.25, inst = inst, type = type, title = NULL, nsim = nsim, cf = cf),
    fig_compare_coverage_mult(dt_list = dt_list$r50, dt_labels = dt_labels, r0 = 0.5, inst = inst, type = type, title = NULL, nsim = nsim, cf = cf),
    fig_compare_coverage_mult(dt_list = dt_list$r75, dt_labels = dt_labels, r0 = 0.75, inst = inst, type = type, title = NULL, nsim = nsim, cf = cf)
  )
  
  # Width plots (row 3)
  pw <- list(
    fig_compare_width_mult(dt_list = dt_list$r10, dt_labels = dt_labels, inst = inst, type = type, title = NULL, breaks_n = c(100, 500, 2000, 20000)),
    fig_compare_width_mult(dt_list = dt_list$r25, dt_labels = dt_labels, inst = inst, type = type, title = NULL, breaks_n = c(100, 500, 2000, 20000)),
    fig_compare_width_mult(dt_list = dt_list$r50, dt_labels = dt_labels, inst = inst, type = type, title = NULL, breaks_n = c(100, 500, 2000, 20000)),
    fig_compare_width_mult(dt_list = dt_list$r75, dt_labels = dt_labels, inst = inst, type = type, title = NULL, breaks_n = c(100, 500, 2000, 20000))
  )
  
  # Get legend from the first bias plot
  legend <- get_legend(pb[[1]])
  
  # Combine plots into a list, flattening rows
  plots <- c(pb, pc, pw)
  
  # Apply consistent theming
  plots_nolegend <- map2(plots, seq_along(plots), function(p, i) {
    row_idx <- ceiling(i / 4)
    col_idx <- ((i - 1) %% 4) + 1
    
    is_bottom_row <- row_idx == 3
    is_left_col <- col_idx == 1
    
    p + theme(
      legend.position = "none",
      axis.title.x = if (!is_bottom_row) element_blank() else element_blank(),
      axis.title.y = if (!is_left_col) element_blank() else element_text(size = 24),
      axis.text.x = element_text(size = 22),
      axis.text.y = element_text(size = 22),
      plot.title = element_text(size = 24),
      plot.margin = margin(t = 5, r = 25, b = 5, l = 5)  
    )
  })
  
  # Plot grid: 3 rows (Bias, Coverage, Width) × 4 columns (rho)
  # Combine panels into 3 row blocks (each 4-column)
  row1 <- plot_grid(plotlist = plots_nolegend[1:4], ncol = 4, align = "h")
  row2 <- plot_grid(plotlist = plots_nolegend[5:8], ncol = 4, align = "h")
  row3 <- plot_grid(plotlist = plots_nolegend[9:12], ncol = 4, align = "h")
  
  # Stack row blocks and label the rows A / B / C
  main_panel <- plot_grid(
    row1, row2, row3,
    ncol = 1,
    labels = c("A", "B", "C"),
    label_x = 0,          # flush left
    label_y = 1,          # start at top of row
    hjust = -0.2,         # move label slightly outside plot
    vjust = 0.8,
    label_size = 24
  )
  
  
  
  # Add shared legend below
  
  xaxis_label <- ggplot() + 
    theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = "Model Sample Size", size = 8, hjust = 0.5)
  
  # Combine grid + legend + x-axis label
  final_plot <- (main_panel / xaxis_label / legend) + 
    plot_layout(heights = c(1, 0.05, 0.1))
  
  return(final_plot)
}



## Functions for Supplemental figures
get_extended_palette_shades <- function(n) {
  base_pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7", "#999999")
  base_len <- length(base_pal)
  
  if (n <= base_len) {
    return(base_pal[1:n])
  } else {
    reps <- ceiling(n / base_len)
    pal <- character(0)
    
    for (i in 0:(reps-1)) {
      # Alternate darken/lighten per repetition
      if (i %% 2 == 0) {
        # Darken progressively for even repeats (first, third, etc.)
        shade_factor <- 1 - i * 0.25
        pal <- c(pal, darken(base_pal, amount = max(0, 1 - shade_factor)))
      } else {
        # Lighten progressively for odd repeats (second, fourth, etc.)
        shade_factor <- 0.15 + (i-1) * 0.25
        pal <- c(pal, lighten(base_pal, amount = min(1, shade_factor)))
      }
    }
    
    pal[1:n]  # truncate to exactly n colors
  }
}


# Comparing different transfromations
fig_compare_bias_trans = function(dt_list, dt_labels, r0, inst = TRUE, ylims = NA,
                                 types = "pq", title = "Bias", source_name = "source",
                                 ns = NULL, breaks_n = c(100, 5000, 10000, 20000), ...) {
  estimator_order <- c("p", types)  # this defines the order used in both plots
  if (is.null(names(dt_list))) names(dt_list) <- dt_labels
  if (length(dt_list) != length(dt_labels)) names(dt_list) <- dt_labels
  
  types <- as.character(types)
  est_cols <- paste0(types, ifelse(inst, "_inst", "_hold"))
  p_col <- ifelse(inst, "p_inst", "p_hold")
  
  plot_data = bind_rows(lapply(names(dt_list), function(name) {
    dt = dt_list[[name]]
    dt %>%
      group_by(n) %>%
      summarize(
        p = mean(.data[[p_col]]),
        across(all_of(est_cols), ~ mean(.x), .names = "{.col}"),
        .groups = "drop"
      ) %>%
      rename_with(.fn = function(x) {
        sapply(x, function(xx) {
          idx <- match(xx, est_cols)
          if (!is.na(idx)) return(types[idx]) else return(xx)
        }, USE.NAMES = FALSE)
      }, .cols = all_of(est_cols)) %>%
      mutate(source = name)
  }))
  
  if (!is.null(ns)) plot_data <- plot_data %>% filter(n %in% ns)
  
  long_cols <- c("p", types)
  plot_data_long = tidyr::pivot_longer(
    plot_data,
    cols = all_of(long_cols),
    names_to = "estimator",
    values_to = "value"
  )
  
  plot_data_long$estimator <- factor(plot_data_long$estimator, levels = estimator_order)
  
  
  sources <- unique(plot_data_long$source)
  shape_vals <- setNames( ((seq_along(sources) - 1) %% 6) + 1, sources)
  linetype_vals <- setNames(rep_len(1:3, length(sources)), sources)
  
  # Custom palette
  # Estimators (excluding 'truth')
  # Color mapping for bias plot
  pal_colors <- get_extended_palette_shades(length(estimator_order))
  estimator_colors <- setNames(pal_colors, estimator_order)
  estimator_colors <- c("truth" = "#000000", estimator_colors)  # add truth manually
  
  
  
  est_labels <- c("truth" = TeX("$\\rho_{target}$"))
  est_labels_rest <- setNames(
    c(TeX("$\\hat{\\rho}_{P}$"), sapply(types, function(t) TeX(sprintf("$\\hat{\\rho}_{%s}$", t)))),
    c("p", types)
  )
  est_labels <- c(est_labels, est_labels_rest[setdiff(names(est_labels_rest), names(est_labels))])
  
  # Define linetype mapping: make "p" and "pq" solid, others dashed
  linetype_vals <- c(
    "p" = "solid",
    "pq" = "solid"
  )
  # Any other estimators not listed default to "dashed"
  unique_ests <- unique(plot_data_long$estimator)
  missing_ests <- setdiff(unique_ests, names(linetype_vals))
  linetype_vals[missing_ests] <- "dashed"
  
  p <- ggplot(plot_data_long, aes(x = n, y = value, color = estimator, linetype = estimator)) +
    geom_hline(aes(yintercept = r0, color = "truth"), size = 1.5) +
    #geom_point(size = 3) +
    geom_line(size = 1.2, alpha = 0.75) +
    scale_x_continuous(name = "Model Sample Size", breaks = breaks_n) +
    scale_y_continuous(name = "Correlation") +
    scale_color_manual(name = NULL, values = estimator_colors, labels = est_labels, guide = guide_legend(ncol = 1)) +
    #scale_shape_manual(name = source_name, values = shape_vals) +
    scale_linetype_manual(values = linetype_vals, guide = "none") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.title = element_text(size = 22),
      legend.text = element_text(size = 22),
      legend.key.size = unit(40, "pt"),
      legend.position = "bottom"
    ) +
    ggtitle(title)
  
  if (!is.na(ylims[1])) p <- p + ylim(ylims)
  
  p
}



fig_compare_coverage_trans = function(dt_list, dt_labels, r0, cf = "inst", inst = TRUE,
                                     types = "pq", nsim = 500, title = "Coverage",
                                     source_name = "source", ns = NULL,
                                     breaks_n = c(100, 5000, 10000, 20000), ...) {
  #browser()
  estimator_order <- estimator_levels <- c("p", types)  # this defines the order used in both plots
  if (is.null(names(dt_list))) names(dt_list) <- dt_labels
  if (length(dt_list) != length(dt_labels)) names(dt_list) <- dt_labels
  
  types <- as.character(types)
  ci_type <- ifelse(inst, "inst", "hold")
  
  all_ns <- sort(unique(unlist(lapply(dt_list, function(dt) unique(dt$n)))))
  if (!is.null(ns)) all_ns <- intersect(all_ns, ns)
  
  ribbon_df <- data.frame(
    n = all_ns,
    ymin = qbinom(0.025, nsim, 0.95) / nsim,
    ymax = qbinom(0.975, nsim, 0.95) / nsim
  )
  
  coverage_data <- bind_rows(lapply(names(dt_list), function(name) {
    dt <- dt_list[[name]]
    
    cp_os_all <- get_cp(dt, types, r0) %>%
      filter(est %in% types, ci == cf) %>%
      transmute(n, coverage = coverage_probability, estimator = est, source = name)
    
    cp_p <- get_cp(dt, "p", r0) %>%
      filter(est == "p", ci == ci_type) %>%
      transmute(n, coverage = coverage_probability, estimator = "p", source = name)
    
    bind_rows(cp_p, cp_os_all)
  }))
  
  if (!is.null(ns)) coverage_data <- coverage_data %>% filter(n %in% ns)
  
  
  
  sources <- unique(coverage_data$source)
  shape_vals <- setNames(((seq_along(sources) - 1) %% 6) + 1, sources)
  linetype_vals <- setNames(rep_len(1:3, length(sources)), sources)
  
  coverage_data$estimator <- factor(coverage_data$estimator, levels = estimator_order)
  
  # Color mapping for coverage plot
  estimator_colors <- setNames(get_extended_palette_shades(length(estimator_order)), estimator_order)
  
  
  
  est_labels <- c("p" = TeX("$\\hat{\\rho}_{P}$"))
  for (t in types) est_labels[[t]] <- TeX(sprintf("$\\hat{\\rho}_{%s}$", t))
  est_labels <- est_labels[estimator_levels]
  
  linetype_vals <- c(
    "p" = "solid",
    "pq" = "solid"
  )
  # Any other estimators not listed default to "dashed"
  unique_ests <- unique(coverage_data$estimator)
  missing_ests <- setdiff(unique_ests, names(linetype_vals))
  linetype_vals[missing_ests] <- "dashed"
  
  ggplot(coverage_data, aes(x = n, y = coverage, color = estimator, linetype = estimator)) +
    geom_ribbon(data = ribbon_df, aes(x = n, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = "grey92") +
    geom_hline(yintercept = 0.95, color = "black", linetype = 3, size = 1.2) +
    #geom_point(size = 3) +
    geom_line(size = 1.2, alpha = 0.75) +
    scale_x_continuous(name = "Model Sample Size", breaks = breaks_n) +
    scale_y_continuous(name = "Coverage Probability", limits = c(0,1)) +
    scale_color_manual(name = NULL, values = estimator_colors, labels = est_labels) +
    #scale_shape_manual(name = source_name, values = shape_vals) +
    scale_linetype_manual(values = linetype_vals, guide = "none") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.text = element_text(size = 14),
      legend.key.size = unit(20, "pt"),
      legend.position = "bottom"
    ) +
    ggtitle(title)
}


# Panel for transformations (bias/coverage)
plot_panel_trans <- function(dt_suffixes, dt_labels, inst = TRUE, types = "pq", cf = "inst",
                                 nsim = 500, ylims = c(0.3, 1), ...) {
  # Create list of all data.frames
  dt_list = list(r10 = list(), r25 = list(), r50 = list(), r75 = list())
  i = 1
  for (s in dt_suffixes){
    dt_list$r10[[i]] = get(paste0("allResults_10_", s))
    dt_list$r25[[i]] = get(paste0("allResults_25_", s))
    dt_list$r50[[i]] = get(paste0("allResults_50_", s))
    dt_list$r75[[i]] = get(paste0("allResults_75_", s))
    i = i + 1
  }
  names(dt_list$r10) = names(dt_list$r25) = names(dt_list$r50) = names(dt_list$r75) = dt_labels
  
  # Titles for columns
  titles <- list(
    r10 = "MAPA = 0.1",
    r25 = "MAPA = 0.25",
    r50 = "MAPA = 0.5",
    r75 = "MAPA = 0.75"
  )
  
  # Bias plots (row 1)
  pb <- list(
    fig_compare_bias_trans(dt_list = dt_list$r10, dt_labels = dt_labels, r0 = 0.1,
                          inst = inst, types = types, title = titles$r10, ...),
    fig_compare_bias_trans(dt_list = dt_list$r25, dt_labels = dt_labels, r0 = 0.25,
                          inst = inst, types = types, title = titles$r25, ...),
    fig_compare_bias_trans(dt_list = dt_list$r50, dt_labels = dt_labels, r0 = 0.5,
                          inst = inst, types = types, title = titles$r50, ...),
    fig_compare_bias_trans(dt_list = dt_list$r75, dt_labels = dt_labels, r0 = 0.75,
                          inst = inst, types = types, title = titles$r75, ...)
  )
  
  # Coverage plots (row 2)
  pc <- list(
    fig_compare_coverage_trans(dt_list = dt_list$r10, dt_labels = dt_labels, r0 = 0.1,
                              inst = inst, types = types, title = NULL, nsim = nsim, cf = cf),
    fig_compare_coverage_trans(dt_list = dt_list$r25, dt_labels = dt_labels, r0 = 0.25,
                              inst = inst, types = types, title = NULL, nsim = nsim, cf = cf),
    fig_compare_coverage_trans(dt_list = dt_list$r50, dt_labels = dt_labels, r0 = 0.5,
                              inst = inst, types = types, title = NULL, nsim = nsim, cf = cf),
    fig_compare_coverage_trans(dt_list = dt_list$r75, dt_labels = dt_labels, r0 = 0.75,
                              inst = inst, types = types, title = NULL, nsim = nsim, cf = cf)
  )
  
  # Get legend from first bias plot
  legend <- get_legend(pb[[1]] + theme(legend.position = "right"))
  
  # Combine plots into a single list (Bias + Coverage)
  plots <- c(pb, pc)
  
  # Apply consistent theming
  plots_nolegend <- map2(plots, seq_along(plots), function(p, i) {
    row_idx <- ceiling(i / 4)
    col_idx <- ((i - 1) %% 4) + 1
    
    is_bottom_row <- row_idx == 2
    is_left_col <- col_idx == 1
    
    p + theme(
      legend.position = "none",
      axis.title.x = if (!is_bottom_row) element_blank() else element_blank(),
      axis.title.y = if (!is_left_col) element_blank() else element_text(size = 28),
      axis.text.x = element_text(size = 22),
      axis.text.y = element_text(size = 24),
      plot.title = element_text(size = 30),
      plot.margin = margin(5, 25, 5, 5)
    )
  })
  
  # Arrange in 2 rows × 4 columns
  plot_grid <- wrap_plots(plots_nolegend, nrow = 2, ncol = 4, byrow = TRUE)
  
  # X-axis label
  xaxis_label <- ggplot() +
    theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = "Model Sample Size", size = 8, hjust = 0.5)
  
  # Combine grid + legend + x-axis label
  final_plot <- (plot_grid / xaxis_label / legend) + 
    plot_layout(heights = c(1, 0.05, 0.1))
  #final_plot <- cowplot::plot_grid(plot_grid, legend, ncol = 2, rel_widths = c(4, 1))
  
  return(final_plot)
}

# Compare methods of aggregating estimates across folds
fig_compare_bias_inst_hold <- function(dt_list, dt_labels, r0, inst = TRUE, ylims = NA,
                                   types = "pq", title = "Bias", source_name = "source",
                                   ns = NULL, breaks_n = c(100, 5000, 10000, 20000), ...) {
  #browser()
  types <- as.character(types)
  
  # Create columns for both "_inst" and "_hold"
  est_cols <- c(paste0(types, "_inst"), paste0(types, "_hold"))
  p_cols <- c("p_inst", "p_hold")
  
  plot_data <- bind_rows(lapply(names(dt_list), function(name) {
    dt <- dt_list[[name]]
    dt %>%
      group_by(n) %>%
      summarize(
        across(all_of(c(p_cols, est_cols)), mean, .names = "{.col}"),
        .groups = "drop"
      ) %>%
      pivot_longer(
        cols = -n,
        names_to = "estimator",
        values_to = "value"
      ) %>%
      mutate(
        type = ifelse(grepl("_inst$", estimator), "inst", "hold"),
        estimator = gsub("_(inst|hold)$", "", estimator),
        source = name
      )
  }))
  
  if (!is.null(ns)) plot_data <- plot_data %>% filter(n %in% ns)
  
  estimator_order <- c("p", types)
  plot_data$estimator <- factor(plot_data$estimator, levels = estimator_order)
  
  # Colors
  pal_colors <- get_extended_palette_shades(length(estimator_order))
  estimator_colors <- setNames(pal_colors, estimator_order)
  estimator_colors <- c("truth" = "#000000", estimator_colors)
  
  # Labels
  est_labels <- c("truth" = TeX("$\\rho_{target}$"))
  est_labels_rest <- setNames(
    c(TeX("$\\hat{\\rho}_{P}$"), sapply(types, function(t) TeX(sprintf("$\\hat{\\rho}_{%s}$", t)))),
    c("p", types)
  )
  est_labels <- c(est_labels, est_labels_rest)
  
  # Linetypes: "_inst" = solid, "_hold" = dashed
  linetype_vals <- c("inst" = "solid", "hold" = "dashed")
  
  ggplot(plot_data, aes(x = n, y = value, color = estimator, linetype = type)) +
    geom_hline(aes(yintercept = r0, color = "truth"), size = 1.5) +
    geom_line(size = 1.2, alpha = 0.75) +
    scale_x_continuous(name = "Model Sample Size", breaks = breaks_n) +
    scale_y_continuous(name = "Correlation") +
    scale_color_manual(name = NULL, values = estimator_colors, labels = est_labels) +
    scale_linetype_manual(name = "Aggregation:", values = linetype_vals, labels = c("inst" = "Instant", "hold" = "Hold")) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.title = element_text(size = 22),
      legend.text = element_text(size = 22),
      legend.key.size = unit(40, "pt"),
      legend.position = "bottom"
    ) +
    ggtitle(title) +
    if(!is.na(ylims[1])) ylim(ylims) else NULL
}


plot_panel_inst_hold <- function(dt_suffixes, dt_labels, inst = TRUE, types = "pq", cf = "inst",
                             nsim = 500, ylims = c(0.3, 1), ...) {
  # Create list of all data.frames
  dt_list = list(r10 = list(), r25 = list(), r50 = list(), r75 = list())
  i = 1
  for (s in dt_suffixes){
    dt_list$r10[[i]] = get(paste0("allResults_10_", s))
    dt_list$r25[[i]] = get(paste0("allResults_25_", s))
    dt_list$r50[[i]] = get(paste0("allResults_50_", s))
    dt_list$r75[[i]] = get(paste0("allResults_75_", s))
    i = i + 1
  }
  names(dt_list$r10) = names(dt_list$r25) = names(dt_list$r50) = names(dt_list$r75) = dt_labels
  
  # Titles for columns
  titles <- list(
    r10 = "MAPA = 0.1",
    r25 = "MAPA = 0.25",
    r50 = "MAPA = 0.5",
    r75 = "MAPA = 0.75"
  )
  
  # Bias plots (row 1)
  pb <- list(
    fig_compare_bias_inst_hold(dt_list = dt_list$r10, dt_labels = dt_labels, r0 = 0.1,
                           inst = inst, types = types, title = titles$r10, ...),
    fig_compare_bias_inst_hold(dt_list = dt_list$r25, dt_labels = dt_labels, r0 = 0.25,
                           inst = inst, types = types, title = titles$r25, ...),
    fig_compare_bias_inst_hold(dt_list = dt_list$r50, dt_labels = dt_labels, r0 = 0.5,
                           inst = inst, types = types, title = titles$r50, ...),
    fig_compare_bias_inst_hold(dt_list = dt_list$r75, dt_labels = dt_labels, r0 = 0.75,
                           inst = inst, types = types, title = titles$r75, ...)
  )
  
  
  # Get legend from first bias plot
  legend <- get_legend(pb[[1]] + theme(legend.position = "right"))
  
  # Combine plots into a single list (Bias + Coverage)
  plots <- c(pb)
  
  # Apply consistent theming
  plots_nolegend <- map2(plots, seq_along(plots), function(p, i) {
    col_idx <- ((i - 1) %% 4) + 1
    
    is_left_col <- col_idx == 1
    
    p + theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = if (!is_left_col) element_blank() else element_text(size = 28),
      axis.text.x = element_text(size = 22),
      axis.text.y = element_text(size = 24),
      plot.title = element_text(size = 30),
      plot.margin = margin(5, 25, 5, 5)
    )
  })
  
  # Arrange in 2 rows × 4 columns
  plot_grid <- wrap_plots(plots_nolegend, nrow = 1, ncol = 4, byrow = TRUE)
  
  # X-axis label
  xaxis_label <- ggplot() +
    theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = "Model Sample Size", size = 8, hjust = 0.5)
  
  # Combine grid + legend + x-axis label
  final_plot <- (plot_grid / xaxis_label / legend) + 
    plot_layout(heights = c(1, 0.1, 0.25))
  #final_plot <- cowplot::plot_grid(plot_grid, legend, ncol = 2, rel_widths = c(4, 1))
  
  return(final_plot)
}


