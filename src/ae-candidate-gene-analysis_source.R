
# 0.0_candidate-gene-analysis_source.R
#
# This script provides the source code necessary to 
# check whether recorders are affecting egg laying 
# and make plots to this effect
# The functions below are very specific to this task and
# are very 'dirty'.
# 
# Copyright (c) Andrea Estandia, 2020, except where indicated
# Date Created: 2020-04-30


# --------------------------------------------------------------------------
# REQUIRES
# --------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(tidybayes)
  library(modelr)
  library(Hmisc)
  library(scales)
  library(devtools)
  library(jcolors)
  library(viridis)
  library(readxl)
  library(lubridate)
  library(janitor)
  library(lme4)
  library(ggrepel)
  library(RColorBrewer)
  library(patchwork)
  library(knitr)
  library(cowplot)
  library(genepop)
  library(sf)
  library(raster)
  library(dplyr)
  library(spData)
  library(spDataLarge)
  library(report)
  library(ggsci)
  library(patchwork)
  library(mcp)
  library(ggstance)
  library(clickR)
  library(gginnards)
  library(ggmap)
  library(maps)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(ggforce)
  library(ghibli)
  library(ggspatial)
  library(ape)
  library(ggpubr)
  library(diveRsity)

})

'%!in%' <- function(x,y)!('%in%'(x,y))

# --------------------------------------------------------------------------
# PATHS
# --------------------------------------------------------------------------

data_path <- file.path(getwd(), "data")
reports_path <- file.path(getwd(), "reports")
figures_path <- file.path(getwd(), "reports", "plots")

if (!dir.exists(data_path)) {
  dir.create(data_path, recursive = TRUE)
}

if (!dir.exists(reports_path)) {
  dir.create(reports_path, recursive = TRUE)
}

if (!dir.exists(figures_path)) {
  dir.create(figures_path, recursive = TRUE)
}


# --------------------------------------------------------------------------
# RESOURCES
# --------------------------------------------------------------------------

#' Data needed for map functions (map and geneflow_map)
island_info <- read_csv(file.path(data_path, "island_info.csv"))

sites_global <- island_info %>%
  filter(country %in% c("Australia", "New Zealand"))

sites_melanesia <- island_info %>%
  filter(country %in% c("Vanuatu", "New Caledonia"))

subset_australia <- subset(island_info,
                           country != "Vanuatu" &
                             country != "New Caledonia")

text_size=12

##COLOURS FOR MAP##
island_colour <- "#999999"
background_colour <- "#f0f0f0"
point_colour_global <- "#3d6b80"
point_colour_melanesia <- "#c99732"
myPalette <- colorRampPalette(rev(brewer.pal(9, "YlGnBu")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(1, 8))

# --------------------------------------------------------------------------
# FUNCTIONS
# --------------------------------------------------------------------------

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m - sd(x)
  ymax <- m + sd(x)
  return(c(y = m, ymin = ymin, ymax = ymax))
}

global_map <- function(data) {
  aus_map <- ggplot(data = data) +
  geom_sf(fill = island_colour, color = NA) +
  geom_point(
    data = sites_global,
    aes(x = lon, y = lat),
    size = 7,
    color = point_colour_global,
    alpha = 0.8
  ) +
  coord_sf(xlim = c(140, 185),
           ylim = c(-52,-12),
           expand = FALSE) +
  geom_text(data = island_info,
            aes(x = lon,
                y = lat,
                label = area),
            position = position_nudge(y = -1.5)) +
  theme(
    panel.background = element_rect(fill = background_colour),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")
  ) +
  annotate(
    "segment",
    x = 161.5,
    y = -12.5,
    xend = 161.5,
    yend = -24.6,
    alpha = .5,
    size = 1
  ) +
  annotate(
    "segment",
    y = -24.6,
    x = 172,
    xend = 172,
    yend = -12.5,
    alpha = .5,
    size = 1
  ) +
  annotate(
    "segment",
    y = -12.5,
    x = 172,
    xend = 161.5,
    yend = -12.5,
    alpha = .5,
    size = 1
  ) + annotate(
    "segment",
    y = -24.6,
    x = 172,
    xend = 161.5,
    yend = -24.6,
    alpha = .5,
    size = 1
  ) +
  annotation_scale(
    location = "bl", 
    width_hint = 0.12,
    style ="ticks"
  ) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    pad_x = unit(0.3, "in"),
    pad_y = unit(0.25, "in"),
    height = unit(1, "cm"),
    width = unit(1, "cm"),
    style = north_arrow_fancy_orienteering
  ) 

melanesia_map <- ggplot(data = data) +
  geom_sf(fill = island_colour, color = NA) +
  geom_point(
    data = sites_melanesia,
    aes(x = lon, y = lat),
    size = 7,
    color = point_colour_melanesia,
    alpha = 0.8
  ) +
  coord_sf(
    xlim = c(163, 171),
    ylim = c(-24.6, -12.5),
    expand = FALSE
  ) +
  geom_text_repel(
    data = subset(sites_melanesia,
                  island != "Espiritu Santo" &
                    island != "Malekula"),
    aes(x = lon,
        y = lat,
        label = island),
    segment.size = 0,
    nudge_x = 0.5,
    hjust = 0,
    vjust = 0,
  ) +
  geom_text(
    data = data.frame(
      x = 166.526517458422,
      y = -15.2067049176657,
      label = "Espiritu Santo"
    ),
    mapping = aes(x = x, y = y, label = label),
    size = 3.86605783866058,
    angle = 0L,
    lineheight = 1L,
    hjust = 1L,
    vjust = 0.5,
    colour = "black",
    family = "sans",
    fontface = "plain",
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_text(
    data = data.frame(
      x = 166.9,
      y = -16.1317272759486,
      label = "Malekula"
    ),
    mapping = aes(x = x, y = y, label = label),
    size = 3.86605783866058,
    angle = 0L,
    lineheight = 1L,
    hjust = 1L,
    vjust = 0.5,
    colour = "black",
    family = "sans",
    fontface = "plain",
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  theme(
    panel.background = element_rect(fill = background_colour),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")
  ) + annotation_scale(
    location = "bl", 
    width_hint = 0.12,
    style = "ticks"
  )

both_maps <- aus_map + melanesia_map

both_maps <- both_maps + annotate(
  "segment",
  x = 0.116302215576172,
  y = 0.641714284624372,
  xend = 0.262419862634995,
  yend = 0.384571427481515,
  colour = "black"
) 
return(both_maps)
}

#############################################

geneflow_map <- function(data) {
  map <- ggplot(data = data) +
  geom_sf(fill = island_colour, color = NA) +
  geom_point(
    data = subset_australia,
    aes(x = lon, y = lat),
    size = 7,
    color = point_colour_global,
    alpha = 0.6
  ) +
  coord_sf(xlim = c(140, 184),
           ylim = c(-48,-22),
           expand = FALSE) +
  geom_text(data = subset_australia,
            aes(x = lon,
                y = lat,
                label = area),
            position = position_nudge(y = -1.5)) +
  theme(
    panel.background = element_rect(fill = "#e3e3e3"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")
  ) +
  annotate(
    "curve",
    yend = -27.5,
    xend = 152,
    y = -42.6,
    x = 147,
    size = 0.5,
    curvature = -.5,
    alpha = 0.7,
    colour = "#404040",
    arrow = arrow(length = unit(2, "mm"))
  ) +
  annotate(
    "curve",
    yend = -45.9,
    xend = 169,
    y = -42.8,
    x = 147,
    size = 0.5,
    curvature = 0,
    alpha = 0.7,
    colour = "#404040",
    arrow = arrow(length = unit(2, "mm"))
  ) +
  annotate(
    "curve",
    yend = -44.4,
    xend = 178,
    y = -42.8,
    x = 147,
    size = 0.5,
    curvature = 0,
    alpha = 0.7,
    colour = "#404040",
    arrow = arrow(length = unit(2, "mm"))
  ) +
  geom_text(
    data = data.frame(x = 144, y = -32,
                      label = "0.22"),
    mapping = aes(x = x, y = y, label = label),
    size = 3.6,
    angle = 0L,
    lineheight = 1L,
    hjust = 0.5,
    vjust = 0.5,
    colour = "black",
    family = "sans",
    fontface = "plain",
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_text(
    data = data.frame(
      x = 158.742767273088,
      y = -42.3619526042132,
      label = "0.14"
    ),
    mapping = aes(x = x, y = y, label = label),
    size = 3.6,
    angle = 0L,
    lineheight = 1L,
    hjust = 0.5,
    vjust = 0.5,
    colour = "black",
    family = "sans",
    fontface = "plain",
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_text(
    data = data.frame(x = 158, y = -45.5,
                      label = "0.14"),
    mapping = aes(x = x, y = y, label = label),
    size = 3.6,
    angle = 0L,
    lineheight = 1L,
    hjust = 0.5,
    vjust = 0.5,
    colour = "black",
    family = "sans",
    fontface = "plain",
    inherit.aes = FALSE,
    show.legend = FALSE
  )
  
  map2 <- ggplot(data = data) +
    geom_sf(fill = island_colour, color = NA) +
    geom_point(
      data = sites_melanesia,
      aes(x = lon, y = lat),
      size = 7,
      color = point_colour_melanesia,
      alpha = 0.8
    ) +
    coord_sf(
      xlim = c(163, 171),
      ylim = c(-23, -12.5),
      expand = FALSE
    ) +
    geom_text_repel(
      data = subset(sites_melanesia,
                    island != "Espiritu Santo" &
                      island != "Malekula"),
      aes(x = lon,
          y = lat,
          label = island),
      segment.size = 0,
      nudge_x = 0.5,
      hjust = 0,
      vjust = 0,
    ) +
    geom_text(
      data = data.frame(
        x = 166.526517458422,
        y = -15.2067049176657,
        label = "Espiritu Santo"
      ),
      mapping = aes(x = x, y = y, label = label),
      size = 3.86605783866058,
      angle = 0L,
      lineheight = 1L,
      hjust = 1L,
      vjust = 0.5,
      colour = "black",
      family = "sans",
      fontface = "plain",
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    geom_text(
      data = data.frame(
        x = 166.9,
        y = -16.1317272759486,
        label = "Malekula"
      ),
      mapping = aes(x = x, y = y, label = label),
      size = 3.86605783866058,
      angle = 0L,
      lineheight = 1L,
      hjust = 1L,
      vjust = 0.5,
      colour = "black",
      family = "sans",
      fontface = "plain",
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    theme(
      panel.background = element_rect(fill = "#e3e3e3"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")
    ) +
    annotate(
      "segment",
      yend = -18.8,
      xend = 169,
      y = -17.8,
      x = 168.3,
      size = 0.5,
      alpha = 0.7,
      colour = "#404040",
      arrow = arrow(length = unit(2, "mm"))
    ) +
    annotate(
      "segment",
      yend = -20.9,
      xend = 167.2,
      y = -20.4,
      x = 166.7,
      size = 0.5,
      alpha = 0.7,
      colour = "#404040",
      arrow = arrow(length = unit(2, "mm"))
    ) +
    #Ambrym
    annotate(
      "segment",
      yend = -16.1,
      xend = 167.85,
      y = -15.3,
      x = 167.2,
      size = 0.5,
      alpha = 0.7,
      colour = "#404040",
      arrow = arrow(length = unit(2, "mm"))
    ) +
    #Malekula
    annotate(
      "segment",
      yend = -16,
      xend = 167.25,
      y = -15.3,
      x = 167.2,
      size = 0.5,
      alpha = 0.7,
      colour = "#404040",
      arrow = arrow(length = unit(2, "mm"))
    ) +
    #Gaua
    annotate(
      "segment",
      yend = -14.3,
      xend = 167.55,
      y = -15.3,
      x = 167.2,
      size = 0.5,
      alpha = 0.7,
      colour = "#404040",
      arrow = arrow(length = unit(2, "mm"))
    ) +
    #Pentecost
    annotate(
      "segment",
      yend = -15.6,
      xend = 168,
      y = -15.3,
      x = 167.2,
      size = 0.5,
      alpha = 0.7,
      colour = "#404040",
      arrow = arrow(length = unit(2, "mm"))
    ) +
    geom_text(
      data = data.frame(
        x = 166.7,
        y = -20.8856387456712,
        label = "0.14"
      ),
      mapping = aes(x = x, y = y, label = label),
      size = 3.6,
      angle = 0L,
      lineheight = 1L,
      hjust = 0.5,
      vjust = 0.5,
      colour = "black",
      family = "sans",
      fontface = "plain",
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    geom_text(
      data = data.frame(
        x = 169,
        y = -18.2227990219934,
        label = "0.16"
      ),
      mapping = aes(x = x, y = y, label = label),
      size = 3.6,
      angle = 0L,
      lineheight = 1L,
      hjust = 0.5,
      vjust = 0.5,
      colour = "black",
      family = "sans",
      fontface = "plain",
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    geom_text(
      data = data.frame(
        x = 166.6,
        y = -15.8,
        label = "0.17"
      ),
      mapping = aes(x = x, y = y, label = label),
      size = 3.6,
      angle = 0L,
      lineheight = 1L,
      hjust = 0.5,
      vjust = 0.5,
      colour = "black",
      family = "sans",
      fontface = "plain",
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    geom_text(
      data = data.frame(
        x = 167.1,
        y = -14.5,
        label = "0.07"
      ),
      mapping = aes(x = x, y = y, label = label),
      size = 3.6,
      angle = 0L,
      lineheight = 1L,
      hjust = 0.5,
      vjust = 0.5,
      colour = "black",
      family = "sans",
      fontface = "plain",
      inherit.aes = FALSE,
      show.legend = FALSE
    ) + 
    geom_text(
      data = data.frame(
        x = 167.8,
        y = -15.7929577741375,
        label = "**"
      ),
      mapping = aes(x = x, y = y, label = label),
      size = 3.86605783866058,
      angle = 0L,
      lineheight = 1L,
      hjust = 0.5,
      vjust = 0.5,
      colour = "black",
      family = "sans",
      fontface = "plain",
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    geom_text(
      data = data.frame(x = 167.6, y = -15.3,
                        label = "*"),
      mapping = aes(x = x, y = y, label = label),
      size = 3.6,
      angle = 0L,
      lineheight = 1L,
      hjust = 0.5,
      vjust = 0.5,
      colour = "black",
      family = "sans",
      fontface = "plain",
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    geom_text(
      data = data.frame(x = 170, y = -13,
                        label = "* = 0.14\n** = 0.13"),
      mapping = aes(x = x, y = y, label = label),
      size = 3.6,
      angle = 0L,
      lineheight = 1L,
      hjust = 0.5,
      vjust = 0.5,
      colour = "black",
      family = "sans",
      fontface = "plain",
      inherit.aes = FALSE,
      show.legend = FALSE
    ) 
  geneflow_map <- map + map2
  return(geneflow_map)
}

# ----------------------------------------------------------------------------------- #
# Function adapted from the mcp source code (Liedvogel 2020) by Andrea Estandia 2020
# ----------------------------------------------------------------------------------- #


#' Plot full fits
#'
#' Plot prior or posterior model draws on top of data. Use `plot_pars` to
#' plot individual parameter estimates.
#'
#' @aliases plot plot.mcpfit
#' @param x An \code{\link{mcpfit}} object
#' @param facet_by String. Name of a varying group.
#' @param lines Positive integer or `FALSE`. Number of lines (posterior
#'   draws). FALSE or `lines = 0` plots no lines. Note that lines always plot
#'   fitted values - not predicted. For prediction intervals, see the
#'   `q_predict` argument.
#' @param geom_data String. One of "point" (default), "line" (good for time-series),
#'   or FALSE (don not plot).
#' @param cp_dens TRUE/FALSE. Plot posterior densities of the change point(s)?
#'   Currently does not respect `facet_by`. This will be added in the future.
#' @param q_fit Whether to plot quantiles of the posterior (fitted value).
#'   * \strong{TRUE:} Add 2.5% and 97.5% quantiles. Corresponds to
#'       `q_fit = c(0.025, 0.975)`.
#'   * \strong{FALSE (default):} No quantiles
#'   * A vector of quantiles. For example, `quantiles = 0.5`
#'       plots the median and `quantiles = c(0.2, 0.8)` plots the 20% and 80%
#'       quantiles.
#' @param q_predict Same as `q_fit`, but for the prediction interval.
#' @param rate Boolean. For binomial models, plot on raw data (`rate = FALSE`) or
#'   response divided by number of trials (`rate = TRUE`). If FALSE, linear
#'   interpolation on trial number is used to infer trials at a particular x.
#' @param prior TRUE/FALSE. Plot using prior samples? Useful for `mcp(..., sample = "both")`
#' @param which_y What to plot on the y-axis. One of
#'
#'   * `"ct"`: The central tendency which is often the mean after applying the
#'     link function (default).
#'   * `"sigma"`: The variance
#'   * `"ar1"`, `"ar2"`, etc. depending on which order of the autoregressive
#'     effects you want to plot.
#'
#' @param ... Currently ignored.
#' @details
#'   `plot()` uses `fit$simulate()` on posterior samples. These represent the
#'   (joint) posterior distribution.
#' @return A \pkg{ggplot2} object.
#' @encoding UTF-8
#' @author Jonas Kristoffer Lindeløv \email{jonas@@lindeloev.dk}
#' @importFrom ggplot2 ggplot aes aes_string geom_line geom_point facet_wrap
#' @importFrom magrittr %>%
#' @importFrom rlang !! :=
#' @importFrom dplyr .data
#' @export
#' @examples
#' # Typical usage. ex_fit is an mcpfit object.
#' plot(ex_fit)
#' plot(ex_fit, prior = TRUE)  # The prior
#'
#' \donttest{
#' plot(ex_fit, lines = 0, q_fit = TRUE)  # 95% HDI without lines
#' plot(ex_fit, q_predict = c(0.1, 0.9))  # 80% prediction interval
#' plot(ex_fit, which_y = "sigma", lines = 100)  # The variance parameter on y
#' }
#'
#' # Show a panel for each varying effect
#' # plot(fit, facet_by = "my_column")
#'
#' # Customize plots using regular ggplot2
#' library(ggplot2)
#' plot(ex_fit) + theme_bw(15) + ggtitle("Great plot!")
#'
plot.mcpfit.andrea = function(x,
                       facet_by = NULL,
                       lines = 25,
                       geom_data = "point",
                       cp_dens = TRUE,
                       q_fit = FALSE,
                       q_predict = FALSE,
                       rate = TRUE,
                       prior = FALSE,
                       which_y = "ct",
                       reverse_xaxis = FALSE,
                       ...) {
  
  # Just for consistent naming in mcp
  fit = x
  
  # Check arguments
  if (class(fit) != "mcpfit")
    stop("Can only plot mcpfit objects. x was class: ", class(fit))
  
  if (!coda::is.mcmc.list(fit$mcmc_post) & !coda::is.mcmc.list(fit$mcmc_prior))
    stop("Cannot plot an mcpfit without prior or posterior samples.")
  
  if (lines != FALSE) {
    check_integer(lines, "lines", lower = 1)
  } else {
    lines = 0
  }
  
  if (lines > 1000) {
    lines = 1000
    warning("Setting `lines` to 1000 (maximum).")
  }
  
  if (!geom_data %in% c("point", "line", FALSE))
    stop("`geom_data` has to be one of 'point', 'line', or FALSE.")
  
  if (!is.logical(cp_dens))
    stop("`cp_dens` must be TRUE or FALSE.")
  
  if (all(q_fit == TRUE))
    q_fit = c(0.025, 0.975)
  
  if (all(q_predict == TRUE))
    q_predict = c(0.025, 0.975)
  
  if (!is.logical(q_fit) & !is.numeric(q_fit))
    stop("`q_fit` has to be TRUE, FALSE, or a vector of numbers.")
  
  if (is.numeric(q_fit) & (any(q_fit > 1) | any(q_fit < 0)))
    stop ("All `q_fit` have to be between 0 (0%) and 1 (100%).")
  
  if (!is.logical(q_predict) & !is.numeric(q_predict))
    stop("`q_predict` has to be TRUE, FALSE, or a vector of numbers.")
  
  if (is.numeric(q_predict) & (any(q_predict > 1) | any(q_predict < 0)))
    stop ("All `q_predict` have to be between 0 (0%) and 1 (100%).")
  
  if (!is.logical(rate))
    stop("`rate` has to be TRUE or FALSE.")
  
  if (!is.logical(prior))
    stop("`prior` must be either TRUE or FALSE.")
  
  # Is facet_by a random/nested effect?
  varying_groups = logical0_to_null(unique(stats::na.omit(fit$.other$ST$cp_group_col)))
  if (!is.null(facet_by)) {
    if (!facet_by %in% varying_groups)
      stop("`facet_by` is not a data column used as varying grouping.")
  }
  
  
  # R CMD Check wants a global definition of ".". The formal way of doing it is
  # if(getRversion() >= "2.15.1") utils::globalVariables(".")
  # but that makes the tests fail.
  . = "ugly fix to please R CMD check"
  
  # Select posterior/prior samples
  samples = get_samples(fit, prior = prior)
  
  # General settings
  xvar = rlang::sym(fit$pars$x)
  yvar = rlang::sym(fit$pars$y)
  simulate = fit$simulate  # To make a function call work later
  dens_threshold = 0.001  # Do not display change point densities below this threshold.
  dens_height = 0.2  # proportion of plot y-axis that the change point density makes up
  if (all(q_fit == FALSE) & all(q_predict == FALSE)) {
    HDI_SAMPLES = lines
  } else {
    HDI_SAMPLES = 1000 # Maximum number of draws to use for computing HDI
    HDI_SAMPLES = min(HDI_SAMPLES, length(samples) * nrow(samples[[1]]))
  }
  is_arma = length(fit$pars$arma) > 0
  if (is_arma & (all(q_fit != FALSE) | all(q_predict != FALSE)))
    message("Plotting ar() with quantiles can be slow. Raise an issue at GitHub (or thumb-up existing ones) if you need this.")
  if (!is.null(facet_by)) {
    n_facet_levels = length(unique(fit$data[, facet_by]))
  } else {
    n_facet_levels = 1
  }
  
  #################
  # GET PLOT DATA #
  #################
  regex_pars_pop = paste0(fit$pars$population, collapse="|")
  
  # No faceting
  if (is.null(facet_by)) {
    samples = samples %>%
      tidybayes::spread_draws(!!rlang::sym(regex_pars_pop), regex = TRUE)
    
  } else {
    # Prepare for faceting
    # Read more about this weird syntax at https://github.com/mjskay/tidybayes/issues/38
    varying_by_facet = stats::na.omit(fit$.other$ST$cp_group[stringr::str_detect(fit$.other$ST$cp_group, paste0("_", facet_by, "$"))])
    varying_by_facet = paste0(varying_by_facet, collapse="|")
    
    samples = samples %>%
      tidybayes::spread_draws(!!rlang::sym(regex_pars_pop),
                              (!!rlang::sym(varying_by_facet))[!!rlang::sym(facet_by)],
                              regex = TRUE)
  }
  
  # Remove some samples
  samples = tidybayes::sample_draws(samples, n = HDI_SAMPLES)  # TO DO: use spread_draws(n = draws) when tidybayes 1.2 is out
  
  # Get x-coordinates to evaluate simulate (etc.) at
  eval_at = get_eval_at(fit, facet_by)
  
  # First, let's get all the predictors in shape for simulate
  if (fit$family$family != "binomial") {
    samples = samples %>%
      tidyr::expand_grid(!!xvar := eval_at)  # correct name of x-var
  } else if (fit$family$family == "binomial") {
    if (!is.null(facet_by) & rate == FALSE)
      stop("Plot with rate = FALSE not implemented for varying effects (yet).")
    
    # Interpolate trials for binomial at the values in "eval_at"
    # to make sure that there are actually values to plot
    #interpolated_trials = suppressWarnings(stats::approx(x = fit$data[, fit$pars$x], y = fit$data[, fit$pars$trials], xout = eval_at)$y)
    interpolated_trials = stats::approx(x = fit$data[, fit$pars$x], y = fit$data[, fit$pars$trials], xout = eval_at)$y %>%
      suppressWarnings() %>%
      round()  # Only integers
    
    samples = samples %>%
      tidyr::expand_grid(!!xvar := eval_at) %>%  # correct name of x-var
      dplyr::mutate(!!fit$pars$trials := rep(interpolated_trials, nrow(samples)))
  }
  
  # For ARMA prediction, we need the raw data
  # We know that eval_at is the same length as nrow(data), so we can simply add corresponding data$y for each draw
  if (is_arma) {
    samples = dplyr::mutate(samples, ydata = rep(fit$data[, fit$pars$y], HDI_SAMPLES * n_facet_levels))
  }
  
  # Predict y from model by adding fitted/predicted draws (vectorized)
  if (lines > 0 | (any(q_fit != FALSE))) {
    samples = samples %>%
      dplyr::mutate(!!yvar := rlang::exec(simulate, !!!., type = "fitted", rate = rate, which_y = which_y, add_attr = FALSE))
  }
  if (any(q_predict != FALSE)) {
    samples = samples %>%
      dplyr::mutate(predicted_ = rlang::exec(simulate, !!!., type = "predict", rate = rate, add_attr = FALSE))
  }
  
  
  
  ###########
  # PLOT IT #
  ###########
  # If this is a binomial rate, divide by the number of trials
  if (fit$family$family == "binomial" & rate == TRUE) {
    fit$data[, fit$pars$y] = fit$data[, fit$pars$y] / fit$data[, fit$pars$trials]
  }
  
  # Initiate plot.
  gg = ggplot(fit$data,
              aes_string(x = xvar,
                         y = yvar)) +
    geom_jitter(
      aes(group = cluster,
        shape =cluster),
      size=2,
      width = 0.002,
      height = 0,
      alpha = 0.6,
      colour= point_colour_melanesia)

  
  if (which_y == "ct") {
    if (geom_data == "line")
      gg = gg + geom_line()
  }

  
  # Add lines?
  if (lines > 0) {
    # Sample the desired number of lines
    data_lines = samples %>%
      tidybayes::sample_draws(250) %>%
      dplyr::mutate(
        # Add line ID to separate lines in the plot.
        line = !!xvar == min(!!xvar),
        line = cumsum(.data$line)
      )
    
    # Plot it
    gg = gg + geom_line(aes(group = .data$line), data = data_lines, alpha = 0.1, color = grDevices::rgb(0.5, 0.5, 0.5, 0.4))
  }
  
  # Add quantiles?
  if ((any(q_fit != FALSE))) {
    samples_fit = dplyr::mutate(samples, y_quant = !!yvar)
    gg = gg + geom_quantiles(samples_fit, q_fit, xvar, facet_by, color = "red")
  }
  if (any(q_predict != FALSE)) {
    samples_predict = dplyr::mutate(samples, y_quant = .data$predicted_)
    gg = gg + geom_quantiles_andrea(samples_predict, q_predict, xvar, facet_by, alpha=0.4,
                                    size=0.8,
                                    color = "#4c6173")
  }
  
  # Add change point densities?
  if (cp_dens == TRUE & length(fit$model) > 1) {
    # The scale of the actual plot (or something close enough)
    if (which_y == "ct" & geom_data != FALSE) {
      y_data_max = max(fit$data[, fit$pars$y])
      y_data_min = min(fit$data[, fit$pars$y])
    } else if (any(q_predict != FALSE)) {
      y_data_max = max(samples$predicted_)
      y_data_min = min(samples$predicted_)
    } else if (as.character(yvar) %in% names(samples)) {
      y_data_max = max(dplyr::pull(samples, as.character(yvar)))
      y_data_min = min(dplyr::pull(samples, as.character(yvar)))
    } else {
      stop("Failed to draw change point density for this plot. Please raise an error on GitHub.")
    }
    
    # Function to get density for each grouping in the dplyr pipes below.
    density_xy = function(x) {
      tmp = stats::density(x)
      df = data.frame(x_dens = tmp$x, y_dens = tmp$y) %>%
        dplyr::filter(.data$y_dens > dens_threshold) %>%
        return()
    }
    
    # Get and group samples to be used for computing density
    cp_regex = "^cp_[0-9]+$"
    cp_dens_xy = get_samples(fit, prior = prior) %>%  # Use all samples for this
      tidybayes::gather_draws(!!rlang::sym(cp_regex), regex = TRUE) %>%
      dplyr::group_by(.data$.chain, add = TRUE) %>%
      
      # Get density as x-y values by chain and change point number.
      dplyr::summarise(dens = list(density_xy(.data$.value))) %>%
      tidyr::unnest(cols = c(.data$dens)) %>%
      
      # Then scale to plot.
      dplyr::mutate(y_dens = y_data_min + .data$y_dens * (y_data_max - y_data_min) * dens_height / max(.data$y_dens))
    
    # Add cp density to plot
    gg = gg + ggplot2::geom_line(
      data = cp_dens_xy,
      mapping = aes(
        x = .data$x_dens,
        y = .data$y_dens,
        color = .data$.chain,
        group = interaction(.data$.variable, .data$.chain))) +
      ggplot2::theme(legend.position = "none") +
      sc 
    
  }
  
  #Reverse x axis
  if (reverse_xaxis == TRUE) {
    gg = gg + scale_x_reverse()
  return(gg)
  }
  
  # Add faceting?
  if (!is.null(facet_by)) {
    gg = gg + facet_wrap(paste0("~", facet_by))
  }
  
  # Add better y-labels
  if (fit$family$family == "bernoulli" | (fit$family$family == "binomial" & rate == TRUE))
    gg = gg + ggplot2::labs(y = paste0("Probability of success for ", fit$pars$y))
  if (which_y != "ct")
    gg = gg + ggplot2::labs(y = which_y)
  
  return(gg)
}

### A custom function adapted for re-designing the output plot in mcp ###

geom_quantiles_andrea = function(samples, q, xvar, facet_by, ...) {
  # Trick to declare no facet = common group for all
  if (length(facet_by) == 0)
    facet_by = xvar
  
  # First: add quantiles column
  data_quantiles = samples %>%
    tidyr::expand_grid(quant = q) %>%
    
    # Now compute the quantile for each parameter, quantile, and (optionally) facet:
    dplyr::group_by(!!xvar, .data$quant) %>%
    dplyr::group_by(!!rlang::sym(facet_by), add = TRUE) %>%
    dplyr::summarise(
      y = stats::quantile(.data$y_quant, probs = .data$quant[1])
    )
  
  # Return geom
  geom = ggplot2::geom_line(
    mapping = aes(
      y = .data$y,
      group = .data$quant
    ),
    data = data_quantiles,
    linetype = "dotted",
    ...)
  return(geom)
}
environment(plot.mcpfit.andrea) <- asNamespace('mcp')
environment(geom_quantiles_andrea) <- asNamespace('mcp')
