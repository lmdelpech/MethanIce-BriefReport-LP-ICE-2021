# ggPCA() function
#
# A function to draw a PCA biplot (scaling 1 or scaling 2) from an object
# of class "rda" containing PCA results, computed with vegan's rda() function.
# A circle of equilibrium contribution is drawn on scaling type 1 biplots.
#
# This function is adapted from the cleanplot.pca() function from Legendre & Legendre, 2012
# to use {ggplot}
# It uses the same scalings and streching of species scores
#
# ARGUMENTS
#
# ##### General parameters
# pca.out          An rda{vegan} object. This object contains a PCA result if a
#                  single matrix was used in the rda(). The function will still
#                  operate with aRDA result but a Warning will be printed.
# ax1, ax2         Canonical axes to be drawn as abscissa and ordinate.
#                  Defaults: 1 and 2.
# scaling          Scaling type: only 1 or 2 are supported. Default: 1.
#
# ##### Items to be plotted
# plot.sites       If TRUE, the sites will be plotted as small circles.
# plot.spe         If TRUE, the species (or other response variables) will be
#                  plotted.
# label.sites      If TRUE, labels are added to the site symbols.
# label.spe        If TRUE, labels are added to the species arrows.
#
#
# ##### Multipliers, selection of species to be plotted
# mult.spe         Multiplier for length of the species arrows. Default: 1.
# select.spe       Vector containing a selection of the species numbers to be
#                  drawn in the biplot, e.g. c(1, 2, 5, 8, 12). Draw all species
#                  if select.spe = NULL (default value). The species that are
#                  well represented in the RDA plot can be identified using
#                  goodness(pca.output.object, display = "species").
#
# ##### Position of the plot in frame, margins
# mar.percent      Factor to expand plot size to accomodate all objects and
#                  labels. Positive values increase the margins around the plot,
#                  negative values reduce them.
# optimum          If TRUE, the longest species arrow is stretched to
#                  a length equal to the distance to the origin of the site
#                  farthest from the origin of the plot of (ax1, ax2). This is
#                  an optimal combined representation of the sites and species.
#                  The lengths of the species arrows can be further modified
#                  using the arguments mult.spe.
# move.origin      Move plot origin right-left and up-down.
#                  Default: move.origin = c(0,0).
#                  Ex. move.origin = c(-1, 0.5) moves origin by 1 unit left and
#                  0.5 unit up. # Not implemented in this function
#
# ##### Colors and shape of items
# color.vec        Is the vector to give matching the number of sites, with
#                  the variable used to color the sites (argument in aes(col)
#                  in ggplot)
# shape.vec        Is the vector to give matching the number of sites, with
#                  the variable used to give a shape to the sites (argument in
#                  aes(shape) in ggplot)
# ...              Additional arguments passed on to geom_point() for plotting sites
#
# ##### Varia
# silent           If FALSE, intermediate computation steps are printed.
#                  Default: TRUE. Interpretation: examine the code lines
#                  controlled by argument "silent".
#
# Reference
# Legendre, P. & L. Legendre. 2012. Numerical ecology, 3rd English edition.
#    Elsevier Science BV, Amsterdam.
# Adapted by Lisa-Marie Delpech



require(vegan)
require(ggplot2)
require(ggrepel)
# require(ggforce)


'ggPCA' <-
  function(pca.out,
           ax1 = 1,
           ax2 = 2,
           scaling = 1,
           plot.sites = TRUE,
           plot.spe = TRUE,
           label.sites = TRUE,
           label.spe = TRUE,
           draw.circle = TRUE,
           mult.spe = 1,
           select.spe = NULL,
           # mar.percent = 0.1,
           color.vec = NULL,
           shape.vec = NULL,
           fill.vec = NULL, # If needed when using filled symbols
           optimum = TRUE,
           spe.names = function(x){x}, # A function to rename variables
           # move.origin = c(0, 0),
           silent = TRUE,
           plot=TRUE,
           ...) {


####### Internal function start #######

    'stretch' <-
      function(sites, mat, ax1, ax2, n, silent = silent) {
        # Compute stretching factor for the species arrows
        # First, compute the longest distance to centroid for the sites
        tmp1 <- rbind(c(0, 0), sites[, c(ax1, ax2)])
        D <- dist(tmp1)
        target <- max(D[1:n])
        # Then, compute the longest distance to centroid for the species arrows
        if (is.matrix(mat)) {
          p <- nrow(mat)   # Number of species to be drawn
          tmp2 <- rbind(c(0, 0), mat[, c(ax1, ax2)])
          D <- dist(tmp2)
          longest <- max(D[1:p])
        } else {
          tmp2 <- rbind(c(0, 0), mat[c(ax1, ax2)])
          longest <- dist(tmp2)
          # print(tmp2)
        }  # If a single row left in 'mat'
        #
        if (!silent)
          cat("target =",
              target,
              " longest =",
              longest,
              " fact =",
              target / longest,
              "\n")
        fact <- target / longest
      }

    'pca.circle' <- function(pca,
                             mult.spe = 1,
                             fact.spe = 1,
                             center = c(0,0), # if we move the origin, we can modify this too
                             npoints = 100){
      eigenv <- pca.out$CA$eig
      p <- length(eigenv)
      radius <- (2 / p) ^ 0.5 * mult.spe * fact.spe
      tt <- seq(0,2*pi,length.out = npoints)
      xx <- center[1] + radius * cos(tt)
      yy <- center[2] + radius * sin(tt)
      return(data.frame(x = xx, y = yy))
    }

####### Internal function end ########


    if (!class(pca.out)[1] == "rda")
      stop("The input file is not a vegan output object of class 'rda'",
           call. = FALSE)
    if (!(is.null(pca.out$CCA)))
      stop(
        "The input file contains an RDA, not a PCA result. ",
        "Use function triplot.rda from the NEwR (2018) book to produce an RDA triplot."
      )
    if (scaling != 1 &
        scaling != 2)
      stop("Function only available for scaling 1 or 2", call. = FALSE)


    k <- length(pca.out$CA$eig)         # number of PCA eigenvalues (or axes)
    n.sp <- length(pca.out$colsum)      # number of species

    # 'vec' will contain the selection of species to be drawn
    if (is.null(select.spe)) {
      vec <- 1:n.sp
    } else {
      vec <- select.spe
    }


# df.pca <- data.frame(scores(pca.out, choices = c(1:3), display = c("sites"), scaling = 1)) # Get the PC for each site as a data frame
# df.pca.species <- data.frame(scores(pca.out, choices = c(1:3), display = c("species"), scaling = 1)) # Get the PC for each species as a data frame


    # Scaling 1: the species scores have norms of 1
    # Scaling 1: the site scores are scaled to variances = can.eigenvalues
    # Scaling 2: the species scores have norms of sqrt(can.eigenvalues)
    # Scaling 2: the site scores are scaled to variances of 1

    # This version reconstructs and uses the original RDA output of L&L 2012,
    # Section 11.1.3

    Tot.var = pca.out$tot.chi         # Total variance in response data Y
    eig.val = pca.out$CA$eig          # Eigenvalues of Y-hat
    Lambda = diag(eig.val)            # Diagonal matrix of eigenvalues
    eig.val.rel = eig.val / Tot.var   # Relative eigenvalues of Y-hat
    Diag = diag(sqrt(eig.val.rel))    # Diagonal matrix of sqrt(relative eigenvalues)
    #const = ((n-1)*Tot.var)^(1/4)     # The constant used in the scaling from the vegan package
    U.sc1 = pca.out$CA$v              # Species scores, scaling=1 # Why no multiplying by the constant sqrt(n-1)?
    U.sc2 = U.sc1 %*% Lambda ^ (0.5)  # Species scores, scaling=2
    n = nrow(pca.out$CA$u)            # Number of observations
    Z.sc2 = pca.out$CA$u * sqrt(n - 1)# Site scores, scaling=2 NB. scaling in vegan is different
    Z.sc1 = Z.sc2 %*% Lambda ^ (0.5)  # Site scores, scaling=1

    if (scaling == 1) {
      sit.sc <- Z.sc1
      spe.sc <- U.sc1[vec,]
    } else {
      # For scaling=2
      sit.sc <- Z.sc2
      spe.sc <- U.sc2[vec,]
    }

    fact.spe <- 1
    if (optimum & (plot.spe | label.spe)) { # Compute stretching factor for the species arrows (only if plotting species)
      fact.spe <-
        stretch(sit.sc[, 1:k], spe.sc[, 1:k], ax1, ax2, n, silent = silent) # k = number of PCA eigenvalues (or axes)
    }
    spe.sc <- spe.sc * fact.spe * mult.spe

    # Draw the main plot

    df.pca <- data.frame(sit.sc[, 1:k])
    df.pca.sp <- data.frame(spe.sc[, 1:k])
    
    if (!plot){
      return(list(df.pca=df.pca,df.pca.sp=df.pca.sp))
    }

    plot <- ggplot(data = df.pca,
                   aes(x = df.pca[,ax1], y = df.pca[,ax2])) # Plot ax1 and ax2

    # Draw the site scores
    if (plot.sites) {
      plot <- plot + geom_point(aes(col = color.vec,
                                    shape = shape.vec,
                                    fill = fill.vec),
                                size = 4,
                                ...) # Draw the site scores
      if (label.sites)
        plot <- plot + geom_text_repel(aes(label = rownames(df.pca)),
                                       size = 3)
    } else {
      if (label.sites) {
        plot <- plot + geom_text_repel(aes(label = rownames(df.pca)),
                                       size = 3)
      }
    }

    # Draw the species scores
    if (plot.spe) {
      plot <- plot + geom_segment(data = df.pca.sp, # Draw the species scores
                                  aes(x = 0, y = 0, xend = df.pca.sp[,ax1], yend = df.pca.sp[,ax2]),
                                  arrow=arrow(length=unit(0.01,"npc")),
                                  color = "#C20000")
      if (label.spe)
        plot <- plot + geom_text_repel(data = df.pca.sp,
                                       aes(x = df.pca.sp[,ax1], y = df.pca.sp[,ax2], label = spe.names(rownames(df.pca.sp))),
                                       color = "#C20000",
                                       size = 4,
                                       fontface = "italic")
    } else {
      if (label.spe)
        plot <- plot + geom_text_repel(data = df.pca.sp,
                                       aes(x = df.pca.sp[,ax1], y = df.pca.sp[,ax2], label = spe.names(rownames(df.pca.sp))),
                                       color = "#C20000",
                                       size = 4,
                                       fontface = "italic")
    }

    # If scaling = 1 draw circle of equilibrium contribution
    if (scaling == 1 & draw.circle) {
      plot <- plot +
        geom_path(data = pca.circle(pca.out, mult.spe, fact.spe),
                  aes(x, y),
                  color = "#C20000")
    }

    plot <- plot +
      geom_hline(yintercept=0, linetype="dotted") +
      geom_vline(xintercept=0, linetype="dotted") +
      labs(title = paste("PCA Scaling", scaling),
           x = paste0("PC", ax1, " (",
                      round(100 * eig.val.rel[ax1], 1), "%)"),
           y = paste0("PC", ax2, " (",
                      round(100 * eig.val.rel[ax2], 1), "%)")) +
      theme_linedraw() +
      theme(panel.grid = element_blank(),
            legend.background = element_blank(),
            legend.box.background= element_rect(colour="black"),
            legend.text = element_text(size = 14),
            legend.title = element_text(size = 16),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18))

    # Or with geom_circle
    # require(ggforce)
    # geom_circle(data = NULL,
    #             inherit.aes = FALSE,
    #             aes(x0 = 0, y0 = 0, r = radius),
    #             n= 1000)

    #plot + geom_path(data = data.frame(x = (0 + radius *cos(seq(0,2*pi,length.out = 100))), y = (0 + radius *sin(seq(0,2*pi,length.out = 100)))), aes(x, y))

    plot
}

