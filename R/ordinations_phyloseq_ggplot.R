#' Create a ggplot object of the NMDS from a phyloseq object. Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and plots the NMDS of a treatment or set of treatments.
#' @useDynLib phylosmith
#' @usage nmds_phyloseq_ggplot(phyloseq_obj, treatment, colors = 'Spectral', circle = TRUE)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package (must contain \code{\link[phyloseq:sample_data]{sample_data()}}).
#' @param treatment Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param colors Name of a color set from the \link[=RColorBrewer]{RColorBrewer} package.
#' @param circle If TRUE, add elipses around each treatment.
#' @import phyloseq
#' @import ggplot2
#' @import RColorBrewer
#' @import vegan
#' @export

nmds_phyloseq_ggplot <- function(phyloseq_obj, treatment, colors = "Spectral", circle = TRUE){
  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0)
  treatment <- paste(treatment, collapse = '.')
  getPalette = grDevices::colorRampPalette(brewer.pal(8, colors)); colorCount = 1 + length(unlist(unique(phyloseq_obj@sam_data[[treatment]]))); colors = getPalette(colorCount); theme_set(theme_bw())

  MDS <- metaMDS(t(phyloseq_obj@otu_table), autotransform = FALSE, distance = "bray", k=3, trymax=100)
  # plot(MDS, display = c("sites", "species"), choices = c(1,2), type = "p");abline(h=0,lty=2);abline(v=0,lty=2)
  # stressplot(MDS)

  Treatment <- phyloseq_obj@sam_data[[treatment]]
  NMDS1 <- data.table(scores(MDS))$NMDS1
  NMDS2 <- data.table(scores(MDS))$NMDS2
  NMDS <- data.table(NMDS1,NMDS2,Treatment)
  NMDS.narm <- subset(NMDS, !is.na(Treatment))

  p <- ggplot(data = NMDS.narm, aes(NMDS1, NMDS2, color = NMDS.narm$Treatment)) +
    # coord_fixed(xlim = c(floor(min(NMDS.narm[,c(1,2)])), ceiling(max(NMDS.narm[,c(1,2)]))),
    #             ylim = c(floor(min(NMDS.narm[,c(1,2)])), ceiling(max(NMDS.narm[,c(1,2)])))) +
    geom_point(aes(color = Treatment), size=1.5, alpha=0.75) +
    scale_color_manual(values=colors) +
    theme_classic() +
    theme(aspect.ratio=1,
          axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=1, linetype='solid'),
          axis.text.x=element_text(size=10, face = "bold"),
          axis.text.y=element_text(size=10, face = "bold"),
          axis.title.x=element_text(size=12, face= "bold"),
          axis.title.y=element_text(size=12, face= "bold"),
          legend.title=element_blank(),
          legend.text=element_text(size=11, face= "bold"),
          legend.background = element_rect(fill=(alpha = 0))
    )
  if(circle == TRUE){
    p <- p + stat_ellipse(geom="polygon", type="norm", size=.6, linetype = 1, alpha=0.0, aes(fill=ord$Treatment))
  }
  return(p)
}

#' Create a ggplot object using t-SNE from a phyloseq object. Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and plots the t-SNE of a treatment or set of treatments.
#' @useDynLib phylosmith
#' @usage tsne_phyloseq_ggplot(phyloseq_obj, treatment, perplexity = 10, colors = 'Spectral', circle = TRUE)
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \code{\link[=phyloseq]{phyloseq}} package (must contain \code{\link[phyloseq:sample_data]{sample_data()}}).
#' @param treatment Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param perplexity similar to selecting the number of neighbors to consider in decision making (should not be bigger than 3 * perplexity < nrow(X) - 1, see \code{\link[=Rtsne]{Rtsne}} for interpretation)
#' @param colors Name of a color set from the \code{\link[=RColorBrewer]{RColorBrewer}} package.
#' @param circle If TRUE, add elipses around each treatment.
#' @import phyloseq
#' @import ggplot2
#' @import RColorBrewer
#' @import vegan
#' @import Rtsne
#' @seealso \code{\link[=Rtsne]{Rtsne}}
#' @export

tsne_phyloseq_ggplot <- function (phyloseq_obj, treatment, perplexity = 10, colors = "Spectral", circle = TRUE){
  if (is.numeric(treatment)) {
    treatment <- colnames(phyloseq_obj@sam_data[, treatment])
  }
  phyloseq_obj <- taxa_filter(phyloseq_obj, treatment, frequency = 0)
  treatment <- paste(treatment, collapse = ".")
  getPalette = grDevices::colorRampPalette(brewer.pal(8, colors))
  colorCount = 1 + length(unlist(unique(phyloseq_obj@sam_data[[treatment]])))
  colors = getPalette(colorCount)

  tsne <- Rtsne(vegdist(t(phyloseq_obj@otu_table), method = 'bray'), dims = 2, theta = 0.0, perplexity=perplexity)

  Treatment <- phyloseq_obj@sam_data[[treatment]]
  tsne1 <- tsne$Y[,1]
  tsne2 <- tsne$Y[,2]
  ord <- data.table(tsne1, tsne2, Treatment)
  ord <- subset(ord, !is.na(Treatment))

  p <- ggplot(data = ord, aes(tsne1, tsne2, color = ord$Treatment)) +
    geom_point(aes(color = ord$Treatment), size=1.5, alpha=1) +
    scale_color_manual(values=colors) +
    theme_classic() +
    theme(aspect.ratio=1,
          axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=1, linetype='solid'),
          axis.text.x=element_text(size=10, face = "bold"),
          axis.text.y=element_text(size=10, face = "bold"),
          axis.title.x=element_text(size=12, face= "bold"),
          axis.title.y=element_text(size=12, face= "bold"),
          legend.title=element_blank(),
          legend.text=element_text(size=11, face= "bold"),
          legend.background = element_rect(fill=(alpha = 0))
    )
  if(circle == TRUE){
    p <- p + stat_ellipse(geom="polygon", type="norm", size=.6, linetype = 1, alpha=0.0, aes(fill=ord$Treatment))
  }
  return(p)
}



