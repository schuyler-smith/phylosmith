#' Create a ggplot object of the NMDS from a phyloseq object. Function from the phylosmith-package.
#'
#' This function takes a \code{\link[phyloseq]{phyloseq-class}} object and plots the NMDS of a treatment or set of treatments.
#' @useDynLib phylosmith
#' @usage nmds_phyloseq_ggplot(phyloseq_obj, treatment, colors = 'Spectral')
#' @param phyloseq_obj A \code{\link[phyloseq]{phyloseq-class}} object created with the \link[=phyloseq]{phyloseq} package (must contain \code{\link[phyloseq:sample_data]{sample_data()}}).
#' @param treatment Column name or number, or vector of, in the \code{\link[phyloseq:sample_data]{sample_data}}.
#' @param colors Name of a color set from the \link[=RColorBrewer]{RColorBrewer} package.
#' @import phyloseq
#' @import ggplot2
#' @import RColorBrewer
#' @import vegan
#' @export

nmds_phyloseq_ggplot <- function(phyloseq_obj, treatment, colors = "Spectral"){
  if(is.numeric(treatment)){treatment <- colnames(phyloseq_obj@sam_data[,treatment])}
  phyloseq_obj <- find_generalists(phyloseq_obj, treatment, frequency = 0)
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

  veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100){
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))}

  df_ell <- data.frame()
  for(trt in unique(NMDS.narm[[3]])){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS.narm[NMDS.narm[[3]]==trt,], veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))),group=trt))}

  p <- ggplot(data = NMDS.narm, aes(NMDS1, NMDS2)) +
    coord_fixed(xlim = c(floor(min(NMDS.narm[,c(1,2)])), ceiling(max(NMDS.narm[,c(1,2)]))),
                ylim = c(floor(min(NMDS.narm[,c(1,2)])), ceiling(max(NMDS.narm[,c(1,2)])))) +
    geom_point(aes(color = Treatment), size=1.5, alpha=0.75) +
    geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2, colour=group), size=1.5, linetype=1) +
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
  return(p)
}
