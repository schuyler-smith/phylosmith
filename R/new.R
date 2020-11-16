# phyloseq_obj = soil_column
# classification = 'Phylum'
# treatment = c('Matrix', 'Treatment')
# subset = NULL
# transformation = 'none'
# colors = 'default'
# wrap_by = NULL
# sep = '_'
# bar_width = 0.52
#
# # taxa_abundance_bars <-
# #   function(phyloseq_obj,
# #            classification = NULL,
# #            treatment,
# #            subset = NULL,
# #            transformation = 'none',
# #            colors = 'default',
# #            wrap_by = NULL) {
#
#     phyloseq_obj <-
#       taxa_filter(phyloseq_obj,
#                   treatment,
#                   frequency = 0,
#                   subset = subset)
#     if (!(is.null(classification))) {
#       phyloseq_obj <- conglomerate_taxa(phyloseq_obj,
#                                         classification, hierarchical = TRUE)
#     } else {
#       classification <- 'OTU'
#     }
#     treatment_name <- paste(treatment, collapse = sep)
#     if(!(is.null(wrap_by)) && !(wrap_by %in% treatment)){
#       treatment <- c(treatment, wrap_by)
#     }
#
#     graph_data <- melt_phyloseq(phyloseq_obj)
#     graph_data[[classification]] <-
#       factor(graph_data[[classification]],
#              levels = unique(graph_data[[classification]]))
#
#     set(graph_data, which(is.na(graph_data[[classification]])),
#         classification, 'Unclassified')
#
#     treatments <- levels(graph_data[[treatment_name]])
#     color_count <- length(treatments)
#     graph_colors <- create_palette(color_count, colors)
#
#     means <- graph_data[, mean(Abundance), by = c(treatment_name, treatment, classification)]
#     groups <- unique(means[V1 > 0.01*max(V1)][[classification]])
#     graph_data <- graph_data[get(classification) %in% groups]
#     g <-
#       ggplot(graph_data,
#              aes_string(x = classification, y = abundance, fill = treatment_name))
#     g <- g + stat_summary(geom = "bar",
#                      fun.y = eval(transformation),
#                      position = position_dodge2(padding = 3.5),
#                      size = 0.45,
#                      color = 'black',
#                      alpha = 0.85,
#                      width = bar_width) +
#       guides(colour = guide_legend(ncol = ceiling(length(
#         unique(graph_data[[classification]])
#       ) / 25))) +
#       scale_fill_manual(values = graph_colors,
#                         aesthetics = c('color', 'fill'))
#     g <- g + stat_summary(geom = "errorbar",
#                      fun.data = mean_se,
#                      position = "dodge",
#                      size = 0.4,
#                      width = 0.52)
#     g <- g + theme_light() +
#       theme(
#         axis.line.x = element_line(
#           colour = 'black',
#           size = 1,
#           linetype = 'solid'
#         ),
#         axis.line.y = element_line(
#           colour = 'black',
#           size = 1,
#           linetype = 'solid'
#         ),
#         axis.text.x = element_text(
#           size = 10,
#           vjust = 1,
#           hjust = 1,
#           angle = 30
#         ),
#         axis.text.y = element_text(size = 10),
#         axis.title.x = element_text(size = 10, face = "bold"),
#         axis.title.y = element_text(size = 10, face = "bold"),
#         legend.text = element_text(size = 8),
#         legend.title = element_text(size = 10, face = "bold"),
#         legend.background = element_rect(fill = (alpha = 0)),
#         legend.key.size = unit(4, "mm"),
#         legend.spacing.x = unit(0.005, 'npc'),
#         strip.text.x = element_text(size = 10, face = 'bold', color = 'black'),
#         strip.background = element_rect(colour = 'black', size = 1.4, fill = 'white')
#       ) +
#       scale_y_continuous(expand = expand_scale(mult = c(0.0025, 0.002)))
#
#
#
#     buffer <- ggplot_build(g)$layout$panel_scales_y[[1]]$range$range[2] * 0.05
#     ttest <- setNames(data.table(matrix(ncol = 5, nrow = 0)), c(classification, "p_value", "start", "end", "y"))
#     for (taxon in  groups) {
#       sample <-  graph_data[get(classification) == taxon]
#       t <- pairwise.wilcox.test(sample$Abundance, sample[[treatment_name]])[['p.value']]
#       for(r in rownames(t)){
#         for(c in colnames(t)){
#           ttest <- rbind(ttest,
#                          setNames(data.frame(taxon, t[r,c],
#                                              which(groups == taxon), ## this is going to set it at whichever taxon it is at.. then adjust to the 0.52 bar width based on which treatment it is!!!!
#                            c,
#                            max(means[get(classification) == taxon & get(treatment_name) %in% c(r,c)]$V1) + buffer),
#                            c(classification, "p_value", "start", "end", "y")))
#         }
#       }
#     }
#
#     ttest[, significance := cut(p_value,  breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "NS"))]
#     ttest <- ttest[significance %in% c("***", "**", "*")]
#     g + ggsignif::geom_signif(data = ttest,
#                               aes(xmin = start,
#                                   xmax = end,
#                                   y_position = y,
#                                   annotations = significance),
#                               step_increase = 0.12,
#                               vjust = 1,
#                               size = 0.2,
#                               tip_length = 0.01,
#                               manual = TRUE)
#
#
#
#
#     if(!is.null(wrap_by)){
#       g <- g + facet_wrap(reformulate(wrap_by))
#     }
#     return(g)
