library(ggh4x)

transfamlabel_df <- as.data.frame(trans_families_labels)
transfamlabel_df$trans_families_labels <- factor(transfamlabel_df$trans_families_labels, 
                                                 levels=c('Delta Method', 'GLM Residuals', 'Lat. Expr.', 'Count', 'Neg.'))
transfamlabel_df$trans_families_labels <- paste0("~ ", transfamlabel_df$trans_families_labels)
transfamlabel_df$family <- rownames(transfamlabel_df)
trans_families <-trans_families %>%
  left_join(transfamlabel_df)

res_main <- res %>%
  filter(alpha %in% c("TRUE", "FALSE")) %>%
  tidylog::inner_join(parameter_choices)  %>%
  mutate(knn_recovery = overlap / knn) %>%
  group_by(dataset, replicate, knn) %>%
  mutate(knn_recovery = knn_recovery / mean(knn_recovery)) %>%
  tidylog::left_join(trans_families)



# consistency_pl <- 
res_main %>%
  filter(benchmark == "consistency") %>%
  # make_main_plot_panel(add_group_label = FALSE) +
  ggplot(aes(x= knn_recovery, y = interaction(transformation, trans_families_labels, sep = "~ "), color = family)) +
  geom_vline(xintercept = 1, size = 0.3, linetype = 2) +
  ggbeeswarm::geom_quasirandom(color = "grey", size = 0.3, alpha = 0.7, groupOnX = FALSE) +
  stat_summary(geom = "point", position = position_dodge2(width = 0.3), fun.data = mean_cl_boot) +
  scale_y_discrete(guide = guide_axis_nested(delim = "~ ")) + 
  scale_x_continuous(breaks = c(0.5, 1, 1.5)) +
  coord_cartesian(xlim = c(0.2, 1.8)) +
  scale_color_manual(values = trans_families_colors, labels = trans_families_labels, guide = "none") +
  theme_grouped_axis(axis.grouping.line_padding = unit(5, "pt"), axis.grouping.line_height = unit(10, "pt"), 
                     axis.grouping.label.y = element_text(size = font_size_small, angle = 90)) +
  theme(axis.title.y = element_blank(),  plot.title.position = "plot",
        ggh4x.axis.nesttext.y = element_text(angle = 90, hjust = 0.5)) + 
  labs(x = "Relative $k$-NN Overlap",
       subtitle = "10X gene subset 1 vs. gene subset 2")
  
# TODO: order the plot 
# TODO: add some space between the grouping axis and original axis
# TODO: figure out LaTeX in plot axis



