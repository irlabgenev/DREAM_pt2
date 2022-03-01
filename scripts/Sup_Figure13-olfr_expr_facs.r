require(tidyverse)
require(patchwork)

setwd("D:/Projects/dream2-figures/FigureS13/")
source("N:/Novell/__articles in progress__/DREAM transcriptomic adaptations/figures/R scripts/common_layout.r")

or_list <- c(ape::read.tree("Mus_musculus.GRCm38.olfrs.clustalo_full.trim.phyml.rooted.20210916.txt")$tip.label,
             "Olfr151")

d.list <- mapply(d =read_rds("FACSseq_M71_20170405_MOR23_20170112_gene_counts_per_sample-remapping20210310.rds"),
                 MoreArgs = list("or_list" = or_list),
                 SIMPLIFY = F,
                 FUN = function(d, or_list) {
                   d <- rownames_to_column(d, "gene_name") %>%
                     pivot_longer(cols = where(is.numeric), 
                                  values_to = "counts", 
                                  names_to = "condition") %>%
                     group_by(condition) %>%
                     mutate(cpm = counts / (sum(counts) / 1e+06)) %>%
                     filter(gene_name %in% or_list) %>%
                     group_by(gene_name) %>%
                     summarise(cpm.se = sd(cpm)/sqrt(length(cpm)),
                               cpm.mean = mean(cpm)) %>%
                     arrange(-cpm.mean) %>%
                     mutate(rank = 1:n(),
                            gene_name.label = factor(gene_name, levels = unique(gene_name)))
                   return(d)
                   })


# bar plot
p.list <- mapply(d = d.list,
       condition = str_extract(names(d.list), "(?<=_)[^_]+(?=_)"),
       ymax = lapply(d.list, function(d) { max(d$cpm.mean) + max(d$cpm.mean)*0.2 }),
       SIMPLIFY = F,
       FUN = function(d, condition, ymax) {
         p <- ggplot(d, aes(gene_name.label, cpm.mean)) +
           geom_errorbar(aes(ymin = cpm.mean,
                             ymax = cpm.mean + cpm.se),
                         width = 0.2) +
           geom_col() +
           scale_y_continuous(limits = c(0, ymax), expand = c(0,0)) +
           scale_x_discrete(limits = filter(d, rank <= 40)$gene_name, expand = c(1/40, 1/40)) +
           common_layout +
           theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
           ggtitle(condition) +
           xlab("top 40 most expressed Olfrs") +
           ylab("mean CPM")
          
       })

# dot plot
p.list <- mapply(d = d.list,
                 condition = str_extract(names(d.list), "(?<=_)[^_]+(?=_)"),
                 ymax = lapply(d.list, function(d) { max(d$cpm.mean) + max(d$cpm.mean)*0.2 }),
                 SIMPLIFY = F,
                 FUN = function(d, condition, ymax) {
                   p <- ggplot(d, aes(gene_name.label, cpm.mean)) +
                     geom_point() +
                     geom_segment(aes(xend = gene_name.label,
                                      y = cpm.mean - cpm.se,
                                      yend = cpm.mean + cpm.se),
                                   width = 0.2) +
                     scale_y_continuous(limits = c(0, ymax)) +
                     scale_x_discrete(limits = filter(d, rank <= 40)$gene_name, expand = c(1/40, 1/40)) +
                     common_layout +
                     theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                     ggtitle(condition) +
                     xlab("top 40 most expressed Olfrs") +
                     ylab("mean CPM")
                   
                 })

p <- wrap_plots(p.list, ncol = 1)

ggsave(filename = "fig_supp9-olfr_expr_facs-dots.pdf",
       plot = p,
       device = "pdf", 
       units = "cm",
       width = 20, 
       height = 20, 
       useDingbats=FALSE)

