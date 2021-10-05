
# common layout
common_layout <- theme_classic(base_size = 10, base_family = "ArialMT") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    legend.key.size = unit(0.7, "lines")
  )

clusters.fill.moe <- c("GBC" = "#FCB519", 
                       "INP" = "#ED7953", 
                       "iOSN1" = "#EA71B7",     
                       "iOSN2" = "#4F99D3",      
                       "mOSN" = "#CA404A",  
                       "Sus"="#47C16E",       
                       "Mv2"="#24868E",       
                       "HBC1"="#3B528B",   
                       "HBC2"="#0D0887")

clusters.fill.osn <- c("Ventral"="#47C16E",
                        "Calb2+ V."="#24868E", 
                        "Dlg2+ V."="#4F99D3",    
                        "Cd36+ V."="#0D0887",     
                        "Cd55+ V."="#481F70",  
                        "Dorsal"="#FCA108",       
                        "Calb2+ D."="#ED7953",    
                        "Dlg2+ D."="#EA71B7",   
                        "Cd36+ D."="#932B80",     
                        "Cd55+ D."="#CA404A")

dist.colors <- c("dist.S100a5" = "skyblue1",
                 "dist.Pcp4l1" = "deeppink1",
                 "dist.trans" = "yellow",
                 "dist.aa.tm" = "red",
                 "dist.aa.trimmed" = "green")

phylo.annotation <- c("P-element controlled genes (strong)" = "#DA4850",
                      "P-element controlled genes (moderate)" = "#FFBB74",
                      "P-cluster genes" = "#FFF5AF",
                      "9(chr7)" = "#9CD694")

### FUNCTIONS
choose.cluster <- function(clusters) {
  counts <- table(clusters)
  counts.max <- subset(counts, counts == max(counts))
  return(sample(names(counts.max), 1))
}


