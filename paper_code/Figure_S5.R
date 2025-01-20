#!/usr/bin/env Rscript

# Author: Robert Wang (Xing Lab)
# Date: 2025.01.01
# Supplementary Figure 5a

# Expression ratios for haplotypes carrying pathogenic variants based on patient TEQUILA-seq reads. The expression ratio 
# for a haplotype is calculated as the number of reads assigned to that haplotype divided by the total number of reads 
# assigned across both haplotypes. Variants located in the -3 to +6 region of the donor splice site or the -20 to +3 region
# of the acceptor splice site were classified as "splice site region" variants. 

# =====================================================================================================================
#                                                      LIBRARIES 
# =====================================================================================================================

# Load required libraries
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))

# =====================================================================================================================
#                                                        MAIN
# =====================================================================================================================

# Establish working directories and file paths
workdir <- "/mnt/isilon/lin_lab_share/STRIPE"
outfile <- file.path(workdir, "manuscript/Revisions/20241219/Supplementary_Figures/Figure_S5/Figure_S5.pdf")

# Construct dataframe summarizing all heterozygous pathogenic variants and their haplotype expression ratios
df <- tibble(
    Variant = c("chr16:Δ23525826-23551747", "NM_018127.7:c.297-2_297-1delinsT", 
        "NM_001278716.2:c.1838T>A, p.Val613Glu", "NM_001278716.2:c.1389+3_1389+6del",
        "NM_000143.4:c.40dup, p.Leu14ProfsTer42", "NM_004870.4:c.19G>T, p.Gly7Ter",
        "NM_018297.4:c.1533_1536del, p.Asn511LysfsTer51",
        "NM_025152.3:c.693+1G>A", "NM_176787.5:c.1674+1G>C",
        "NM_176787.5:c.1434+5G>A", "NM_176787.5:c.2126G>A, p.Arg709Gln",
        "NM_176787.5:c.1674+1G>C", "NM_176787.5:c.932T>G, p.Leu311Trp",
        "NM_176787.5:c.1291G>A, p.Gly431Arg", "NM_176787.5:c.963G>A, p.Gln321=",
        "NM_176787.5:c.1434+5G>A", "NM_176787.5:c.1759C>T, p.Arg587Ter",
        "NM_176787.5:c.2679C>G, p.Ser893Arg", "NM_176787.5:c.932T>G, p.Leu311Trp",
        "NM_004204.5:c.942+1G>A", "NM_004204.5:c.1199_1201del, p.Tyr400del",
        "NM_000303.3:c.415G>A, p.Glu139Lys", "NM_000303.3:c.422G>A, p.Arg141His",
        "NM_000303.3:c.357C>A, p.Phe119Leu", "NM_000303.3:c.422G>A, p.Arg141His",
        "chr16:Δ8805603-8820989", "NM_000303.3:c.338C>T, p.Pro113Leu",
        "NM_005660.3:c.211G>A, p.Val71Met", "NM_018389.5:c.503_505del, p.Phe168del",
        "NM_017646.6:c.22C>T, p.Arg8Ter", "chr22:Δ46334556-46336380",
        "NM_018006.5:c.1084G>A, p.Ala362Thr"),
    Class = factor(c("Structural deletion", "Splice site region", "Missense", "Splice site region",
        "Frameshift variant", "Nonsense", "Frameshift variant", "Splice site region", "Splice site region",
        "Splice site region", "Missense", "Splice site region", "Missense", "Missense", "Splice site region",
        "Splice site region", "Nonsense", "Missense", "Missense", "Splice site region", "Inframe deletion",
        "Missense", "Missense", "Missense", "Missense", "Structural deletion", "Missense", "Missense",
        "Inframe deletion", "Nonsense", "Structural deletion", "Missense")),
    Ratio = c(0.574, 0.318, 0.440, 0.560, 0.366, 0.458, 0.183, 0.655, 0.365, 0.635, 0.567, 0.416,
        0.584, 0.770, 0.230, 0.495, 0.187, 0.454, 0.546, 0.329, 0.671, 0.297, 0.703, 0.421, 0.579,
        0.502, 0.498, 0.474, 0.559, 0.782, 0.000, 1.000)
)
df <- df %>% group_by(Variant, Class) %>% summarise(Ratio = mean(Ratio))
df$Variant <- factor(df$Variant, levels = df %>% arrange(desc(Class), desc(Ratio)) %>% pull(Variant))

# Plot expression ratios for haplotypes carrying heterozygous pathogenic variants
p <- ggplot(df, aes(x = Ratio * 100, y = Variant, fill = Class)) + geom_bar(stat = "identity") + theme_classic() +
    geom_vline(xintercept = 50, linetype = "dashed", linewidth = 0.5) + theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text = element_text(color = "black", 
    size = 6), axis.title.x = element_text(color = "black", size = 7), axis.title.y = element_blank(), 
    axis.ticks = element_line(color = "black", linewidth = 0.25), axis.line = element_line(color = "black", linewidth = 0.25), 
    legend.title = element_blank(), legend.text = element_text(color = "black", size = 6), legend.key.size = unit(0.4, 'cm')) + 
    xlab("Haplotype expression ratio (%)") + scale_fill_manual(values = c("#CBB3D4", "#FCAD65", "#66C2A5", "#DC626E",
    "#388BBF", "#446C68"))

# Save p to outfile
ggsave(outfile, plot = p, width = 5, height = 4)