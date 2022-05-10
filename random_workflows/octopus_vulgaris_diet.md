library(tidyverse)
library(reshape2)
library(viridis)


data <- read.table("species_counts.txt", header=T, sep="\t")
data_2 <- melt(data, id.vars = c("Species_group", "Species", "Count", "Reads", "Reads_percent"))

ggplot(data_2) + geom_tile(aes(variable,Species,fill=log10(value))) + scale_fill_viridis() + scale_y_discrete(limits = rev)



# no facet
ggplot(data_2) + geom_tile(aes(variable,paste0(Species_group," | ",Species),fill=log10(value+1), group=Species_group)) + scale_fill_viridis() + scale_y_discrete(limits = rev)


# with a facet
ggplot(data_2) +
     geom_tile(aes(variable, Species, fill=log10(value+1), group=Species_group)) + scale_fill_viridis(option="plasma") +
     scale_y_discrete(limits = rev) +
     facet_grid(Species_group~., scales="free", space = "free") + theme(panel.spacing.y = unit(0, "lines"), strip.text.y = element_text(angle = 0)) +
     labs(y="Prey species", x="Sampling location and conditions", fill="log10(read_count)")



library(tidyverse)
library(reshape2)
library(viridis)
library(patchwork)

data <- read.table("species_counts_summary.txt", header=T, sep="\t")
data_2 <- melt(data, id.vars = c("Group"))


# plot
plot_1 <- ggplot(data_2, aes(x=variable, y=value, fill=Group)) +
     geom_bar(position="fill", stat="identity") +
     scale_fill_viridis_d(option="plasma")

# log10+1 transformed
plot_2 <- ggplot(data_2, aes(x=variable, y=log10(value+1), fill=Group)) +
     geom_bar(position="fill", stat="identity") +
     scale_fill_viridis_d(option="plasma")

plot_1 + plot_2 + plot_layout(ncol=1, guides = 'collect')
