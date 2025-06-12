# This code helps exporting graphs to images

# Exporting data
setwd("/Users/ducciorocchini/Downloads")
setwd("~/Desktop")
# Windowds users: C://comp/Downloads
# \
# setwd("C://nome/Downloads")

getwd()

plot(gr)

png("greenland_output.png")
plot(gr)
dev.off()

pdf("greenland_output.pdf")
plot(gr)
dev.off()

pdf("difgreen.pdf")
plot(grdif)
dev.off()

jpeg("difgreen.jpeg")
plot(grdif)
dev.off()

# anche ggsave per salvare pdf come nel corso di R statistica dove possibile
# usa anche im.export per mostrare che esporti i file .tif
