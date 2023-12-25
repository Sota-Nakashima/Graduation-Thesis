#パッケージのインストール
library(ape)
library(ggtree)
library(tidyverse)

#newick形式ファイルの読み込み
tree <- read.tree("data/nwk/pre_mammal_ref.nwk")

#描写
g <- ggtree(tree) + geom_tiplab()
#描写位置調整
g <- g + scale_x_continuous(limits = c(0,layer_data(g)$xend + 40))

#保存
ggsave("output/presentation/mammal_ref_presentation.pdf",plot = g)
