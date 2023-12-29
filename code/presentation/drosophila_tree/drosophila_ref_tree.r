#ライブラリの読み込み
library(ape)
library(ggtree)
library(tidyverse)
library(RColorBrewer)

#newick形式ファイルの読み込み
tree <- read.tree("data/nwk/drosophia_ref.nwk")
#描写範囲取得
tree_limit <- plot(tree)$x.lim
#ラベルのズレを調整
tree$tip.label <- c(
    "D.wil     ",
    "D.ana     ",
    "D.sec     ",
    "D.mel     ",
    "D.yak     ",
    "D.pse     "
)

#色の指定
cols <- brewer.pal(6, "Pastel1")
color_list <- list(
    "wil" = cols[1],
    "ana" = cols[2],
    "sec" = cols[3],
    "mel" = cols[4],
    "yak" = cols[5],
    "pse" = cols[6]
)

#phyloオブジェクトに色指定の列を追加
color_df <- tibble(
    node=1:(Nnode(tree) + Ntip(tree)),
    #ノードの行はnonで埋める
    color = c("wil","ana","sec","mel","yak","pse",c(rep("non", Nnode(tree))))
    )
tree <- full_join(tree, color_df, by="node")

#描写
g <- ggtree(tree) +
    geom_tiplab(
        #色指定
        aes(fill = color),
        size = 18,
        geom = "label",
        #余白を大きく
        label.padding = unit(0.8, "lines"),
        #枠線無くす
        label.size = 0,
        #角を丸くしない
        label.r =  unit(0, "lines"),
        show.legend = FALSE) +
    scale_x_continuous(limits = c(0,tree_limit[2] * 1.4)) +
    scale_fill_manual(values = color_list)

#保存
ggsave("output/presentation/drosophila_tree/drosophila_ref.pdf",plot = g)
