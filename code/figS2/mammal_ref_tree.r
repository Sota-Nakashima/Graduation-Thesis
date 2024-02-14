#初期化
rm(list = ls())

#パッケージのインストール
library(ape)
library(ggtree)
library(tidyverse)
library(RColorBrewer)

#画像のディレクトリ
image_dir <- "data/picture/"

#newick形式ファイルの読み込み
tree <- read.tree("data/nwk/mammal_ref.nwk")
#描写範囲取得
tree_limit <- plot(tree)$x.lim
#ラベルのズレを調整
tree$tip.label <- c(
    "cow  ",
    "rat  ",
    "mouse  ",
    "human  ",
    "chimp  ",
    "macaque  "
)

#色の指定
cols <- brewer.pal(6, "Pastel1")
color_list <- list(
    "human" = cols[1],
    "chimp" = cols[2],
    "macaque" = cols[3],
    "cow" = cols[4],
    "mouse" = cols[5],
    "rat" = cols[6]
)

#phyloオブジェクトに色指定の列を追加
color_df <- tibble(
    node=1:(Nnode(tree) + Ntip(tree)),
    #ノードの行はnonで埋める
    color = c("cow","rat","mouse","human","chimp","macaque",c(rep("non", Nnode(tree))))
    )
tree <- full_join(tree, color_df, by="node")

#描写
g <- ggtree(tree) +
    geom_tiplab(
        #色指定
        aes(fill = color),
        size = 12,
        geom = "label",
        #余白を大きく
        label.padding = unit(0.5, "lines"),
        #枠線無くす
        label.size = 0,
        #角を丸くしない
        label.r =  unit(0, "lines"),
        show.legend = FALSE) +
    scale_x_continuous(limits = c(0,tree_limit[2] * 1.1)) +
    scale_fill_manual(values = color_list)

#保存
ggsave(
    "output/figS2/mammal_ref.pdf",plot = g,
    width = 10.1,height = 10.1)
