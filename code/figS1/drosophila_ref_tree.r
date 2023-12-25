#ライブラリの読み込み
library(ape)
library(ggtree)
library(tidyverse)

#newick形式ファイルの読み込み
tree <- read.tree("data/nwk/drosophia_ref.nwk")
#ラベルのズレを調整
tree$tip.label <- c(
    "D.wil  ",
    "D.ana  ",
    "D.sec  ",
    "D.mel  ",
    "D.yak  ",
    "D.pse  "
)

#色を指定
#ggplot標準の色直す
col_pallete <-  rainbow(6)

#指定する用のリスト
color_manual_df <- list(
    "wil" = col_pallete[1],
    "ana" = col_pallete[2],
    "sec" = col_pallete[3],
    "mel" = col_pallete[4],
    "yak" = col_pallete[5],
    "pse" = col_pallete[6]
)

#phyloオブジェクトに色指定の列を追加
color_df <- tibble(
    node=1:(Nnode(tree) + Ntip(tree)),
    #ノードの行はnonで埋めてる
    color = c("wil","ana","sec","mel","yak","pse",c(rep("non", Nnode(tree))))
    )
tree <- full_join(tree, color_df, by="node")

#描写
g <- ggtree(tree) +
    geom_tiplab(
        #色指定
        aes(fill = color),
        size = 9,
        geom = "label",
        #余白を大きく
        label.padding = unit(0.5, "lines"),
        #枠線無くす
        label.size = 0,
        #角を丸くしない
        label.r =  unit(0, "lines"),
        show.legend = FALSE) 
#端まで描写＆色の指定
g <- g + scale_x_continuous(limits = c(0,layer_data(g)$xend + 8)) +
    scale_fill_manual(values = color_manual_df)

#保存
ggsave("output/figS1/drosophila_ref.pdf",plot = g)
