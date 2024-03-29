#初期化
rm(list = ls())

#パッケージのインストール
library(ape)
library(ggtree)
library(tidyverse)

#画像のディレクトリ
image_dir <- "data/picture/"

#newick形式ファイルの読み込み
tree <- read.tree("output/fig1/nwk/mammal_mRNA_brain.nwk")

tree$tip.label <- paste0(tree$tip.label, " ")
#描写範囲取得
tree_limit <- plot(tree)$x.lim
#phyloオブジェクトに画像指定の列を追加
image_df <- tibble(
    node=1:(Nnode(tree) + Ntip(tree)),
    #ノードの行はnonで埋める
    image_path = c(
        "human","chimp","mouse","rat",
        "macaque","cow",
        c(rep("non", Nnode(tree))))
    ) %>%
    #パス名に変換
    mutate(image_path = paste0(image_dir,image_path, ".jpeg"))
tree <- full_join(tree, image_df, by="node")

#論文用のツリーの描写
g <- ggtree(tree,size = 1.5) + 
    #写真を付加
    geom_tiplab(
    aes(image = image_path),
    geom = "image",
    size=.24,
    #画像の位置調整
    offset = .27,
    align = T) +
    #ブートストラップ値
    geom_nodelab(
    size = 7,
    nudge_x = -tree_limit[2] * 0.1,#数字の位置調整
    nudge_y = tree_limit[2] * 0.4) +
    #名前を付加
    geom_tiplab(
        size = 18,
        geom = "label",
        #余白を大きく
        label.padding = unit(0.5, "lines"),
        #枠線無くす
        label.size = 0,
        #角を丸くしない
        label.r =  unit(0, "lines"),
        show.legend = FALSE) +
    geom_highlight(
        node = 5,fill = "#FF5050",
        extend = .30
        ) + 
    geom_highlight(
        node = 5,fill = "#FF5050",
        extend = -.11
        ) +
    scale_x_continuous(
        limits = c(-tree_limit[2] * 0.11,tree_limit[2]*2.0)
        ) + #描画範囲調整
    annotate(
            "text",x = -Inf,y = Inf,label = "brain",
            hjust = -.2,vjust = 2,size = 14
        ) #器官のタイトル

#保存
ggsave(
    "output/presentation/mammal_picture_tree/mammal_mRNA_brain.pdf",
    plot = g,width = 10.1,height = 10.1
    )
