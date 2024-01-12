#初期化
rm(list = ls())

#パッケージのインストール
library(ape)
library(ggtree)
library(tidyverse)

#画像のディレクトリ
image_dir <- "data/picture/"

#newick形式ファイルの読み込み
tree <- read.tree("output/fig1/nwk/mammal_mRNA_kidney.nwk")

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
g <- ggtree(
    tree,size = 1.5,branch.length = "none" #末端を揃える
    ) + 
    #写真を付加
    geom_tiplab(
    aes(image = image_path),
    geom = "image",
    size=.24,
    #画像の位置調整
    offset = .1,
    align = T) +
    scale_x_continuous(
        limits = c(0,4.5)
        ) + #描画範囲調整
    annotate(
            "text",x = -Inf,y = Inf,label = "kidney",
            hjust = -.2,vjust = 2,size = 9
        ) #器官のタイトル

#保存
ggsave(
    "output/presentation/mammal_picture_tree/mammal_mRNA_kidney.pdf",
    plot = g,width = 7,height = 7
    )