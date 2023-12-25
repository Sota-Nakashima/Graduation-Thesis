#パッケージのインストール
library(ape)
library(ggtree)
library(tidyverse)

#画像のディレクトリ
image_dir <- "data/picture/"

#newick形式ファイルの読み込み
tree <- read.tree("data/nwk/mammal_ref.nwk")

#phyloオブジェクトに画像指定の列を追加
image_df <- tibble(
    node=1:(Nnode(tree) + Ntip(tree)),
    #ノードの行はnonで埋めてる
    image_path = c(
        "cow","rat","mouse","human",
        "chimp","macaque",
        c(rep("non", Nnode(tree))))
    ) %>%
    #パス名に変換
    mutate(image_path = paste0(image_dir,image_path, ".jpeg"))
tree <- full_join(tree, image_df, by="node")


#論文用のツリーの描写
g <- ggtree(tree) + 
    #写真を付加
    geom_tiplab(
    aes(image = image_path),
    geom = "image",
    size=.24,
    #画像の位置調整
    offset = 65) +
    #名前を付加
    geom_tiplab(size = 9)
#描写位置調整
g <- g + scale_x_continuous(limits = c(0,layer_data(g)$xend + 80))

#保存
ggsave("output/figS1/mammal_ref.pdf",plot = g)
