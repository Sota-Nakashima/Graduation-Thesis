#パッケージのインストール
library(ape)
library(ggtree)
library(tidyverse)

#newick形式ファイルの読み込み
tree <- read.tree("data/nwk/pre_mammal_ref.nwk")
#描写範囲取得
tree_limit <- plot(tree)$x.lim

#phyloオブジェクトに画像指定の列を追加
image_df <- tibble(
    node=1:(Nnode(tree) + Ntip(tree)),
    #ノードの行はnonで埋める
    image_path = c(
        "human",
        "gorilla",
        "mouse",
        "opossum",
        "platypus",
        "chiken",
        "macaque",
        "orangutan",
        "bonobo",
        "chimp",
        c(rep("non", Nnode(tree))))
    ) %>%
    #パス名に変換
    mutate(image_path = paste0(image_dir,image_path, ".jpeg"))
tree <- full_join(tree, image_df, by="node")

#論文用のツリーの描写
#反対にしているのであとで反転処理
g <- ggtree(tree) + 
    #写真を付加
    geom_tiplab(
    aes(image = image_path),
    geom = "image",
    size=.14,
    #画像の位置調整
    offset = 31) +
    scale_x_continuous(limits = c(0,tree_limit[2])) #描画範囲調整

#保存
ggsave("output/presentation/mammal_picture_tree/mammal_preref.pdf",plot = g)
