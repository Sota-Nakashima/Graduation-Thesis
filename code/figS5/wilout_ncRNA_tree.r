#ライブラリの読み込み
library(ape)
library(ggtree)
library(tidyverse)
library(parallel)
library(RColorBrewer)

#シード値の設定
set.seed(1234)

tree_function <- function(x) nj(dist(x)) #ブートストラップ計算関数
nboot <- 1000 #試行回数
organism_list <- c("abdomen","digestive","gonad","head","thorax")
#色の指定
cols <- brewer.pal(6, "Pastel1")
color_list <- list(
    "D.wil  " = cols[1],
    "D.ana  " = cols[2],
    "D.sim  " = cols[3],
    "D.mel  " = cols[4],
    "D.yak  " = cols[5],
    "D.pse  " = cols[6],
    "OK" = "blue",
    "DN" = "green",
    "NG" = "red"
)

#データの読み込み
#あとでMySQL仕様に変更
df <- read_csv("data/tmp/ncRNA.csv")
df <- df %>% filter(Organism != "D.wil") %>% #wilを除く
    filter(source_name %in% organism_list) %>%
    select(where(~all(!is.na(.))))

for (organism in organism_list){
    #Organismごとにデータの切り出し
    df_organism <- df %>% filter(source_name == organism)
    #各分類群ごとに遺伝子発現量の中央値を算出
    df_group_median <- df_organism %>% select(-c(Run,source_name,sex)) %>%
        group_by(Organism) %>% summarise_all(median)
    ###########
    #ここに正規化処理を書く
    ###########
    #分類群の名前
    df_index <- df_group_median %>% select(Organism) %>%
        mutate(Organism = paste0(Organism, "  ")) #ラベルの位置調整
    #発現量
    df_value <- df_group_median %>% select(-Organism)

    #ブート無しのツリー作成
    tree <- tree_function(df_value)
    #root指定
    tree <- root(
        tree,as.character(which(df_index == "D.pse  ")),
        resolve.root = T)

    #ブートの計算
    bs_tree <- boot.phylo(
        tree,
        df_value,
        tree_function,
        B = nboot, #試行回数
        block = nrow(df_value) - 1, #ジャックナイフサンプリング
        jumble = FALSE, #ジャックナイフサンプリング
        mc.cores = detectCores()/2, #コア数
        trees = TRUE
    )$trees
    #boot値を計算
    boot_values <- as_tibble(round(prop.clades(tree, bs_tree)/(nboot / 100)))
    #boot値によって違う印をつける(treeの描写には関係無し)
    boot_values <- boot_values %>%
        mutate(
            status = case_when(
                value >= 90 ~ "OK",
                value >= 70 ~ "DN",
                TRUE ~ "NG"
            )
        )

    #ラベルを代入
    tree$tip.label <- df_index$Organism
    tree$node.label <- boot_values$value
    #枝長の正規化
    branch_sum <- sum(tree$edge.length)
    tree$edge.length <- tree$edge.length/branch_sum
    #描写範囲取得
    tree_limit <- plot(tree)$x.lim
    #プレゼン用にnewick形式として保存
    write.tree(
        tree,
        paste0("output/figS5/nwk/wilout_ncRNA_",organism,".nwk")
    )

    #ラベル背景の色指定用dataframe
    color_df <- tibble(
        node=1:(Nnode(tree) + Ntip(tree)),
        #ノードの行はbootの値のstatusで埋める
        color = c(df_index$Organism,boot_values$status)
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
            show.legend = FALSE) +
        geom_nodelab(
            size = 7,
            #数字の位置調整
            nudge_x = -tree_limit[2] * 0.04,
            nudge_y = tree_limit[2] * 0.2) +
        scale_x_continuous(
            limits = c(-tree_limit[2] * 0.05,tree_limit[2] * 1.1)
            ) + #端まで描写
        scale_fill_manual(values = color_list) + #色の指定
        annotate(
            "text",x = -Inf,y = Inf,label = organism,
            hjust = -.2,vjust = 2,size = 9
            ) #器官のタイトル
    #保存
    ggsave(
        paste0("output/figS5/tree/wilout_ncRNA_",organism,".pdf"),
        plot = g
        )

}