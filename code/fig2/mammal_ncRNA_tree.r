#初期化
rm(list = ls())

#ライブラリの読み込み
library(ape)
library(ggtree)
library(tidyverse)
library(parallel)
library(RColorBrewer)
library(RMySQL)
library(DBI)
library(ConfigParser)

#シード値の設定
set.seed(1234)

tree_function <- function(x) nj(dist(x)) #ブートストラップ計算関数
nboot <- 1000 #試行回数
organism_list <- c("brain","heart","kidney","liver","testes")
#色の指定
cols <- brewer.pal(6, "Pastel1")
color_list <- list(
    "human  " = cols[1],
    "chimp  " = cols[2],
    "macaque  " = cols[3],
    "cow  " = cols[4],
    "mouse  " = cols[5],
    "rat  " = cols[6],
    "OK" = "blue",
    "DN" = "green",
    "NG" = "red"
)

#db接続
password_config <- read.ini('password.ini') #パスワードファイルの読み込み
password <- password_config$development$password
con <- dbConnect(
        user = 'nakashima',
        MySQL(),
        password = password,
        dbname = 'nakashima_db',
        host = 'localhost',
        port = 3306
)

#コマンド作成
command <- readLines('data/sql/fig2.sql')
command <- command[!grepl("^-",command)]
command <-  paste(command,collapse = "")

#データの取得&前処理
df_raw <- as_tibble(dbGetQuery(con,command)) %>%
    pivot_longer(
        cols = -c("id","species"),
        names_to = "Organism"
        )

dbDisconnect(con) #dbとの接続解除

for (organism in organism_list) {
    df_organism <- df_raw %>% filter(Organism == organism) %>%
    select(-Organism) %>% pivot_wider(names_from = id,values_from = value)
    df_index <- df_organism %>% select(species) %>%
        mutate(
            species = case_when(
                species == "hg19" ~ "human  ",
                species == "panTro3" ~ "chimp  ",
                species == "rheMac2" ~ "macaque  ",
                species == "bosTau6" ~ "cow  ",
                species == "mm9" ~ "mouse  ",
                species == "rn4" ~ "rat  ",
                TRUE ~ "unknown  "
            )
        )
    df_value <- df_organism %>% select(-c(species))

    #ブート無しのツリー作成
    tree <- tree_function(df_value)
    #root指定
    tree <- root(
        tree,as.character(which(df_index == "cow  ")),
        resolve.root = T
        )

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
    tree$tip.label <- df_index$species
    tree$node.label <- boot_values$value
    #枝長の正規化
    branch_sum <- sum(tree$edge.length)
    tree$edge.length <- tree$edge.length/branch_sum
    #描写範囲取得
    tree_limit <- plot(tree)$x.lim #X11が起動しちゃうかも
    #プレゼン用にnewick形式として保存
    write.tree(
        tree,
        paste0("output/fig2/nwk/mammal_ncRNA_",organism,".nwk")
    )

    #ラベル背景の色指定用dataframe
    color_df <- tibble(
        node=1:(Nnode(tree) + Ntip(tree)),
        #ノードの行はbootの値のstatusで埋める
        color = c(df_index$species,boot_values$status)
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
            nudge_y = tree_limit[2] * 0.4) +
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
        paste0("output/fig2/tree/mammal_ncRNA_",organism,".pdf"),
        plot = g,width = 7,width = 7
        )
}

