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

#normailze関数
normalize_func <- function(x) {
        min_val <- min(x)
        max_val <- max(x)

        return((x - min_val)/(max_val - min_val))
}

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
idx_command <- readLines('data/sql/fig3_idx.sql')
idx_command <- idx_command[!grepl("^-",idx_command)]
idx_command <-  paste(idx_command,collapse = "")

df_command <- readLines('data/sql/fig3_df.sql')
df_command <- df_command[!grepl("^-",df_command)]
df_command <-  paste(df_command,collapse = "")

#発現データ取得＆前処理
df_raw <- as_tibble(dbGetQuery(con,df_command)) %>% 
    distinct() %>% group_by(`Gene ID`,Run) %>%
    summarize(TPM = sum(TPM),.groups = 'drop') %>%
    pivot_wider(names_from = `Gene ID`,values_from = TPM)

#indexデータ取得&前処理
idx_raw <- as_tibble(dbGetQuery(con,idx_command)) %>%
    mutate(source_name = case_when(
        source_name ==
        "abdomen without digestive or reproductive system"
        ~ "abdomen",
        source_name ==
        "digestive plus excretory system"
        ~ "digestive",
        source_name ==
        "reproductive system without gonad and genitalia; reproductive system without gonad"
        ~"reproductive system",
        source_name ==
        "reproductive system without gonad and genitalia"
        ~ "reproductive system",
        source_name ==
        "thorax without digestive system"
        ~ "thorax",
        source_name ==
        "3rd instar larvae antennal disc"
        ~ "antenna",
        source_name ==
        "8hr APF pupal antennal disc"
        ~ "antenna",
        source_name ==
        "40hr APF pupal antenna"
        ~ "antenna",
        source_name ==
        "reproductive system without gonad"
        ~ "reproductive system",
        source_name ==
        "adult antenna"
        ~ "antenna",
        source_name ==
        "testis" 
        ~ "gonad",
        TRUE ~ source_name
        )
        ) %>%
        mutate(Organism = case_when(
        Organism == "Drosophila willistoni" ~ "D.wil",
        Organism == "Drosophila ananassae" ~ "D.ana",
        Organism == "Drosophila sechellia" ~"D.sec",
        Organism == "Drosophila simulans" ~ "D.sim",
        Organism == "Drosophila melanogaster" ~ "D.mel",
        Organism == "Drosophila yakuba" ~ "D.yak",
        Organism == "Drosophila pseudoobscura" ~ "D.pse",
        TRUE ~ Organism
        )
        #使うデータだけ取る
        ) %>% filter(source_name %in% organism_list)

df <- inner_join(df_raw,idx_raw,by = "Run") %>% select(-Run)

#dbとの接続解除&メモリの開放
dbDisconnect(con)
rm(idx_raw,df_raw)

for (organism in organism_list){
    #Organismごとにデータの切り出し
    df_organism <- df %>% filter(source_name == organism) %>%
        select(-source_name) %>% select(where(~all(!is.na(.))))
    
    #各分類群ごとに遺伝子発現量の中央値を算出
    df_group_median <- df_organism %>% group_by(Organism) %>%
    summarise_all(median)

    #分類群の名前
    df_index <- df_group_median %>% select(Organism) %>%
        mutate(Organism = paste0(Organism, "  ")) #ラベルの位置調整
    #発現量
    df_value <- df_group_median %>% select(-Organism)

    #ブート無しのツリー作成
    tree <- tree_function(df_value)
    #root指定
    tree <- root(
        tree,as.character(which(df_index == "D.wil  ")),
        resolve.root = TRUE)

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
        paste0("output/fig3/nwk/drosophila_mRNA_",organism,".nwk")
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

    ggsave(
        paste0("output/fig3/tree/drosophila_mRNA_",organism,".pdf"),
        plot = g,width = 7,height = 7
        )
}
