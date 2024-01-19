#初期化
rm(list = ls())

#ライブラリの読み込み
library(tidyverse)
library(RMySQL)
library(DBI)
library(ConfigParser)
library(scatterplot3d)
library(animation)

#シード値の設定
set.seed(1234)

#色の指定
colors <- scale_color_hue()$palette(9)

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
idx_command <- readLines('data/sql/presentation/drosophila_PCA_idx.sql')
idx_command <- idx_command[!grepl("^-",idx_command)]
idx_command <-  paste(idx_command,collapse = "")

df_command <- readLines('data/sql/presentation/drosophila_PCA_df.sql')
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
        ))

#データの結合&前処理
df <- inner_join(df_raw,idx_raw,by = "Run") %>% #select(-Run) %>%
    select(where(~all(!is.na(.)))) #欠損値を除く

#dbとの接続解除&メモリの開放
dbDisconnect(con)
rm(idx_raw,df_raw)

#データの抽出
df_index <- df %>% select(c(Organism,source_name,Run))
df_value <- df %>% select(-c(Organism,source_name,Run)) %>%
    apply(1,normalize_func) %>% t() %>% as_tibble() %>% #正規化
    mutate_all(~ifelse(is.na(.), 0, .)) #NAを0にする

#pcaの計算
pca_result <- prcomp(df_value)

#後処理
pca_df <- as_tibble(pca_result$x) %>% #結果の抽出
    mutate_all(~sign(.) * log1p(abs(.))) %>% #neg-log変換
    bind_cols(df_index) %>% #index列の結合
    mutate(
        colors = case_when(
            source_name == "antenna" ~ colors[1],
            source_name == "gonad" ~ colors[2],
            source_name == "abdomen" ~ colors[3],
            source_name == "digestive" ~ colors[4],
            source_name == "head" ~colors[5],
            source_name == "reproductive system" ~ colors[6],
            source_name == "thorax" ~ colors[7],
            source_name == "whole body" ~ colors[8],
            source_name == "genitalia" ~ colors[9],
            TRUE ~ source_name
        )
    ) #色のために番号分け

#scatterplot3d(
#    pca_df[,1:3],
#    pch = 16,
#    color = pca_df$colors,
#    angle = 30
#)
plot3d360 <- function() {
    for (i in 1:360)  {
        scatterplot3d(
            pca_df[,1:3],
            pch = 16,
            color = pca_df$colors,
            main = i,
            angle = i)
    }
}

saveGIF(
    plot3d360(),
    movie.name = "output/presentation/animation_mRNA.gif",
    interval=0.05,movietype="gif")
