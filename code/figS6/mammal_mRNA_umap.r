#初期化
rm(list = ls())

#ライブラリの読み込み
library(tidyverse)
library(umap)
library(RColorBrewer)
library(RMySQL)
library(DBI)
library(ConfigParser)

#シード値の設定
set.seed(1234)

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
command <- readLines('data/sql/figS6_mammal_mRNA.sql')
command <- command[!grepl("^-",command)]
command <-  paste(command,collapse = "")

#データの取得&前処理
df_raw <- as_tibble(dbGetQuery(con,command)) %>%
    pivot_longer(
        cols = -c("id","species"),
        names_to = "Organism"
        ) %>%
    pivot_wider(names_from = id,values_from = value) %>%
    select(where(~all(!is.na(.))))

dbDisconnect(con) #dbとの接続解除

df_index <- df_raw %>% select(Organism,species)
df_value <- df_raw %>% select(-Organism,-species) %>% select_if(~sum(.) != 0)

#計算
umap_result <- umap(df_value)

#結果を取り出してlegendフレームと結合
umap_df <- as_tibble(umap_result$layout) %>%
    bind_cols(df_index)

#描写
g <- ggplot(umap_df,aes(x = V1,y = V2,color = Organism)) +
    geom_point(size = 3) +
    stat_ellipse() + #確率楕円
    theme_classic() +
    labs(
        x = "UMAP1",
        y = "UMAP2",
        title = "Mammal",
        caption = "mRNA"
    ) +
    theme(
        text = element_text(size = 18),
        plot.title = element_text(hjust = 0.5)
    )
#保存
ggsave(
    "output/figS6/mammal_mRNA.pdf",plot = g,
    width = 7,height = 7)