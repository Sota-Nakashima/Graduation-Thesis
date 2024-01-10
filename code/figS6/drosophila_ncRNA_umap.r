#初期化
rm(list = ls())

#ライブラリの読み込み
library(tidyverse)
library(RMySQL)
library(DBI)
library(ConfigParser)

#シード値の設定
set.seed(1234)

organism_list <- c("abdomen","digestive","gonad","head","thorax")

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
idx_command <- readLines('data/sql/figS6_drosophila_idx.sql')
idx_command <- idx_command[!grepl("^-",idx_command)]
idx_command <-  paste(idx_command,collapse = "")

df_command <- readLines('data/sql/figS6_drosophila_ncRNA.sql')
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
        "reproductive system without gonad and genitalia;
        reproductive system without gonad"
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

df <- inner_join(df_raw,idx_raw,by = "Run") %>% select(-Run) %>%
    select(where(~all(!is.na(.)))) %>%
    #器官・種ごとにまとめて中央値を計算
    group_by(Organism,source_name) %>%
    summarise(across(where(is.numeric), median),.groups = "drop") #グループ化解除

#dbとの接続解除&メモリの開放
dbDisconnect(con)
rm(idx_raw,df_raw)

df_index <- df %>% select(c(Organism,source_name))
df_value <- df %>% select(-c(Organism,source_name)) %>%
    select(where(~all(!is.na(.)))) %>%
    apply(1,normalize_func) %>% t() %>% as_tibble() #正規化

#計算
umap_result <- umap(df_value)

#結果を取り出してlegendフレームと結合
umap_df <- as_tibble(umap_result$layout) %>%
    bind_cols(df_index)

#描写
g <- ggplot(umap_df,aes(x = V1,y = V2,color = source_name)) +
    geom_point(size = 3) +
    stat_ellipse() + #確率楕円
    theme_classic() +
    labs(
        x = "UMAP1",
        y = "UMAP2",
        color = "Organism",
        title = "Drosophila",
        caption = "ncRNA"
    ) +
    theme(
        text = element_text(size = 18),
        plot.title = element_text(hjust = 0.5)
    )
#保存
ggsave("output/figS6/drosophila_ncRNA.pdf",plot = g)