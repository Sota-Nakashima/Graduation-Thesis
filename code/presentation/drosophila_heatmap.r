#初期化
rm(list = ls())

#ライブラリの読み込み
library(tidyverse)
library(RMySQL)
library(DBI)
library(ConfigParser)

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
idx_command <- readLines('data/sql/presentation/drosophila_idx.sql')
idx_command <- idx_command[!grepl("^-",idx_command)]
idx_command <-  paste(idx_command,collapse = "")

df_command <- readLines('data/sql/presentation/drosophila_df.sql')
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
df <- inner_join(df_raw,idx_raw,by = "Run") %>%
    select(where(~all(!is.na(.)))) #欠損値を除く

#dbとの接続解除&メモリの開放
dbDisconnect(con)
rm(idx_raw,df_raw)

#indexとdata部分を分離
df_index <- df %>% select(Run) #正規化のためにindex列を退避
df_value <- df %>% select(-c(Organism,source_name,Run)) %>%
    apply(1,normalize_func) %>% t() %>% as_tibble() %>% #正規化
    mutate_all(~ifelse(is.na(.), 0, .)) #NAを0にする

#発現上位10の列を抜き出す
top_10_columns <- df_value %>%
    summarise(across(everything(), sum)) %>%
    pivot_longer(everything()) %>%
    arrange(desc(value)) %>%
    head(10) %>%
    pull(name)

#heatmap用にデータを整形
df_final <- bind_cols(df_index,df_value) %>% #データの結合
    select(Run,top_10_columns) %>% slice(1:10) #10x10のデータを切り出す

#カラム名の変更
colnames(df_final) <- c(
    "Run","RpL3","RpS25","RpS3","RpL10","RpS20",
    "RpL13A","RpL27","RpL6","RpS8","MtnA"
    )

#melt
df_final <- df_final %>% pivot_longer(cols = -Run)


g <- ggplot(df_final,aes(x = Run,y = name,fill = value)) +
    geom_tile(color = "white",linewidth = 1) +
    coord_equal() + #マス目を正方形に
    scale_fill_gradient2(
        high = "#0071BC",
        low = "#FF5050",
        mid = "white",
        midpoint = 0.5,
        limit = c(0, 1),
        name = "Value"
    ) + #色合いの設定
    geom_text(
        aes(label = sprintf("%0.3f", value)),
        color = "black",
        size = 2.5
    )+
    theme_minimal() +
    theme(
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(
            angle = 30,size = 6,hjust = .5
            ),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 6),
        panel.grid=element_blank()
        ) + #テーマの手動設定
    labs(
        x = "Sample",
        y = "Gene"
    )

#保存
ggsave("output/presentation/drosophila_heatmap.pdf",plot = g)