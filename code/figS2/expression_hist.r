#初期化
rm(list = ls())

#ライブラリの読み込み
library(tidyverse)
library(RMySQL)
library(DBI)
library(ConfigParser)

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
command <- readLines('data/sql/figS2.sql')
command <- command[!grepl("^-",command)]
command <-  paste(command,collapse = "")

#データの取得
df_raw <- as_tibble(dbGetQuery(con,command))

dbDisconnect(con) #dbとの接続解除

df <- df_raw %>%
    #置き換え
    mutate(gene_type = case_when(
        gene_type == "protein_coding" ~ "mRNA",
        gene_type == "lincRNA" ~ "ncRNA",
        TRUE ~ gene_type
        )
        )


g <- ggplot(
    df,
    aes(
        x = testes,
        y = after_stat(density), #正規化
        fill = gene_type
        )
    ) +
    geom_histogram(binwidth = 0.1, position = "identity", alpha = 0.6) +
    labs(
        x = "TPM",
        fill = "") +
    scale_fill_manual(values = c("#0071BC","#FF5050")) +
    #theme要調整
    theme_classic() +
    theme(
        text = element_text(size = 18)
    )

#保存
ggsave(
    "output/figS2/expression_hist.pdf",plot = g,
    width = 7,height = 7)
