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
command <- readLines('data/sql/presentation/mammal_organ_expression.sql')
command <- command[!grepl("^-",command)]
command <-  paste(command,collapse = "")

#データの取得
df_raw <- as_tibble(dbGetQuery(con,command))

dbDisconnect(con) #dbとの接続解除

df <- df_raw %>%
    #melt
    pivot_longer(cols = -id,names_to = "organ",values_to = "TPM")


g <- ggplot(df, aes(y = organ, x = TPM,fill = organ)) +
    #バイオリンプロットとボックスプロットを重ねがけ
    geom_violin(show.legend = FALSE) +
    geom_boxplot(width = .1, fill = "white",outlier.color = NA) +
    #theme要調整
    theme_classic() +
    theme(
        axis.title.y = element_blank(),
        text = element_text(size = 18)
        ) +
    scale_fill_brewer(palette = "Pastel1")

ggsave(
    "output/presentation/mammal_organ_expression.pdf",plot = g,
    width = 7,height = 7)