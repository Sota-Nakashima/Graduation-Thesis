#ライブラリの読み込み
library(tidyverse)

#データの読み込み
#あとでSQLから読み込める用に変更
df <- read_delim("data/tmp/mammal.csv",delim = ";")

df <- df %>%
    #ヒトの睾丸のデータを抜き出す
    filter(species == "hg19" & set != "random") %>%
    #RNAの種類と発現量のデータのみを抜き出してくる
    select(set, testes) %>%
    #置き換え
    mutate(set = case_when(
        set == "coding" ~ "mRNA",
        set == "lincs" ~ "ncRNA",
        TRUE ~ set
        )
        )


g <- ggplot(
    df,
    aes(
        x = testes,
        y = after_stat(density), #正規化
        fill = set
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
ggsave("output/figS2/expression_hist.pdf",plot = g)
