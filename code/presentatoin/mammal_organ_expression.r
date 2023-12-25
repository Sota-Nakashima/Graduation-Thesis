#ライブラリの読み込み
library(tidyverse)

#データの読み込み
#あとでSQLから読み込める用に変更
df <- read_delim("data/tmp/mammal.csv",delim = ";")

df <- df %>%
    #ヒトデータを抜き出す
    filter(species == "hg19" & set == "coding") %>%
    select("id","brain","heart","kidney","liver","testes") %>%
    #melt
    pivot_longer(cols = -id,names_to = "organ",values_to = "TPM")


g <- ggplot(df, aes(y = organ, x = TPM,fill = organ)) +
    #バイオリンプロットとボックスプロットを重ねがけ
    geom_violin(show.legend = FALSE) +
    geom_boxplot(width = .1, fill = "white",outlier.color = NA) +
    #theme要調整
    theme_classic() +
    theme(
        axis.title.y = element_blank()
        )

ggsave("output/presentation/mammal_organ_expression.pdf",plot = g)