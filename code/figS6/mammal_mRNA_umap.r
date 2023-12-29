#ライブラリの読み込み
library(tidyverse)
library(RColorBrewer)

#シード値の設定
set.seed(1234)

#データの読み込み
#あとでMySQL仕様に変更
df <- read_delim("data/tmp/mammal.csv",delim = ";")
df <- df %>% filter(set == "coding") %>%
    select(c(id,species,brain,heart,kidney,liver,testes)) %>%
    pivot_longer(cols = -c(id,species),names_to = "Organism") %>%
    pivot_wider(names_from = id,values_from = value) %>%
    select(where(~all(!is.na(.))))

df_index <- df %>% select(Organism,species)
df_value <- df %>% select(-Organism,-species) %>% select_if(~sum(.) != 0)

#計算
umap_result <- umap(df_value)

#結果を取り出してlegendフレームと結合
umap_df <- as_tibble(umap_result$layout) %>% bind_cols(df_index)

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
ggsave("output/figS6/mammal_mRNA.pdf",plot = g)