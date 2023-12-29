#ライブラリの読み込み
library(tidyverse)
library(RColorBrewer)
library(umap)

#シード値の設定
set.seed(1234)

#データの読み込み
#あとでMySQL仕様に変更
df <- read_csv("data/tmp/mRNA.csv")
df <- df %>% filter(source_name %in% organism_list) %>%
    select(where(~all(!is.na(.)))) %>% select(-sex) %>%
    #器官・種ごとにまとめて中央値を計算
    group_by(Organism,source_name) %>%
    summarise(across(where(is.numeric), median),.groups = "drop") #グループ化解除

#計算
umap_result <- umap(df %>% select(-c(Organism,source_name)))

#結果を取り出してlegendフレームと結合
umap_df <- as_tibble(umap_result$layout) %>%
    bind_cols(df %>% select(c(Organism,source_name)))

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
        caption = "mRNA"
    ) +
    theme(
        text = element_text(size = 18),
        plot.title = element_text(hjust = 0.5)
    )
#保存
ggsave("output/figS6/drosophila_mRNA.pdf",plot = g)