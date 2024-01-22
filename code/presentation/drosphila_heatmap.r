library(tidyverse)

set.seed(123)  # 再現性のためにシードを設定
my_data <- as_tibble(matrix(runif(100), nrow = 10, ncol = 10)) %>%
    mutate(index = 1:10) %>%
    mutate(index = as.character(index)) %>%
    pivot_longer(cols = -index)

# データフレームを表示
print(my_data)


g <- ggplot(my_data,aes(x = index,y = name,fill = value)) +
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
        aes(label = sprintf("%0.2f", value)),
        color = "black",
        size = 2.5
    )+
    theme_minimal() +
    theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 18),
        panel.grid=element_blank()
        ) + #テーマの手動設定
    labs(
        x = "Sample",y = "Gene"
    ) #軸ラベル名の変更

#保存
ggsave("output/presentation/drosophila_heatmap.pdf",plot = g)