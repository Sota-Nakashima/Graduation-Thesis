library(tidyverse)

organism_list <- c("brain","heart","kidney","liver","testes")
year <- c(0,6.4,28.8,94,87.1,87)
df_raw <- read_delim(
    "/Users/nakashimasota/Downloads/article/data/supplementary-data.csv",
    delim = ";"
    )

mammal_mRNA <- df_raw %>% filter(gene_type == "protein_coding") %>%
    select(id,species,brain,heart,kidney,liver,testes) %>%
    pivot_longer(
        cols = -c("id","species"),
        names_to = "Organism"
        )

mammal_ncRNA <- df_raw %>% filter(gene_type == "lincRNA") %>%
    select(id,species,brain,heart,kidney,liver,testes) %>%
    pivot_longer(
        cols = -c("id","species"),
        names_to = "Organism"
        )

calc_spearman <- function(
    df, #データフレーム
    organism_list, #器官のリスト
    taxon, #mammal or drosophila
    kind, #mRNA or ncRNA
    year #分岐年代のリスト
){
    #空のデータフレームの作成
    result_df <- tibble(
        value = double(0),
        Organism = character(0),
        RNAkind = character(0),
        species = character(0),
        year = double(0)
    )

    for (organism in organism_list) {
        #器官ごとにデータを分割＆前処理
        df_organism <- df %>% filter(Organism == organism) %>%
        select(-Organism) %>% pivot_wider(names_from = id,values_from = value) %>%
        mutate_all(~ifelse(is.na(.), 0, .))

        #データをデータ部とインデックス部に分割
        df_index <- df_organism %>% select(species)
        df_value <- df_organism %>% select(-species)

        #データをcor()に入れるための前処理
        df_calc <- df_value %>% t()
        colnames(df_calc) <- df_index$species #転置でカラム名が消えるので再度代入

        #計算
        tmp_result <- as_tibble(cor(df_calc,method = "spearman")[,1])

        #付属情報を付加
        tmp_result <- tmp_result %>%
        mutate(
            taxon = taxon,
            species = df_index$species,
            RNAkind = kind,
            Organism = organism,
            year = year
            )
        
        #他の結果との結合
        result_df <- bind_rows(result_df,tmp_result)
    }
    return(result_df)
}

print(
    calc_spearman(df = mammal_mRNA,organism_list = organism_list,taxon = "mammal",kind = "mRNA",year = year),
    n = 200)

print(
    calc_spearman(df = mammal_ncRNA,organism_list = organism_list,taxon = "mammal",kind = "ncRNA",year = year),
    n = 200)

#tryCatch(
#    {
#
#    },
#    error = function(e){
#        
#    }
#)
