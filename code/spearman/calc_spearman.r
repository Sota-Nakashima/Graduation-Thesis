#初期化
rm(list = ls())

#ライブラリの読み込み
library(tidyverse)
library(RMySQL)
library(DBI)
library(ConfigParser)

#コマンド作成
mammal_mRNA_command <- readLines('data/sql/spearman/mammal_mRNA.sql')
mammal_mRNA_command <- mammal_mRNA_command[!grepl("^-",mammal_mRNA_command)]
mammal_mRNA_command <-  paste(mammal_mRNA_command,collapse = "")

mammal_ncRNA_command <- readLines('data/sql/spearman/mammal_ncRNA.sql')
mammal_ncRNA_command <- mammal_ncRNA_command[!grepl("^-",mammal_ncRNA_command)]
mammal_ncRNA_command <-  paste(mammal_ncRNA_command,collapse = "")

drosophila_mRNA_command <- readLines('data/sql/spearman/drosophila_mRNA.sql')
drosophila_mRNA_command <- drosophila_mRNA_command[!grepl("^-",drosophila_mRNA_command)]
drosophila_mRNA_command <-  paste(drosophila_mRNA_command,collapse = "")

drosophila_ncRNA_command <- readLines('data/sql/spearman/drosophila_ncRNA.sql')
drosophila_ncRNA_command <- drosophila_ncRNA_command[!grepl("^-",drosophila_ncRNA_command)]
drosophila_ncRNA_command <-  paste(drosophila_ncRNA_command,collapse = "")

drosophila_idx_command <- readLines('data/sql/spearman/drosophila_idx.sql')
drosophila_idx_command <- drosophila_idx_command[!grepl("^-",drosophila_idx_command)]
drosophila_idx_command <-  paste(drosophila_idx_command,collapse = "")

#db接続用のパスワード読み込み
password_config <- read.ini('password.ini')
password <- password_config$development$password

#付加情報の定義
mammal_organism_list <- c("brain","heart","kidney","liver","testes")
mammal_year <- c(0,6.4,28.8,94,87.1,87)

drosophila_organism_list <- c("abdomen","digestive","gonad","head","thorax")
drosophila_year <- c(33,0,36,3.95,54,11.6)
drosophila_year_sub <- c(33,0,36,54,11.6)

#変更用カラム名
col_list <- c("species","Organism","id","value")

#spearman計算関数
calc_spearman <- function(
    df, #データフレーム
    pivot, #基準となる種
    organism_list, #器官のリスト
    taxon, #mammal or drosophila
    kind, #mRNA or ncRNA
    year, #分岐年代のリスト
    year_sub #drosophila用のオプション
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
        tmp_result <- as_tibble(
            cor(df_calc,method = "spearman")[,which(colnames(df_calc) == pivot)]
            ) #基準となる種の列を抽出

        #付属情報を付加
        if (nrow(tmp_result) == 6){
            tmp_result <- tmp_result %>%
            mutate(
                taxon = taxon,
                species = df_index$species,
                RNAkind = kind,
                Organism = organism,
                year = year #6種ならそのまま代入
                )
        } else {
            tmp_result <- tmp_result %>%
            mutate(
                taxon = taxon,
                species = df_index$species,
                RNAkind = kind,
                Organism = organism,
                year = year_sub #5種なら減らしたものを代入
                )
        }   

        
        #他の結果との結合
        result_df <- bind_rows(result_df,tmp_result)
    }
    return(result_df)
}

tryCatch(
    {
        #db接続
        con <- dbConnect(
            user = 'nakashima',
            MySQL(),
            password = password,
            dbname = 'nakashima_db',
            host = 'localhost',
            port = 3306
        )

        #mammal x mRNA
        #データの取得&前処理
        df_raw <- as_tibble(dbGetQuery(con,mammal_mRNA_command)) %>%
            pivot_longer(
                cols = -c("id","species"),
                names_to = "Organism"
                )
        
        #計算
        result_mammal_mRNA <- calc_spearman(
            df = df_raw,
            pivot = "hg19",
            organism_list = mammal_organism_list,
            taxon = "mammal",
            kind = "mRNA",
            year = mammal_year
        )
        
        #mammal x ncRNA
        #データの取得＆前処理
        df_raw <- as_tibble(dbGetQuery(con,mammal_ncRNA_command)) %>%
            pivot_longer(
                cols = -c("id","species"),
                names_to = "Organism"
                )
        
        #計算
        result_mammal_ncRNA <- calc_spearman(
            df = df_raw,
            pivot = "hg19",
            organism_list = mammal_organism_list,
            taxon = "mammal",
            kind = "ncRNA",
            year = mammal_year
        )

        #drosophilaのindexデータ取得&前処理
        idx_raw <- as_tibble(dbGetQuery(con,drosophila_idx_command)) %>%
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
                )
                #使うデータだけ取る
                ) %>% filter(source_name %in% drosophila_organism_list)

        #drosophila x mRNA
        #発現データ取得＆前処理
        df_raw <- as_tibble(dbGetQuery(con,drosophila_mRNA_command)) %>% 
            distinct() %>% group_by(`Gene ID`,Run) %>%
            summarize(TPM = sum(TPM),.groups = 'drop') %>%
            pivot_wider(names_from = `Gene ID`,values_from = TPM)

        #データの結合
        df_raw <- inner_join(df_raw,idx_raw,by = "Run") %>% select(-Run) %>%
            #器官・種ごとにまとめて中央値を計算
            group_by(Organism,source_name) %>%
            summarise(
                across(where(is.numeric), median),.groups = "drop"
                ) %>% #グループ化解除
            pivot_longer(
                cols = -c(Organism,source_name),
                names_to = "id",
                values_to = "value"
            )
        #計算用にカラム名変更
        colnames(df_raw) <- col_list

        #計算
        result_drosphila_mRNA <- calc_spearman(
            df = df_raw,
            pivot = "D.mel",
            organism_list = drosophila_organism_list,
            taxon = "drosophila",
            kind = "mRNA",
            year = drosophila_year,
            year_sub = drosophila_year_sub
        )

        #drosophila x ncRNA
        #発現データ取得＆前処理
        df_raw <- as_tibble(dbGetQuery(con,drosophila_ncRNA_command)) %>% 
            distinct() %>% group_by(`Gene ID`,Run) %>%
            summarize(TPM = sum(TPM),.groups = 'drop') %>%
            pivot_wider(names_from = `Gene ID`,values_from = TPM)

        #データの結合
        df_raw <- inner_join(df_raw,idx_raw,by = "Run") %>% select(-Run) %>%
            #器官・種ごとにまとめて中央値を計算
            group_by(Organism,source_name) %>%
            summarise(
                across(where(is.numeric), median),.groups = "drop"
                ) %>% #グループ化解除
            pivot_longer(
                cols = -c(Organism,source_name),
                names_to = "id",
                values_to = "value"
            )
        #計算用にカラム名変更
        colnames(df_raw) <- col_list

        #計算
        result_drosphila_ncRNA <- calc_spearman(
            df = df_raw,
            pivot = "D.mel",
            organism_list = drosophila_organism_list,
            taxon = "drosophila",
            kind = "ncRNA",
            year = drosophila_year,
            year_sub = drosophila_year_sub
        )
        #データベース切断
        dbDisconnect(con)

    },
    #エラーが起きた場合エラーを出力してdb接続を切る
    error = function(e){
        print(e)
        dbDisconnect(con)
    }
)

#結果の結合
result <- bind_rows(
    result_mammal_mRNA,
    result_mammal_ncRNA,
    result_drosphila_mRNA,
    result_drosphila_ncRNA
)

#出力
write_csv(result,"data/csv/spearman2.csv")
