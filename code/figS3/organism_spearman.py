#ライブラリのインポート
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

#データの読み込み
df = pd.read_csv("data/csv/spearman.csv")

#yearをlog1p変換
df.loc[
    df["taxon"] == "mammal", 'year'
    ] = np.log1p(df.loc[df["taxon"] == "mammal", 'year'])

#mammal x mRNAのデータのみ抜き出す
df_pict = df[(df["RNAkind"] == "mRNA") & (df["taxon"] == "mammal")]

#折れ線グラフ
sns.relplot(
    x = "year",y = "value",
    data=df_pict,kind="line",hue="Organism",
    palette=["#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3"])

#細かい調整
plt.xticks(np.log1p([0,6.4,28.8,87]),[0,6.4,28.8,87],rotation=45) 
plt.xlabel("Divergence Time (MYA)",fontsize = "14")
plt.ylabel("Spearmans'ρ",fontsize = "14")
plt.title("mRNA",fontsize="18")
#保存
#bbox_inchesはタイトルの見切れ防止
plt.savefig("output/figS3/organism_spearman.pdf",bbox_inches='tight')
