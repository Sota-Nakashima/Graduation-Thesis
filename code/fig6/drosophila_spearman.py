#ライブラリのインポート
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

#データの読み込み
df = pd.read_csv("data/csv/spearman.csv")

#yearをlog1p変換
df.loc[
    df["taxon"] == "drosophila", 'year'
    ] = np.log1p(df.loc[df["taxon"] == "drosophila", 'year'])

#drosophilaのデータのみ抜き出す
df_pict = df[df["taxon"] == "drosophila"]

#信頼区間つき折れ線グラフ
sns.relplot(
    x = "year",y = "value",
    data=df_pict,kind="line",hue="RNAkind",palette=["#0071BC","#FF5050"]
    )

#細かい調整
plt.xticks(np.log1p([36,0,3.95,54,11.6]),[36,0,3.95,54,11.6],rotation=45)
plt.xlabel("Divergence Time (MYA)",fontsize = "14")
plt.ylabel("Spearmans'ρ",fontsize = "14")
plt.title("Drosophila",fontsize="18")
#保存
#bbox_inchesはタイトルの見切れ防止
plt.savefig("output/fig6/drosophila_spearman.pdf",bbox_inches='tight')
