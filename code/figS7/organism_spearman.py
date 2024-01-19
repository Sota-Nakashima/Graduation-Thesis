#ライブラリのインポート
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#データの読み込み
df = pd.read_csv("data/csv/spearman.csv")

#mammal x mRNAのデータのみ抜き出す
df_pict = df[(df["RNAkind"] == "mRNA") & (df["taxon"] == "mammal")]

#折れ線グラフ
sns.relplot(
    x = "year",y = "value",
    data=df_pict,kind="line",hue="Organism",
    palette=["#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3"])

#細かい調整
plt.xticks([94,0,6.4,28.8,87],rotation=45)
plt.xlabel("MYA",fontsize = "14")
plt.ylabel("Spearmans'ρ",fontsize = "14")
plt.title("mRNA",fontsize="18")
#保存
#bbox_inchesはタイトルの見切れ防止
plt.savefig("output/figS7/organism_spearman.pdf",bbox_inches='tight')
