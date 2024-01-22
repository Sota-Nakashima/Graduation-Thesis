#ライブラリのインポート
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#データの読み込み
df = pd.read_csv("data/csv/spearman.csv")

#mammalのデータのみ抜き出す
df_pict = df[df["taxon"] == "mammal"]

#信頼区間つき折れ線グラフ
sns.relplot(
    x = "year",y = "value",
    data=df_pict,kind="line",hue="RNAkind",palette=["#0071BC","#FF5050"])

#細かい調整
plt.xticks([94,0,6.4,28.8,87],rotation=45) 
plt.xlabel("MYA",fontsize = "14")
plt.ylabel("Spearmans'ρ",fontsize = "14")
plt.title("Mammal",fontsize="18")
#保存
#bbox_inchesはタイトルの見切れ防止
plt.savefig("output/fig5/mammal_spearman.pdf",bbox_inches='tight')
