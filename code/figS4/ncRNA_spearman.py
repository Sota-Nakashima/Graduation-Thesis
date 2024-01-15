#ライブラリのインポート
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#データの読み込み
df = pd.read_csv("data/csv/spearman.csv")

#ncRNAのデータのみ抜き出す
df_pict = df[df["RNAkind"] == "ncRNA"]

#信頼区間つき折れ線グラフ
sns.relplot(
    x = "year",y = "value",
    data=df_pict,kind="line",hue="taxon",palette=["#0071BC","#FF5050"])

#細かい調整
plt.xticks([36,0,33,3.95,54,11.6,94,6.4,28.8,87],rotation=45)
plt.xlabel("MYA",fontsize = "14")
plt.ylabel("Spearmans'ρ",fontsize = "14")
plt.title("ncRNA",fontsize="18")
#保存
#bbox_inchesはタイトルの見切れ防止
plt.savefig("output/figS4/ncRNA_spearman.pdf",bbox_inches='tight')
