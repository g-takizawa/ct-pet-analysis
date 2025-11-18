# 1all データ解析

`data/raw/1all_data.csv` と `data/raw/2before_after.csv` を対象に、CT と PET の散布図＋回帰線、Pearson 相関を出力するスクリプトと、符号検定／Wilcoxon 符号付順位検定（片側：PET > CT）を行う別スクリプトを用意しています。

## 実行方法

```sh
# 1all データ: 散布図＋回帰・相関
cd path/to/project
Rscript R/scripts/analyze_1all.R

# 1all データ: 符号検定と Wilcoxon 検定（図表つき）
cd path/to/project
Rscript R/scripts/tests_1all.R

# 2before_after データ: 散布図＋回帰・相関
cd path/to/project
Rscript R/scripts/analyze_2before_after.R

# 2before_after データ: 符号検定と Wilcoxon 検定（図表つき）
cd path/to/project
Rscript R/scripts/tests_2before_after.R

# 2データセットをまとめた Markdown レポート
cd path/to/project
Rscript R/scripts/report_summary.R
```

## 出力
### 1all データ
- 図: `output/figures/1all_scatter.{png,pdf,svg}`（geom_count）
- 図: `output/figures/1all_jitter.{png,pdf,svg}`（ドットゆらぎ）
- 相関・回帰ログ: `output/logs/1all_correlation.txt`
- 符号検定／Wilcoxon ログ: `output/logs/1all_nonparametric.txt`
- 符号検定／Wilcoxon 集計表: `output/tables/1all_nonparametric.csv`
- PET-CT 差のヒストグラム: `output/figures/1all_difference_hist.{png,pdf,svg}`
- PET-CT 差の症例別棒グラフ: `output/figures/1all_difference_bar.{png,pdf,svg}`

### 2before_after データ
- 図: `output/figures/2before_after_scatter.{png,pdf,svg}`（geom_count）
- 図: `output/figures/2before_after_jitter.{png,pdf,svg}`（ドットゆらぎ）
- 相関・回帰ログ: `output/logs/2before_after_correlation.txt`
- 符号検定／Wilcoxon ログ: `output/logs/2before_after_nonparametric.txt`
- 符号検定／Wilcoxon 集計表: `output/tables/2before_after_nonparametric.csv`
- PET-CT 差のヒストグラム: `output/figures/2before_after_difference_hist.{png,pdf,svg}`
- PET-CT 差の症例別棒グラフ: `output/figures/2before_after_difference_bar.{png,pdf,svg}`

### レポート
- Markdown サマリー: `output/reports/summary.md`

## 依存パッケージ
- readr
- dplyr
- ggplot2
- tidyr
- svglite（SVG 出力用）

必要に応じて `install.packages()` でインストールしてください。
これ追加
これも追加
