[English version](Example_CDK20.en.md)

# 実例: ヒトCDK20（UniProt Q8IZL9）— AlphaFoldモデル1個から検証済みポケット、ドッキング可能なコンフォメーション・アンサンブルまで

`dd_afpocket`の4段階 —— **fetch → pocket → sample → cluster** —— を実際のターゲット、
ヒトCDK20（Cyclin-dependent kinase 20、遺伝子名`CDK20`/`CCRK`/`CDCH`、
UniProt [Q8IZL9](https://www.uniprot.org/uniprotkb/Q8IZL9)、345残基）に対して
end-to-endで実行した、実データによる完全なワークスルー例です。以下の数値・画像は
すべて実際に実行した結果であり、モックアップではありません。関連ファイル
（インタラクティブHTMLビュー、JSONレポート、CSV、PNGスクリーンショット）は
すべて本リポジトリの[`examples/cdk20/`](examples/cdk20/)にチェックインされています。

## なぜCDK20か

CDK20には広く入手可能な共結晶阻害剤構造がなく、これまで`dd_afpocket`で
検証されたこともないターゲットでした。実験構造や既知リガンドが一切ない
AlphaFoldモデル単独からのポケット検出が、UniProtの独立にキュレーションされた
活性部位アノテーションと照らして本当に生物学的に意味のあるポケットを
見つけられるかどうかを試す、良いテストケースです。

## 1. Fetch: AlphaFold DBモデル → MDグレード構造

```bash
dd_afpocket-fetch Q8IZL9 -o data/prepped
```

```
[q8izl9] chains=A residues=345 chain_breaks=0 disulfides=0
[q8izl9] PDBFixer: 0 missing residue(s) found, 1 missing atom(s) added
[q8izl9] PDBFixer: added atom(s) OXT to A:345 GLU
[q8izl9] pLDDT: very_low=8%, low=6%, confident=21%, very_high=65%; trimmed 0 N-term / 1 C-term residue(s)
[q8izl9] -> data/prepped/q8izl9_md.pdb
```

AlphaFold自身の指標（pLDDT≥90）でモデルの65%が「非常に高信頼度」であり、
トリムが必要だったのはC末の低信頼度残基1個のみ —— ポケット検出の出発点として
十分信頼できる構造です。

## 2. ポケット検出: `dd_afpocket-pocket`は本物のATP部位を見つけるか？

```bash
dd_afpocket-pocket data/prepped/q8izl9_md.pdb -o data/prepped/q8izl9_pocket --visualize
```

fpocketは構造全体から**28個の候補ポケット**を検出しました。Druggability Score順に
並べると、rank 1とそれ以降の差は歴然としています:

![全28候補ポケットのdruggability score分布](examples/cdk20/images/druggability_scores.png)

| Rank | Druggability | Alpha sphere数 | 体積 (Å³) |
|---|---|---|---|
| **1** | **0.646** | 95 | 794.1 |
| 2 | 0.023 | 33 | 369.4 |
| 3 | 0.008 | 33 | 338.7 |
| 4 | 0.005 | 38 | 523.4 |
| 5 | 0.004 | 43 | 338.6 |
| 6 | 0.004 | 38 | 383.2 |
| 7 | 0.004 | 32 | 353.5 |
| 8 | 0.002 | 18 | 91.8 |
| 9 | 0.001 | 21 | 170.6 |
| 10 | 0.001 | 28 | 237.0 |
| 11〜28 | ≤0.001 | — | — |

Rank 1はrank 2よりも**約28倍**druggableです。rank 2以降は実質的にドラッガブルでない、
浅い表面の窪みにすぎません。

### Rank 1をUniProtの公式アノテーションと照合する

`dd_afpocket-pocket`はCDK20の生物学的な知識を一切持たず、幾何学的情報のみを見ています。
そのため本当のテストは「その最上位候補が、この蛋白質について*独立に*知られている
情報と一致するか」です。UniProt Q8IZL9はプロテインキナーゼドメイン（残基4〜288）と
その触媒・結合残基を具体的にアノテーションしています。Rank 1の裏打ち残基24個には、
UniProtが挙げる残基がまさに含まれていました:

| dd_afpocketが検出した残基 | UniProtアノテーション |
|---|---|
| A:10, 11, 12, 13, 14, 15, 18 | **ATP結合部位1** — Gly-richループ、残基10〜18 |
| A:33 | **ATP結合部位2** — 触媒Lys33 |
| A:127 | **活性部位** — 触媒ループのプロトン受容体 |
| A:31, 51, 65, 81, 82, 84, 129, 131, 132, 134, 144, 145, 147, 148, 305 | ヒンジ／触媒ループ近傍（典型的なATPポケット壁面残基） |

これは文献でキュレーションされたATP結合キナーゼポケットとほぼ完全に一致しており、
しかも検出器にはCDK20についての事前知識は一切与えられておらず、純粋に幾何学のみから
得られた結果です。

### 上位3候補を可視化する — 残基だけでなく、実際の空洞ボリュームも

`--visualize`（`dd_afpocket` v0.2.0で追加、詳細は本体の[README](README.jp.md)参照）は
[`pocket_candidates.html`](examples/cdk20/pocket_candidates.html)を書き出します —
サーバー不要でブラウザに直接開けるインタラクティブなビューで、上位3候補ポケットを
比較できます: 各ポケットの裏打ち残基を色分けしたスティックで、fpocketのアルファ
スフィア群を単色の球クラスタとして（裏打ち残基だけでなく*実際に空洞が占める体積*を
表現）、さらに残基ごとに1つのラベルを表示します。

![上位3候補ポケット: 裏打ち残基・空洞ボリューム・残基ラベル](examples/cdk20/images/pocket_candidates.png)

青（rank 1、druggability 0.646）は他よりはるかに大きく、キナーゼフォールドの
2本のβストランドの間に本物の閉じた空洞を形成しています。オレンジ（rank 2、0.023）と
緑（rank 3、0.008）は小さく浅く、部分的にしか囲まれていません — 上記のdruggability
scoreの差を視覚的に裏付けています。

Rank 1単体をUniProtアノテーションで色分けしたクローズアップ（Gly-richループは緑、
触媒Lys33は赤、活性部位Asp127はマゼンタ、その他の裏打ち残基はオレンジ）:

![Rank 1ポケットのクローズアップ（UniProtアノテーションで色分け）](examples/cdk20/images/pocket_view.png)

インタラクティブ版: [`pocket_candidates.html`](examples/cdk20/pocket_candidates.html) ・ [`pocket_view.html`](examples/cdk20/pocket_view.html)

生データ: [`pocket_report.json`](examples/cdk20/pocket_report.json) ・ [`pocket_box.json`](examples/cdk20/pocket_box.json)

## 3. Restrained-MDサンプリング: このポケットはどれほど柔軟か

```bash
dd_afpocket-sample data/prepped/q8izl9_md.pdb data/prepped/q8izl9_pocket \
    -o data/sample/q8izl9 --platform OpenCL
```

デフォルトプリセット: 4個の独立レプリカ、各20ps平衡化＋2ns本サンプリング、
GBn2 implicit solvent、`amber14-all`力場、4fsタイムステップ（HMR）。
ポケット近傍以外の276残基は調和拘束され、69残基（ポケット裏打ち残基＋
ポケット重心から1nm以内の全残基）が完全に自由に動けるようにしました。

**実機での注意点（`--platform`について）**: この計算機にはNVIDIA GTX 1660 Tiが
搭載されていますが、`--platform CUDA`は`CUDA_ERROR_UNSUPPORTED_PTX_VERSION`で
失敗しました —— インストール済みドライバ（560.35.05、CUDA 12.6対応）が、
conda-forgeの`openmm`ビルドがコンパイルされたCUDAツールキット（12.9）より古く、
OpenMMのCUDAプラグインがドライバのサポート範囲を超えるツールキット向けにビルドされた
PTXの読み込みを拒否したためです。`--platform OpenCL`は同じGPUをこのドライバ／
ツールキットのバージョン結合なしに使用でき、問題なく動作しました —— 同様のエラーに
遭遇した場合の参考にしてください。

| レプリカ | フレーム数 | 実時間 |
|---|---|---|
| 1 | 400 | 561秒 |
| 2 | 400 | 570秒 |
| 3 | 400 | 573秒 |
| 4 | 400 | 577秒 |
| **合計** | **1600フレーム（プール）** | **2281秒（約38分）** |

温度は全レプリカを通して290〜308 K（目標300 K）の範囲に収まっており、拘束系が
熱的に安定していること（破綻や大きなドリフトがないこと）を確認できました。

生データ: [`restraint_report.json`](examples/cdk20/restraint_report.json)

## 4. クラスタリングによる代表的コンフォメーション・アンサンブルの生成

```bash
dd_afpocket-cluster data/sample/q8izl9 data/prepped/q8izl9_pocket \
    -o data/clusters/q8izl9 --n-clusters 10 --visualize
```

プールされた1600フレームを（ポケットCA原子のRMSD、average-linkage
階層的クラスタリングで）10個の代表構造（各クラスタのmedoidフレーム、
座標平均は決して使わない）にまとめました:

![クラスタ人口と代表構造のポケット体積](examples/cdk20/images/cluster_populations.png)

| クラスタ | フレーム数 | クラスタ内平均RMSD (Å) | ポケット体積プロキシ (nm³) | 由来レプリカ@時刻 |
|---|---|---|---|---|
| 00 | 1118 | 0.663 | 2.604 | レプリカ3 @ 92 ps |
| 01 | 161 | 0.691 | 2.451 | レプリカ2 @ 357 ps |
| 02 | 129 | 0.635 | 2.596 | レプリカ1 @ 142 ps |
| 03 | 78 | 0.655 | 2.498 | レプリカ2 @ 326 ps |
| 04 | 52 | 0.599 | 2.566 | レプリカ2 @ 311 ps |
| 05 | 37 | 0.684 | 2.613 | レプリカ3 @ 386 ps |
| 06 | 12 | 0.606 | 2.546 | レプリカ4 @ 201 ps |
| 07 | 8 | 0.510 | 2.693 | レプリカ4 @ 33 ps |
| 08 | 3 | 0.427 | 2.558 | レプリカ4 @ 30 ps |
| 09 | 2 | 0.263 | 2.649 | レプリカ2 @ 145 ps |

クラスタ00が圧倒的に支配的な「典型的な」ポケット形状（全プールフレームの70%）で、
残り9クラスタはより稀だが構造的に区別できる側鎖／ループの配置を捉えています —
ポケット体積はアンサンブル全体で2.45〜2.69 nm³の範囲で変動しており、ATP部位周辺で
実際の（控えめではあるが）コンフォメーション多様性が見られます。

10個の代表構造すべてを重ね合わせた図（ポケット裏打ち残基をスティックで、
クラスタごとに色分け、共有の半透明カートゥーンの上に表示）:

![10クラスタの重ね合わせ: ポケット裏打ち残基のコンフォメーション](examples/cdk20/images/cluster_overlay.png)

色ごとに目に見えて異なる側鎖の配置は、restrained MDが単一固定コンフォメーション周りの
熱雑音だけでなく、実際にポケットの柔軟性をサンプリングできていることを裏付けています。

インタラクティブ版: [`cluster_overlay.html`](examples/cdk20/cluster_overlay.html)

生データ: [`cluster_report.csv`](examples/cdk20/cluster_report.csv)

## 結論

UniProtのaccession番号だけから出発し、`dd_afpocket`はCDK20の文献アノテーション済み
ATP結合キナーゼ部位とほぼ完全に一致するポケットを検出し、構造上の他のどの候補よりも
どれだけドラッガブルか（約28倍）を定量化し、実際の（控えめではあるが）ポケット柔軟性を
捉えた10個の代表的受容体コンフォメーションを生成しました —— すべて民生用GPU1枚で
約40分の計算時間で完了しています。

## 次のステップ: アンサンブルドッキング

[`examples/cdk20/`](examples/cdk20/)以下の10個の代表構造は、`pocket_box.json`の
ドッキングボックス（中心・サイズ。apo構造で共結晶リガンドがないためボックスは
ポケット残基から導出）と合わせて、そのまま`dd_docking`のアンサンブルドッキング用
受容体セットとして使用できます —— 具体的なコマンドは本体[README](README.jp.md)の
「Feeding the ensemble into dd_docking」節を参照してください。

## この実例に含まれるファイル

| ファイル | 内容 |
|---|---|
| `examples/cdk20/pocket_report.json` | 選択された（rank 1）ポケット: 残基・中心・druggability score |
| `examples/cdk20/pocket_box.json` | ポケット残基から導出したドッキングボックス（中心・サイズ） |
| `examples/cdk20/pocket_view.html` | インタラクティブ3Dビュー: rank1ポケット、UniProtアノテーションで色分け |
| `examples/cdk20/pocket_candidates.html` | インタラクティブ3Dビュー: 上位3候補ポケット、残基＋空洞ボリューム＋ラベル |
| `examples/cdk20/restraint_report.json` | Restrained-MD実行設定と固定／可動残基数 |
| `examples/cdk20/cluster_report.csv` | 10個の代表構造: フレーム数・RMSD・ポケット体積 |
| `examples/cdk20/cluster_overlay.html` | インタラクティブ3Dビュー: 10個の代表構造すべての重ね合わせ |
| `examples/cdk20/images/*.png` | 本ドキュメントで使用した、上記すべての静的スクリーンショット／グラフ |
