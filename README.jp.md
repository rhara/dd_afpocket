# dd_afpocket — AlphaFoldは静止画でも、ポケットは静止しない

静的なAlphaFold DB予測構造を、アンサンブルドッキングに適した少数の受容体コンフォメーションからなるアンサンブルへと変換する: druggableなポケットを検出し、restrained MD（拘束付きMD、ポケット近傍のみが動き、蛋白質の残りは位置拘束される）でそのポケットの局所的な柔軟性をサンプリングし、得られたトラジェクトリを構造的にクラスタリングして少数の代表的なコンフォメーションにまとめる。特定のターゲットに依存しない再利用可能なパッケージとして設計されている（`dd_prep` / `dd_docking` / `dd_overlay` / `dd_viewer` / `dd_confgen` と同じ思想）——概念的には `dd_prep`（AFDB取得＋MDグレードの修復）→ **`dd_afpocket`**（ポケット検出＋局所restrained-MDサンプリング＋クラスタリング）→ `dd_docking`（生成されたコンフォメーション群に対するアンサンブルドッキング）というパイプラインの中間段階にあたる。`dd_afpocket` は他の`dd_*`プロジェクトへの実行時依存を一切持たない: `prep.py` がfetch/repair段階のAFDB専用ロジックを`dd_prep`からvendor（移植）して自己完結している（AlphaFoldモデルには `dd_md` の自己完結的な受容体前処理が対処すべき実PDB由来の癖が一切なく、`dd_prep`が他の消費者向けに提供する汎用的なPDB/ヘテロ原子分類/ドッキング用修復機構は不要なため、AFDB/MD修復のみの狭い経路だけを取り込んだ）。`dd_docking`/`dd_md` もimportせず、それぞれのドッキングボックス規約とharmonic-restraint（調和拘束）機構を踏襲するのみである（後述の「設計メモ」参照）。

- **Fetch (`dd_afpocket-fetch`)**: UniProtアクセッション番号 -> AlphaFold DBモデル ->
  MDグレードの修復。`prep.py` の `fetch_and_prepare_afdb` を使用する。
- **Pocket (`dd_afpocket-pocket`)**: `fpocket` を実行し、Druggability Scoreの順位で
  ポケットを1つ選択する（デフォルトでは最上位）。選択したポケットの
  裏打ち残基／中心座標を書き出し（`pocket_report.json`）、それらの残基の
  座標から導出したドッキングボックスも出力する（`pocket_box.json`。共結晶
  リガンドが存在しない——apo構造であるため——これらの残基座標から
  ボックスを導出している）。
- **Sample (`dd_afpocket-sample`)**: restrained MD。暗黙的溶媒（implicit solvent、
  デフォルトはGBn2、`--implicit-solvent`）または明示的な周期的水ボックス
  （`--solvent explicit`、`--water-model`）を、選択可能な蛋白質フォースフィールド
  （`--protein-forcefield`）の下で使用する。ポケット近傍の外側にある残基は
  開始位置へ調和拘束（harmonically restrained）され、ポケット裏打ち残基
  （加えてポケット重心から `--mobile-radius-nm` 以内にあるもの）のみが
  自由に動く。複数の独立レプリカ（`--n-replicas`）を並列実行する
  （`--n-jobs`）。
- **Cluster (`dd_afpocket-cluster`)**: 全レプリカのフレームをプールし、ポケット原子
  のRMSD距離行列を計算して、クラスタリング（`AgglomerativeClustering`、
  average linkage）により `--n-clusters` 個（デフォルト10）の代表構造に
  まとめる。各代表構造はそのクラスタのmedoidフレームである（座標平均では
  ない）。`--pocket-expand-only` はポケットがMD前の構造に対して幾何学的に
  開いたフレームのみにクラスタリングを限定する。`--visualize` は全代表構造を
  スタンドアロンのHTMLオーバーレイとして書き出す。
- **Run (`dd_afpocket-run`)**: 上記4段階をエンドツーエンドで実行する。

## インストール

rdkit, numpy, pandas, pdbfixer, openmm, mdtraj, matplotlib, scipy,
scikit-learn, py3Dmol（すべて`pyproject.toml`の`dependencies`に明記済み——
他の`dd_*`パッケージへの依存は一切ない）、さらに `fpocket` CLIバイナリ
（`dd_afpocket-pocket`用、pipでは入らないためconda-forge必須）が必要である。
他の`dd_*`プロジェクトとは独立した、専用のconda env `dd_afpocket`
（python 3.12）を使用する:

```bash
mamba create -n dd_afpocket -c conda-forge python=3.12 rdkit numpy pandas pdbfixer openmm mdtraj matplotlib scipy scikit-learn py3dmol pytest fpocket
/opt/miniforge3/envs/dd_afpocket/bin/pip install --no-deps -e .
```

これにより5つのコンソールコマンドがインストールされる: `dd_afpocket-fetch`、
`dd_afpocket-pocket`、`dd_afpocket-sample`、`dd_afpocket-cluster`、`dd_afpocket-run`。

## 使い方

### 1. 取得 (`dd_afpocket-fetch`)

```bash
dd_afpocket-fetch O60674 -o data/prepped
```

ヒトJAK2（UniProt O60674）の現行AlphaFold DBモデルをダウンロードし、
`prep.py` 経由でMDグレードの修復（pLDDTに基づく末端トリミング、
PDBFixer、pH 7.0でのプロトン化）を行い、`data/prepped/o60674_md.pdb` を
書き出す（`dd_afpocket-fetch` は出力を `<uniprot>/` サブディレクトリの下に
ネストしない——それを行うのは以下の複数段階の出力レイアウトを持つ
`dd_afpocket-run` のみである）。

### 2. ポケット検出 (`dd_afpocket-pocket`)

```bash
dd_afpocket-pocket data/prepped/o60674_md.pdb -o data/prepped/o60674_pocket
```

fpocketを実行し、検出された全ポケットをDruggability Scoreの順に出力する。
上記のJAK2モデルで実測された出力（構造全体のマルチドメイン構造にわたって
87個の候補ポケットが検出された。以下は抜粋）:

```
[pocket 1] score=0.158  druggability=0.829  n_alpha_spheres=70
[pocket 2] score=0.056  druggability=0.497  n_alpha_spheres=54
[pocket 3] score=0.027  druggability=0.455  n_alpha_spheres=111
...
[done] pocket rank 1 (druggability=0.829, 17 residue(s)) -> data/prepped/o60674_pocket/pocket_report.json
```

`--pocket-rank N` で別の順位を選択できる（例えば既知の非最上位アロステリック
サイトなど）。`--pocket-residues A:42,A:87,...` はfpocketの残基検出を完全に
バイパスして手動指定のサイトを使用する。

### 3. Restrained-MDサンプリング (`dd_afpocket-sample`)

```bash
dd_afpocket-sample data/prepped/o60674_md.pdb data/prepped/o60674_pocket \
    -o data/sample/o60674 --n-replicas 4 --n-jobs 2
```

ポケット近傍の外側にあるすべての残基を位置拘束し（harmonic、デフォルトで
`k=1000 kJ/mol/nm^2` —— `dd_md/restraints.py` がこの値で300 Kにおいて
熱揺らぎのRMSが~0.86 Åになることを検証済みであり、ポケット形状の
違いというクラスタリングが探す対象のスケールに比べて小さい）、
`--n-replicas` 個の独立した暗黙的溶媒MDレプリカを実行する。上記のJAK2
ポケットの例では、17個の裏打ち残基とその重心から1 nm以内にあるすべての
残基により、全1097残基のうち49残基が可動になった
（`restraint_report.json` の `residues_mobile`）——これほど大きなターゲット
での現実的な所要時間については後述の「性能」を参照。

`--preset {default,quick}` は `--n-replicas`/`--equil-ps`/`--sample-ns`/
`--report-ps`/`--timestep-fs`/`--implicit-solvent` を1つのフラグにまとめる。
`--preset quick` はCPUのみでも扱いやすい粗い設定（2レプリカ、5 psの平衡化、
300 psの本サンプリング、2.5 psのレポート間隔、`obc2` 暗黙的溶媒）であり、
下流のアンサンブルドッキングにポケット形状の大まかな多様性さえ得られれば
よく、収束したトラジェクトリまでは必要ない場合に用いる——`dd_afpocket` はそもそも
単一の安定性判定ではなくコンフォメーションのアンサンブルを生成するツール
であることを踏まえれば妥当な設計である（「既知の限界」参照）。上記6つの
フラグのいずれかを明示的に指定した場合、そのフラグに関してはpresetの値を
上書きする。例: `--preset quick --n-replicas 4`。`--timestep-fs`（デフォルト
4 fs、両presetともそのまま）は別扱いになっている。`restraints.py` が
デフォルトの `k=1000 kJ/mol/nm^2` 拘束に対して安定性を実測検証した唯一の
値であり、これ以上大きくすることはpreset側の判断に委ねられるトレード
オフではないためである。

`--implicit-solvent {gbn2,gbn,obc2,obc1,hct}` はGBモデルを選択する
（デフォルト `gbn2`）。本プロジェクトの開発機で実測（ヒトリゾチーム、
P61626、~2300原子、CPU/4スレッド、`gbn2` を基準としたms/stepの相対値）:
`hct` は1.75倍高速、`obc1` は1.55倍、`obc2` は1.37倍、`gbn` は~1.02倍
（`gbn2` は `gbn` に「neck」Born半径補正項を加えたものであり、これが
`gbn` に対する追加コストのほぼ全てである）。`--preset quick` は `obc2` を
選ぶ: 最速の `hct` は最も古いGBモデルであり、埋没原子のBorn半径を
過小評価することが知られている。ここではポケット裏打ち側鎖の一部が
部分的に埋没しているため、これは無視できない。`obc2` はわずかな速度上の
コストでこれを解決しており、AMBER系ツール群における事実上の標準GB
モデルである（`igb=5`）。

### フォースフィールドと溶媒モデル

`--protein-forcefield {amber14-all,amber99sbildn,amber19-all,charmm36}`
（デフォルト `amber14-all`）は蛋白質のXMLを選択する。すべての `--implicit-
solvent` モデルとの全組み合わせがビルド可能であることを検証済みである。
`--solvent {implicit,explicit}`（デフォルト `implicit`）は、GB連続溶媒
（本プロジェクトの他のCPUフレンドリーなデフォルトが前提とする設定）と、
実際の周期的水ボックスとを切り替える。`--solvent explicit` には
`--water-model` の指定が必要であり、`--protein-forcefield` ごとに
検証済みの組み合わせがある（`sample.PROTEIN_FORCEFIELDS`）:

| `--protein-forcefield` | 検証済みの `--water-model` の選択肢 |
|---|---|
| `amber14-all` | `tip3p`, `tip3pfb`, `tip4pew`, `tip4pfb`, `spce` |
| `amber99sbildn` | `tip3p`, `tip4pew`, `tip5p` |
| `amber19-all` | `tip3p`, `tip3pfb`, `tip4pew` |
| `charmm36` | `tip3p`, `tip4pew`, `tip5p` |

明示的溶媒ではさらに `--solvent-padding-nm`（溶質からボックス端までの
パディング、デフォルト1.0）、`--ion-concentration-molar`（Na+/Cl-の
中和／背景濃度、デフォルト0.15、生理的濃度）、`--pressure-atm`
（`MonteCarloBarostat` の目標値、デフォルト1.0——明示的溶媒のランはNVTでは
なくNPTである）が追加される。溶媒和（`Modeller.addSolvent`、シード
パラメータを持たない）は `dd_afpocket-sample` の1回の呼び出しにつき、いずれの
レプリカが始まる前にちょうど1回だけ行われ、すべてのレプリカがその
共有済みの、既に周期的な構造を読み込む——レプリカごとに独立して溶媒和
すると、各レプリカのDCDが `dd_afpocket-cluster` が全レプリカに対して読み込む
`complex_top.pdb` と異なる原子数／順序になってしまう。

**明示的溶媒はCPU上では暗黙的溶媒より劇的に遅い**——本プロジェクトの
開発機で実測: ヒトリゾチームを~0.6 nmパディングで溶媒和した場合
（amber14-all + tip3p、水・イオンを含めて合計105,716原子、CPU、
スレッド数固定なし）、1レプリカの正しさ確認ラン（最小化 -> 0.08 ps平衡化
-> 20 ps本サンプリング、`--n-jobs 1`）はwall時間で2,304秒かかり、本
サンプリング自体は最初の数psを過ぎると~87 s/psで安定した——これに対し、
同じ蛋白質自身の暗黙的GBn2での数値（後述の「性能」参照、~2,300
原子、CPU/4スレッド固定、4 fsタイムステップで~0.09 s/step、すなわち
~22.7 s/ps）と比べると、これは厳密なapples-to-apples比較ではないものの
（明示的溶媒系の105,716原子は暗黙的溶媒系の原子数の約46倍、大半が水で
あり、両ランでスレッド固定の扱いも異なる）、psあたりで見てすでに数倍
遅い。実務上の結論はどちらにせよ明確である: 簡単な正しさ確認を超える
`--solvent explicit` のランには `--platform CUDA` を強く推奨する。
CPUのみの環境では `--solvent implicit`（デフォルト）が実用的な選択肢
であり続ける（後述の「性能」参照）。`dd_afpocket-cluster` が
`>99,999` 原子の明示的溶媒 `complex_top.pdb` を読み込む際には
`mdtraj... Need to guess atom number ...` という警告も出る——これは
無害であり、固定幅PDBの原子シリアル番号フィールドが99,999を超えて
オーバーフローしているだけである（mdtrajはこのテキストフィールドでは
なく、インデックスで原子を追跡している）。

### 4. クラスタリング (`dd_afpocket-cluster`)

```bash
dd_afpocket-cluster data/sample/o60674 data/prepped/o60674_pocket -o data/clusters/o60674 --n-clusters 10
```

`cluster_00.pdb`（最大のクラスタ）から `cluster_09.pdb` までを書き出し、
さらに `cluster_report.csv`（母集団サイズ、由来レプリカ／時刻、代表構造
ごとのクラスタ内平均RMSD、および各代表構造の `pocket_volume_nm3`——
以下を参照）も出力する。

#### 事後的なポケット拡張フィルタリング (`--pocket-expand-only`)

`cluster.pocket_volume_proxy` は「ポケットがどれだけ開いているか」の
安価な、再実行不要な幾何学的プロキシを計算する: 各プールされたフレーム
におけるポケット裏打ち原子（`--cluster-atoms`）の凸包体積である——
fpocket自身の空洞体積ではない。fpocketをフレームごとに再実行することに
なり、あまりに遅すぎるためである（「性能」参照）。
`--pocket-expand-only` は、体積が `>= reference * (1 + --pocket-expand-
margin)` であるフレームのみにクラスタリングを限定する。`reference` は
同じプロキシをMD前の構造に対して計算した値である——これはすでに
サンプリング済みの、それ自体は偏りのないアンサンブルに対する後付けの
選別であり（`dd_afpocket-sample` 自体は「拡張」について何も知らない。フレームの
生成のされ方は何も変わらず、どのフレームが代表構造の候補になり得るかが
変わるだけである）。

**`--pocket-expand-margin` は較正が必要である——デフォルト（0.0）では
ほとんど何もフィルタされない。** 凸包はその点群を包む外皮であるため、
実際の「開き」のシグナルとは無関係に、熱的なノイズだけで単一の静的な
参照フレームに対して膨張する傾向がある——実際の（短い）サンプリング
ランで実測したところ、マージン0.0では10,000フレーム中10,000フレームが
そのまま残った。`--pocket-expand-only` は常に `[cluster] pocket-expand-
only: kept X/Y frame(s) ...` を出力する。`--pocket-expand-margin` は
0.05〜0.1あたりから始め、デフォルトが何かをしてくれると期待するのでは
なく、この出力される比率を見ながら調整すること。

#### 視覚的な比較 (`--visualize`)

`--out-dir` に `cluster_overlay.html` を書き出す: サーバ不要の
スタンドアロンなpy3Dmolシーンで、各 `cluster_NN.pdb` がそれぞれ独立した
モデルとして読み込まれ、それぞれに異なる色が割り当てられる（受容体
全体は半透明のcartoon、ポケット裏打ち残基の部分は不透明のstickで表示
される——実際にクラスタ間で変化するのはそこだけであり、それ以外の
部分はサンプリング中 `restraints.py` によって固定されている）。事前の
構造アラインメントは不要である: すべてのクラスタPDBはすでにdd_afpocketの
restrained-MDの基準フレームを共有している。このファイルを開くには
インターネット接続が必要である（3Dmol.jsはCDNから読み込まれる。py3Dmol
の通常のJupyter使用時と同様）。

### エンドツーエンド (`dd_afpocket-run`)

```bash
dd_afpocket-run O60674 -o data --n-replicas 4 --n-jobs 2 --n-clusters 10
```

上記4段階すべてを実行し、`data/o60674/` 以下にすべてを書き出す。

### アンサンブルを `dd_docking` に渡す

N個の代表構造は蛋白質のみ（apo）のPDBであり、そのまま `dd_docking` の
アンサンブルメンバー入力として使用できる。ドッキングボックスを導出する
共結晶リガンドが存在しないため、`pocket_box.json` のボックス（中心／
サイズ）を全メンバーに対して使用すること——それらは同じポケット定義と
残基フレームを共有しているため、ボックスは全メンバーにわたって一貫して
いる。

## 設計メモ

- **プレーンな `openmm.app.ForceField` であり、`SystemGenerator` ではない。**
  `dd_docking/refine_md.py` と `dd_md/system_build.py` が `SystemGenerator`
  を必要とするのは、そのシステムに低分子リガンドが含まれるためである
  （GAFF/SMIRNOFFパラメータ化）。`dd_afpocket` のシステムはapoであり——リガンドを
  含むことは決してない——プレーンな `ForceField(...)` だけで十分である。
  これは実際の環境上の問題も回避する: 明示的な
  `small_molecule_forcefield="gaff-*"` を指定せずに `SystemGenerator` を
  構築すると、`openff.toolkit` のSMIRNOFFフォースフィールド探索が即座に
  トリガーされ、`mpro` env上では `ModuleNotFoundError: No module named
  'pkg_resources'`（`openff-amber-ff-ports` がsetuptools >= 81ではもはや
  提供されない `pkg_resources` APIに依存しているため）で失敗する。これが
  `dd_docking`/`dd_md` に影響しないことは確認済みである。両者とも常に
  `small_molecule_forcefield="gaff-2.11"` を明示的に渡しており、この壊れた
  importに到達しない別の（GAFFの）コードパスを通るためである。
- **精選され、実証済みのフォースフィールド／水モデルのレジストリ
  （`sample.PROTEIN_FORCEFIELDS`）であり、「OpenMMが提供する全XML」では
  ない。** OpenMMはdd_afpocketが公開しているよりもはるかに多くの蛋白質／水／GB
  パラメータファイルを同梱している。未検証の組み合わせ（例えばAMBERで
  パラメータ化されたGBファイルをCHARMMの原子タイプに対して使うなど）は
  完全に失敗するか、あるいは黙って不整合なシステムを構築してしまう
  可能性がある。掲載されているすべての `--protein-forcefield` ×
  `--implicit-solvent` の組み合わせ、および `--protein-forcefield` ×
  `--water-model` の組み合わせは、追加される前に本プロジェクトの開発機上で
  実際にビルドに成功することが（`ForceField(...).createSystem(...)` /
  `Modeller.addSolvent(...)` + `createSystem(...)` により）確認されている。
- **明示的溶媒は `dd_afpocket-sample` の呼び出しにつきちょうど1回だけ溶媒和され、
  全レプリカで共有される。** `Modeller.addSolvent` にはシードパラメータが
  ない——そのイオン配置はランダム化される——そのためレプリカごとに独立に
  溶媒和すると（これは*暗黙的*溶媒のシステム構築、つまり決定論的な処理
  では問題にならない）、各レプリカのDCDが、`dd_afpocket-cluster` が最終的に
  全レプリカのDCDに対して読み込む `complex_top.pdb` と異なる原子数／
  順序を持つことになってしまう。`sample.sample_pocket` は1回だけ溶媒和を
  行い、周期的構造を `solvated_input.pdb` に書き出し、すべてのレプリカが
  そのファイルをプレーンな、決定論的な入力として読み込む——apo
  （暗黙的溶媒）の場合と構造的に同じ扱いである。
- **明示的溶媒の場合のみ `rigidWater=True` と `MonteCarloBarostat` を使う。**
  暗黙的溶媒のシステムには剛体化する水が存在しない（このフラグは効果を
  持たない）。明示的溶媒の水は慣例的に剛体（結合・角度拘束であり、
  水素質量再分配はしない）として扱われ、蛋白質自身の4 fsのHMR
  タイムステップで積分されることはない。barostatは平衡化・本サンプリング
  を通じて `--pressure-atm`（デフォルト1気圧）でNPTを実行し続ける——
  `restraints.py` の調和拘束（文字通り凍結しているわけではない）は
  barostatの周期的な体積リスケーリングの試行を妨げない。単に拘束された
  原子を、他のあらゆる摂動と同様に、その後で元の位置へ引き戻すだけである。
- **`NoCutoff` ではなく `CutoffNonPeriodic`。** GBn2のBorn半径／GBエネルギー
  項はすでに遠方の原子による遠距離静電効果を近似しているが、pairwiseの
  非結合ループ自体はその近似とは無関係に素朴なO(n_atoms^2)評価である。
  大きなマルチドメイン蛋白質（全長JAK2: 1132残基、水素を含めて~17800
  原子）では、`NoCutoff` はCPU上でわずか数psであっても実用に耐えないほど
  遅いことが実測された。`restraints.py` はどのみちポケット近傍~1 nm圏外の
  すべてを凍結するため、1.5 nmの非結合カットオフ（デフォルト、
  `--nonbonded-cutoff-nm`）は、本プロジェクトが実際にサンプリングする
  領域についてはほぼ精度上のコストを伴わない。
- **並列レプリカ間でのCPUスレッドの固定。** CPUプラットフォームの
  デフォルトのスレッド数は通常「すべての論理コア」である。`--n-jobs` で
  複数レプリカを各自のスレッド数を制限せずに並行実行すると、N個の
  プロセスが同じコアを奪い合うことになり、実際には並列に走らない——
  `sample_pocket` は `n_jobs != 1` のとき各レプリカのスレッド数を
  `os.cpu_count() // n_workers` に固定する。これは `dd_docking/
  screening.py` が並列Vinaワーカーに対して使っているのと同じパターンで
  ある（`n_jobs == 1`、デフォルトでは、単一レプリカがマシン全体を使える
  ようスレッド数を固定しない）。

## 性能

暗黙的溶媒（GBn2）のCPU MDは速くはなく、実際に積分されるのはポケット近傍
のみであるにもかかわらず、力の評価コストは総原子数に応じてスケールする
（拘束されている領域であっても毎ステップ力の計算が必要である）。実測
（`dd_afpocket-sample`、`--n-jobs 2` により2レプリカ、レプリカあたり~8個のCPU
スレッドを固定）、本プロジェクトの開発機（共有の16コアデスクトップ、
専用の計算ノードではない）:

- Streptavidin（1STP、121残基、水素を含めて1744原子）: レプリカあたり
  合計5250ステップ（平衡化＋本サンプリング）が384秒wallで完了——
  ~0.073 s/step。
- ヒトリゾチーム（P61626、146残基、~2300原子）: レプリカあたり2100ステップ
  が232秒wallで完了——~0.110 s/step。
- 全長JAK2（O60674、1132残基、~17800原子）は、`restraints.py` により
  実際に動く残基は~50個に限られているにもかかわらず、単に力の評価
  コストが大きいために同じスレッド数で1ステップあたり大幅に遅い——
  短い正しさ確認ラン（525ステップ）は、この共有マシン上で他のジョブと
  コアを奪い合っている状況では、スレッドあたり数CPU分以内には完了
  しなかった。小さな蛋白質の上記の数値が原子数に対して線形にスケール
  すると仮定するのではなく、それに応じて予算を確保すること（あるいは
  GPUを使うこと）。ポケット検出（`dd_afpocket-pocket`、fpocket）は蛋白質サイズの
  影響を実質的に受けない——JAK2上でも（候補ポケット87個、最上位の
  druggability 0.829）数秒で完了した。

**実用的な所要時間のための推奨事項:**
- 簡単な正しさ確認を超える用途にはCUDA GPU（`--platform CUDA`）を優先
  すること——暗黙的溶媒のGBカーネルはCPUよりGPU上で劇的に高速である。
- 全長JAK2のような大きなマルチドメインターゲットについては、ワークフロー
  上可能であればより小さなコンストラクト（例えば孤立したキナーゼドメイン）
  を提供することを検討すること——`dd_afpocket` は現時点では残基範囲の
  スライシングを提供しておらず、鎖単位での選択のみを `prep.py` 経由で
  提供している。
- 探索的なランでは `--sample-ns`/`--equil-ps`/`--n-replicas` を小さくし、
  自分のハードウェア上でパイプラインが想定通りに動作することを確認して
  から元に戻すこと。ここには専用の「まずスクリーニングしてから確定」
  というゲートは存在しない（このプロジェクトはコンフォメーションの
  アンサンブルを作るものであり、単一の安定性判定を行うものではないため、
  `dd_md` のスクリーニング→確定フローのように悪いランを早期に打ち切る
  自然な場所がない）。`--preset quick`（上記「Restrained-MDサンプリング」
  参照）はこの縮小版バンドルである——より安価なGBモデル
  （`--implicit-solvent obc2`）も含め、CPUのみのマシン向けにあらかじめ
  選ばれている。`--n-jobs -1` と組み合わせればすべてのコアも使用できる。

## 既知の限界

- apoであるため、`dd_afpocket` はリガンド誘起のコンフォメーション変化
  （induced fit）を再現できない——計算的に検出されたポケット周りの
  非結合蛋白質自身の局所的な柔軟性のみである。ターゲットにとってこれが
  重要な場合は、下流で `dd_docking` 自身のinduced-fitリファインメント
  （`dd_docking-refine`）と組み合わせること。
- ポケット検出の質は完全にfpocketに依存する。fpocketが見逃したり
  誤ってスコア付けしたりしたポケットは、`--pocket-residues` による
  手動指定を除いて回復できない。
- 残基範囲のスライシングがない: `dd_afpocket` は `prep.py` が取得した鎖を
  そのままサンプリングする。マルチドメイン蛋白質の場合、これは
  対象のポケットに実際に関係するドメインよりもはるかに大きく
  (そして遅く)なりうる。

## モジュール構成

| モジュール | 役割 |
|---|---|
| `pocket.py` | fpocketのsubprocessラッパー、Druggability Scoreによるランキング、ポケット残基／中心座標の抽出、ドッキングボックスの計算 |
| `restraints.py` | ポケット中心に基づく調和位置拘束(可動残基集合の計算 + `CustomExternalForce`) |
| `sample.py` | フォースフィールド／水モデルのレジストリ、暗黙的・明示的溶媒でのシステム構築、複数レプリカのrestrained-MDサンプリング |
| `cluster.py` | プールされたトラジェクトリのRMSDクラスタリング、medoid選択、ポケット体積プロキシ／フィルタリング、代表構造・レポートの出力 |
| `visualize.py` | 各クラスタの代表構造のスタンドアロンpy3Dmol HTMLオーバーレイ |
| `pipeline.py` | 各段階のオーケストレーション関数、および `run_end_to_end` |
| `cli.py` | `dd_afpocket-fetch`/`dd_afpocket-pocket`/`dd_afpocket-sample`/`dd_afpocket-cluster`/`dd_afpocket-run` のargparseエントリポイント |
| `progress.py` | 進捗行の出力(レプリカごと、ポケットごと、クラスタごと、ラン中のOpenMMステップレポーター) |
| `parallel.py` | `ProcessPoolExecutor`ベースの並列map。`dd_docking/parallel.py` からコピー |
