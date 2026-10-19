# GitHub準備（個人情報・患者データを持ち込まない）

1. GitHubにログインし、右上 `+` → `New repository` を選択。
2. 例として `dgl-mva-2026-track1` を入力。公開可否は **Private** を選択（後日公開するには公式Rulesに従い審査に支障がないよう変更）。
3. `Create repository` → `uploading an existing file` から、このZIPをPCで展開した中のフォルダ・ファイルをアップロード。**ZIP自体をアップロードしない**（採点者がコードを見られないため）。`README.md` 等も含める。
4. GitHub上のファイル一覧で `notebooks/`, `scripts/`, `templates/`, `reports/`, `docs/` と `.gitignore` の内容を確認。
5. 患者データ、Hugging Faceアクセストークン、秘密鍵、Google Drive/Colabの出力、患者の症状記述・遺伝子型をアップロードしない。`.gitignore` だけでは安全を保証できないので、コミット前に手作業で確認。
6. 現在のURLを控える。競技フォームのGitHub URL欄には、**最終提出の直前**に入力する。リンクの正確性と公開時期を確認する。

## 現在の完了条件

- [ ] GitHub repository URL確定
- [ ] Public code only / controlled data noneを手動確認
- [ ] Methods Reportを実解析結果に合わせて更新（未着手）
- [ ] 実患者VCF解析（適法な計算環境確認後）
- [ ] 科学的根拠のある公式CSV生成・オフライン事前検証
- [ ] ユーザーの明示的な公式提出承認

GitHubのRepository作成は外部公開・契約/権利に関係するため、実施するのは利用者です。
