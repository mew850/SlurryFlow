# 在 VS Code 上传核心成果

本目录已经初始化 Git，当前分支为 `main`；整理前没有提交、没有远程地址。忽略文件的正确名称是 `.gitignore`，不是 `.ignore`。

## 1. 打开项目并检查文件

在 VS Code 选择“文件 → 打开文件夹”，打开本项目根目录。选择“终端 → 新建终端”，使用 PowerShell。

```powershell
git status --short
git remote -v
git ls-files --others --exclude-standard
git add --dry-run .
```

最后两条只预览待上传文件，不上传。原始数据、PDF/PPT/Word、`output/`、`tmp/`、`.codex-finalizer/`、旧版 UDF、压缩包及编译产物应不在清单中。

`.gitignore` 是白名单：以后新增加的脚本、配置或成果，需要显式增加对应 `!` 规则。它不删除本地文件，也不能撤回已经跟踪或已经上传的文件。不要使用 `git add -f` 绕过规则。

## 2. 本地验证（可选运行环境为已有 ml）

按 Ctrl+Shift+P，执行“Python: Select Interpreter”，选择 Conda 的 `ml`。随后在终端执行：

```powershell
conda activate ml
python -m unittest discover -s tests -v
```

已有环境缺依赖时再执行 `python -m pip install -r requirements.txt`。这组 Python 测试不等于 Fluent 编译或设备算例验证。独立 C 数学测试需要可用的 C 编译器。

## 3. 暂存与本地提交

点击左侧“源代码管理”，检查变更文件，然后点击“全部暂存更改”。也可以执行：

```powershell
git add .
git diff --cached --stat
git diff --cached --name-only
git diff --cached --check
git commit -m "Keep core Model C, Fluent lookup UDF and validation results"
```

如果 Git 提示缺少身份信息，仅为本仓库填写自己的名字和 GitHub 邮箱（可以使用 GitHub 提供的 noreply 邮箱），再提交：

```powershell
git config user.name "你的名字"
git config user.email "你的GitHub提交邮箱"
```

## 4. 新建远程仓库并上传（本次选择）

建议仓库名 `MicroRheology-Core`，选择 Private。不连接之前的 MicroRheology 仓库。

1. 按 Ctrl+Shift+G 打开源代码管理。如果本地已有本次核心提交且无变更，跳过第 3 节的重复提交。
2. 点击“Publish Branch / 发布分支”；也可以按 Ctrl+Shift+P 搜索“GitHub: Publish to GitHub”。
3. 按界面提示在浏览器登录自己的 GitHub 账号，完成认证后返回 VS Code。无需在代码或终端中粘贴访问令牌。
4. 输入 `MicroRheology-Core`；如果已存在同名仓库，使用一个新名称。
5. 选择“Publish to GitHub private repository / 发布到 GitHub 私有仓库”。若出现文件选择，核对其符合本项目白名单。
6. 等待推送完成，点击“Open on GitHub / 在 GitHub 上打开”，检查文件列表。

这个操作会创建远程仓库、设置 origin 并推送当前 main 分支。本地 38 个核心文件约 1 MB，无需 Git LFS。

如果界面发布入口不可用，在 GitHub 右上角 + → New repository 创建同名 Private 空仓库，不初始化 README、.gitignore 或 License。复制网页给出的 HTTPS URL，在终端执行（把示例地址替换为网页上的真实地址）：

```powershell
git remote -v
git remote add origin https://github.com/YOUR_ACCOUNT/MicroRheology-Core.git
git push -u origin main
```

如果 origin 已存在，先核对地址，不要盲目重复添加或修改。如果推送失败，检查认证与远程历史，不使用 `--force`。

## 5. 检查上传结果

```powershell
git status -sb
git log -1 --oneline
git ls-files
```

在 GitHub 网页切换到刚推送的分支，检查 README、模型代码、`udf/generated/model_c_table.h`、验证图和指标，并确认没有原始数据和汇报材料。已有仓库其他分支上的文件不会被本次白名单移除。

## 核心文件选择

- `src/`：微观流变、Model C、训练、表格导出、PBM 与工具代码。
- `configs/`：模型参数和 15-bin 设置。
- `artifacts/model_c/`：训练模型、参数摘要、验证图和指标。
- `artifacts/fluent_table/validation.json`：已有抽样插值误差。
- `udf/model_c_viscosity_udf.c`：基本查表 UDF，使用前按源码接口设置输入。
- `udf/model_c_pbm_final.c`、两个 `cb_final_*.h`、`FINAL_SETUP.md`：独立的 Mixture-PBM 工程版本及配置；不与基础版同时挂接。
- `udf/generated/model_c_table.h`：实际编译所需的离散表。
- `tests/`：Python 测试和 C 数学测试源码。

UDF 的 final 文件名仅标识工程版本，不表示所有设备条件已验证。上传时未更改物理模型或参数。
