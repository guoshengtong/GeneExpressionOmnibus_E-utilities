# GitHub 上传指南 | GitHub Upload Guide

## 📋 当前状态 | Current Status

✅ Git 仓库已初始化  
✅ 所有文件已提交到本地仓库  
✅ 准备上传到 GitHub

## 🚀 上传步骤 | Upload Steps

### 方法 1: 使用 GitHub CLI (推荐) | Using GitHub CLI (Recommended)

如果您已安装 GitHub CLI:

```bash
# 创建仓库并推送
gh repo create GeneExpressionOmnibus_E-utilities --public --source=. --remote=origin --push
```

### 方法 2: 在 GitHub 网站创建仓库 | Create Repository on GitHub Website

#### 步骤 1: 创建新仓库 | Step 1: Create New Repository

1. 访问 https://github.com/new
2. 填写仓库信息：
   - **Repository name**: `GeneExpressionOmnibus_E-utilities`
   - **Description**: `GEO Lung Metastasis Data Mining Pipeline - Automated metadata mining and filtering pipeline for lung metastasis samples from NCBI GEO database`
   - **Visibility**: 选择 Public 或 Private
   - **不要**勾选 "Initialize this repository with a README"（我们已经有了）
3. 点击 "Create repository"

#### 步骤 2: 连接远程仓库并推送 | Step 2: Connect Remote and Push

在终端中运行以下命令（将 `YOUR_USERNAME` 替换为您的 GitHub 用户名）：

```bash
cd /Users/khugjil-devstation/Projects/GeneExpressionOmnibus_E-utilities

# 添加远程仓库
git remote add origin https://github.com/YOUR_USERNAME/GeneExpressionOmnibus_E-utilities.git

# 或者使用 SSH（如果您配置了 SSH key）
# git remote add origin git@github.com:YOUR_USERNAME/GeneExpressionOmnibus_E-utilities.git

# 推送代码到 GitHub
git branch -M main
git push -u origin main
```

#### 步骤 3: 验证 | Step 3: Verify

访问您的 GitHub 仓库页面，确认所有文件都已上传。

## 🔐 认证问题 | Authentication Issues

如果遇到认证问题，您可能需要：

### 使用 Personal Access Token (推荐)

1. 访问 https://github.com/settings/tokens
2. 点击 "Generate new token (classic)"
3. 选择权限：至少勾选 `repo`
4. 生成 token 并复制
5. 推送时使用 token 作为密码：

```bash
git push -u origin main
# Username: 您的GitHub用户名
# Password: 粘贴您的token
```

### 或配置 SSH Key

```bash
# 生成 SSH key（如果还没有）
ssh-keygen -t ed25519 -C "your_email@example.com"

# 添加到 ssh-agent
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519

# 复制公钥到剪贴板
cat ~/.ssh/id_ed25519.pub | pbcopy

# 在 GitHub 设置中添加 SSH key: https://github.com/settings/keys
```

## 📝 后续更新 | Future Updates

上传后，每次更新代码时：

```bash
# 添加更改
git add .

# 提交更改
git commit -m "描述您的更改"

# 推送到 GitHub
git push
```

## ✅ 检查清单 | Checklist

- [ ] GitHub 仓库已创建
- [ ] 远程仓库已连接
- [ ] 代码已成功推送
- [ ] README.md 在 GitHub 上正确显示
- [ ] 所有文件都在正确的目录中

## 🆘 需要帮助？| Need Help?

如果遇到问题，请检查：
- Git 是否正确安装：`git --version`
- 是否已登录 GitHub
- 网络连接是否正常
- 仓库名称是否正确

---

**提示**: 上传完成后，您可以在 GitHub 仓库页面添加：
- Topics/标签（如：`bioinformatics`, `genomics`, `data-mining`, `python`）
- 仓库描述
- 许可证文件（如果需要）
