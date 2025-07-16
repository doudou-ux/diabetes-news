# 部署指南 - 文章收录系统

## 🎯 推荐部署方案：GitHub Pages

### 为什么选择GitHub Pages？
- ✅ 完全免费
- ✅ 自动部署
- ✅ 稳定可靠
- ✅ 支持自定义域名
- ✅ 适合静态网站

## 📋 部署步骤

### 第一步：准备GitHub账号
1. 访问 [github.com](https://github.com) 注册账号
2. 验证邮箱

### 第二步：创建仓库
1. 点击右上角 "+" → "New repository"
2. 仓库名建议：`article-collection` 或 `my-articles`
3. 设置为 Public（公开）
4. 勾选 "Add a README file"
5. 点击 "Create repository"

### 第三步：上传文件
**方法1：网页上传（推荐新手）**
1. 在仓库页面点击 "uploading an existing file"
2. 拖拽所有项目文件到页面
3. 写提交信息：`Initial upload`
4. 点击 "Commit changes"

**方法2：Git命令行**
```bash
git clone https://github.com/你的用户名/仓库名.git
cd 仓库名
# 复制所有项目文件到这个文件夹
git add .
git commit -m "Initial upload"
git push
```

### 第四步：启用GitHub Pages
1. 在仓库页面点击 "Settings"
2. 滚动到 "Pages" 部分
3. Source 选择 "Deploy from a branch"
4. Branch 选择 "main"
5. 点击 "Save"

### 第五步：访问网站
- 等待1-2分钟
- 访问：`https://你的用户名.github.io/仓库名`

## 🔄 日常更新流程

### 添加新文章：
1. 本地打开 `admin.html`
2. 添加新文章
3. 下载生成的 `articles.json`
4. 在GitHub仓库中：
   - 点击 `data/articles.json`
   - 点击编辑按钮（铅笔图标）
   - 替换内容
   - 提交更改

### 自动化更新（可选）：
使用 `sync-data.html` 工具：
1. 打开 `sync-data.html`
2. 粘贴新的文章数据
3. 点击同步
4. 手动更新GitHub上的文件

## 🌐 其他部署选项

### Netlify（拖拽部署）
1. 访问 [netlify.com](https://netlify.com)
2. 注册账号
3. 拖拽项目文件夹到页面
4. 获得网址：`https://随机名称.netlify.app`

### Vercel
1. 访问 [vercel.com](https://vercel.com)
2. 用GitHub账号登录
3. 导入你的GitHub仓库
4. 自动部署

## 💡 更新建议

### 每次更新文章时：
1. **本地编辑** - 使用管理后台添加文章
2. **测试预览** - 确保内容正确
3. **更新线上** - 替换GitHub上的articles.json文件
4. **等待部署** - GitHub Pages会自动更新（1-2分钟）

### 批量更新技巧：
- 可以一次性添加多篇文章
- 使用"清理内容"和"转换图片链接"工具优化
- 定期备份articles.json文件

## 🔧 高级配置

### 自定义域名（可选）：
1. 购买域名（如：yourdomain.com）
2. 在GitHub Pages设置中添加自定义域名
3. 配置DNS记录

### 搜索引擎优化：
- 在每篇文章中添加合适的标签
- 使用描述性的文章标题
- 定期更新内容

## 📞 技术支持

如果遇到问题：
1. 检查浏览器控制台错误
2. 确认所有文件都已上传
3. 等待GitHub Pages部署完成
4. 清除浏览器缓存

## 🎉 完成！

现在你就有了一个完全免费的在线文章收录网站！
- 网址：`https://你的用户名.github.io/仓库名`
- 管理：通过GitHub更新文件
- 维护：每1-2周更新一次文章