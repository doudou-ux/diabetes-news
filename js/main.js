// 主页面JavaScript逻辑
class ArticleManager {
    constructor() {
        this.articles = [];
        this.filteredArticles = [];
        this.currentTag = 'all';
        this.currentView = 'list'; // 'list' 或 'detail'
        this.init();
    }

    async init() {
        await this.loadArticles();
        this.setupEventListeners();
        this.renderArticles();
        this.renderTags();
    }

    async loadArticles() {
        try {
            console.log('开始加载文章数据...');
            
            // 首先尝试从localStorage加载
            const localData = localStorage.getItem('articles');
            if (localData) {
                console.log('从localStorage加载文章数据');
                this.articles = JSON.parse(localData);
                this.filteredArticles = [...this.articles];
                console.log('从localStorage加载的文章数量:', this.articles.length);
                return;
            }
            
            // 如果localStorage没有数据，尝试从文件加载
            const response = await fetch('data/articles.json');
            console.log('响应状态:', response.status, response.statusText);
            
            if (response.ok) {
                const data = await response.json();
                console.log('成功加载文章数据:', data);
                this.articles = data;
                // 保存到localStorage
                localStorage.setItem('articles', JSON.stringify(data));
            } else {
                console.error('加载文章失败:', response.status, response.statusText);
                this.articles = this.getDefaultArticles();
            }
        } catch (error) {
            console.error('加载文章时出错:', error);
            // 如果加载失败，使用默认文章数据
            this.articles = this.getDefaultArticles();
        }
        this.filteredArticles = [...this.articles];
        console.log('最终文章数量:', this.articles.length);
    }

    getDefaultArticles() {
        // 返回默认的示例文章
        return [
            {
                "id": 1752666063384,
                "title": "构建有效的AI智能体",
                "url": "",
                "description": "分享构建AI智能体的实用经验和最佳实践",
                "content": "# 构建有效的AI智能体\n\n我们与数十个跨行业构建LLM智能体的团队合作过。一致的发现是，最成功的实现使用的是简单、可组合的模式，而不是复杂的框架。\n\n## 什么是智能体？\n\n\"智能体\"可以有多种定义。一些客户将智能体定义为完全自主的系统，能够在较长时间内独立运行，使用各种工具来完成复杂任务。\n\n## 何时使用智能体\n\n在使用LLM构建应用程序时，我们建议寻找尽可能简单的解决方案，只有在需要时才增加复杂性。\n\n智能体系统通常会以延迟和成本换取更好的任务性能，您应该考虑这种权衡何时有意义。",
                "tags": ["AI", "智能体", "技术"],
                "rating": 5,
                "date": "2025-07-16T11:41:03.384Z"
            }
        ];
    }

    setupEventListeners() {
        // 搜索功能
        const searchInput = document.getElementById('searchInput');
        searchInput.addEventListener('input', (e) => {
            this.filterArticles(e.target.value, this.currentTag);
        });

        // 标签过滤
        document.addEventListener('click', (e) => {
            if (e.target.classList.contains('tag-filter')) {
                document.querySelectorAll('.tag-filter').forEach(btn => 
                    btn.classList.remove('active')
                );
                e.target.classList.add('active');
                this.currentTag = e.target.dataset.tag;
                this.filterArticles(searchInput.value, this.currentTag);
            }
            
            // 文章卡片点击事件
            if (e.target.closest('.article-card')) {
                const articleId = parseInt(e.target.closest('.article-card').dataset.id);
                this.showArticleDetail(articleId);
            }
        });

        // 返回列表按钮
        const backBtn = document.getElementById('backToList');
        backBtn.addEventListener('click', () => {
            this.showArticleList();
        });
    }

    filterArticles(searchTerm = '', tag = 'all') {
        this.filteredArticles = this.articles.filter(article => {
            const matchesSearch = !searchTerm || 
                article.title.toLowerCase().includes(searchTerm.toLowerCase()) ||
                article.description.toLowerCase().includes(searchTerm.toLowerCase()) ||
                article.tags.some(t => t.toLowerCase().includes(searchTerm.toLowerCase()));
            
            const matchesTag = tag === 'all' || article.tags.includes(tag);
            
            return matchesSearch && matchesTag;
        });
        
        this.renderArticles();
    }

    renderArticles() {
        const container = document.getElementById('articlesContainer');
        const noArticles = document.getElementById('noArticles');
        const loading = document.getElementById('loadingArticles');
        
        // 隐藏加载状态
        loading.style.display = 'none';
        
        if (this.filteredArticles.length === 0) {
            container.style.display = 'none';
            noArticles.style.display = 'block';
            return;
        }
        
        container.style.display = 'grid';
        noArticles.style.display = 'none';
        
        container.innerHTML = this.filteredArticles.map(article => `
            <div class="article-card" data-id="${article.id}" style="cursor: pointer;">
                <h3 class="article-title">${article.title}</h3>
                <p class="article-description">${article.description}</p>
                <div class="article-meta">
                    <span class="article-rating">${this.getRatingStars(article.rating)}</span>
                    <span class="article-date">${this.formatDate(article.date)}</span>
                </div>
                <div class="article-tags">
                    ${article.tags.map(tag => `<span class="tag">${tag}</span>`).join('')}
                </div>
                <div style="margin-top: 1rem; color: #667eea; font-size: 0.9rem;">
                    点击阅读全文 →
                </div>
            </div>
        `).join('');
    }

    showArticleDetail(articleId) {
        const article = this.articles.find(a => a.id === articleId);
        if (!article) return;

        this.currentView = 'detail';
        
        // 隐藏列表视图，显示详情视图
        document.getElementById('articlesListView').style.display = 'none';
        document.getElementById('articleDetailView').style.display = 'block';
        
        // 隐藏搜索和过滤功能
        document.querySelector('.search-box').style.display = 'none';
        document.querySelector('.filter-tags').style.display = 'none';
        
        // 渲染文章详情
        const contentContainer = document.getElementById('articleContent');
        
        // 使用marked.js解析markdown
        let htmlContent;
        if (typeof marked !== 'undefined' && marked.parse) {
            htmlContent = marked.parse(article.content);
        } else {
            console.warn('marked.js未加载，使用简单的文本转换');
            // 简单的markdown转换
            htmlContent = this.simpleMarkdownToHtml(article.content);
        }
        
        contentContainer.innerHTML = `
            <div class="article-meta-detail">
                <div class="meta-item">
                    <span>📅 ${this.formatDate(article.date)}</span>
                </div>
                <div class="meta-item">
                    <span>${this.getRatingStars(article.rating)}</span>
                </div>
                ${article.url ? `<div class="meta-item">
                    <a href="${article.url}" target="_blank" class="original-link">🔗 查看原文</a>
                </div>` : ''}
            </div>
            <div class="article-tags" style="margin-bottom: 2rem;">
                ${article.tags.map(tag => `<span class="tag">${tag}</span>`).join('')}
            </div>
            ${htmlContent}
        `;
        
        // 滚动到顶部
        window.scrollTo(0, 0);
    }

    showArticleList() {
        this.currentView = 'list';
        
        // 显示列表视图，隐藏详情视图
        document.getElementById('articlesListView').style.display = 'block';
        document.getElementById('articleDetailView').style.display = 'none';
        
        // 显示搜索和过滤功能
        document.querySelector('.search-box').style.display = 'block';
        document.querySelector('.filter-tags').style.display = 'block';
    }

    renderTags() {
        const allTags = [...new Set(this.articles.flatMap(article => article.tags))];
        const filterContainer = document.querySelector('.filter-tags');
        
        // 保留"全部"按钮，添加其他标签
        const existingButtons = filterContainer.innerHTML;
        const tagButtons = allTags.map(tag => 
            `<button class="tag-filter" data-tag="${tag}">${tag}</button>`
        ).join('');
        
        filterContainer.innerHTML = existingButtons + tagButtons;
    }

    getRatingStars(rating) {
        const stars = '⭐'.repeat(rating);
        return stars;
    }

    formatDate(dateString) {
        const date = new Date(dateString);
        return date.toLocaleDateString('zh-CN');
    }

    simpleMarkdownToHtml(markdown) {
        console.log('使用增强的Markdown转换');
        
        let html = markdown;
        
        // 0. 预处理 - 修复HTML实体编码
        html = html.replace(/&quot;/g, '"')
                  .replace(/&#34;/g, '"')
                  .replace(/&lt;/g, '<')
                  .replace(/&gt;/g, '>')
                  .replace(/&amp;/g, '&')
                  .replace(/&nbsp;/g, ' ');
        
        // 1. 先处理代码块（避免被其他规则影响）
        html = html.replace(/```[\s\S]*?```/g, (match) => {
            return match.replace(/```(\w*)\n?([\s\S]*?)```/g, '<pre><code class="language-$1">$2</code></pre>');
        });
        
        // 2. 处理行内代码（更严格的匹配）
        html = html.replace(/`([^`\n]+?)`/g, '<code>$1</code>');
        
        // 3. 处理图片（更宽松的匹配，支持更多格式）
        // 处理标准Markdown图片
        html = html.replace(/!\[([^\]]*)\]\(([^)]+)\)/g, '<img src="$2" alt="$1" loading="lazy">');
        
        // 处理可能的图片链接格式
        html = html.replace(/\[图片\]\(([^)]+)\)/g, '<img src="$1" alt="图片" loading="lazy">');
        html = html.replace(/\[image\]\(([^)]+)\)/gi, '<img src="$1" alt="图片" loading="lazy">');
        
        // 处理直接的图片URL（如果单独一行）
        html = html.replace(/^(https?:\/\/[^\s]+\.(jpg|jpeg|png|gif|webp|svg))$/gim, '<img src="$1" alt="图片" loading="lazy">');
        
        // 处理特殊的图片链接格式（如你截图中的情况）
        html = html.replace(/\[图片\]\s*\(([^)]+)\)/g, '<img src="$1" alt="图片" loading="lazy">');
        html = html.replace(/【图片】\s*\(([^)]+)\)/g, '<img src="$1" alt="图片" loading="lazy">');
        
        // 处理可能被错误识别为链接的图片URL
        html = html.replace(/\[(https?:\/\/[^\]]+\.(jpg|jpeg|png|gif|webp|svg)[^\]]*)\]\([^)]*\)/gi, '<img src="$1" alt="图片" loading="lazy">');
        
        // 4. 处理标题（按级别从高到低）
        html = html.replace(/^#### (.*$)/gim, '<h4>$1</h4>');
        html = html.replace(/^### (.*$)/gim, '<h3>$1</h3>');
        html = html.replace(/^## (.*$)/gim, '<h2>$1</h2>');
        html = html.replace(/^# (.*$)/gim, '<h1>$1</h1>');
        
        // 5. 处理粗体和斜体
        html = html.replace(/\*\*\*(.*?)\*\*\*/g, '<strong><em>$1</em></strong>');
        html = html.replace(/\*\*(.*?)\*\*/g, '<strong>$1</strong>');
        html = html.replace(/\*(.*?)\*/g, '<em>$1</em>');
        
        // 6. 处理链接（但先检查是否是图片链接）
        // 先处理可能是图片但被识别为链接的情况
        html = html.replace(/\[([^\]]*)\]\((https?:\/\/[^\s)]+)\)/g, (match, text, url) => {
            // 检查URL是否指向图片
            if (url.match(/\.(jpg|jpeg|png|gif|webp|svg|bmp)(\?[^)]*)?$/i)) {
                return `<img src="${url}" alt="${text || '图片'}" loading="lazy">`;
            }
            // 检查文本是否暗示这是图片
            if (text.match(/(图片|image|图|img)/i) || text === '') {
                return `<img src="${url}" alt="${text || '图片'}" loading="lazy">`;
            }
            // 否则作为普通链接处理
            return `<a href="${url}" target="_blank" rel="noopener noreferrer">${text}</a>`;
        });
        
        // 7. 处理引用
        html = html.replace(/^> (.+)$/gim, '<blockquote><p>$1</p></blockquote>');
        
        // 8. 处理水平线
        html = html.replace(/^---+$/gim, '<hr>');
        html = html.replace(/^\*\*\*+$/gim, '<hr>');
        
        // 9. 处理列表
        const lines = html.split('\n');
        let inList = false;
        let listType = '';
        let processedLines = [];
        
        for (let i = 0; i < lines.length; i++) {
            const line = lines[i];
            const unorderedMatch = line.match(/^[\s]*[\*\-\+•]\s+(.+)$/);
            const orderedMatch = line.match(/^[\s]*\d+\.\s+(.+)$/);
            
            if (unorderedMatch) {
                if (!inList || listType !== 'ul') {
                    if (inList) processedLines.push(`</${listType}>`);
                    processedLines.push('<ul>');
                    inList = true;
                    listType = 'ul';
                }
                processedLines.push(`<li>${unorderedMatch[1]}</li>`);
            } else if (orderedMatch) {
                if (!inList || listType !== 'ol') {
                    if (inList) processedLines.push(`</${listType}>`);
                    processedLines.push('<ol>');
                    inList = true;
                    listType = 'ol';
                }
                processedLines.push(`<li>${orderedMatch[1]}</li>`);
            } else {
                if (inList) {
                    processedLines.push(`</${listType}>`);
                    inList = false;
                    listType = '';
                }
                processedLines.push(line);
            }
        }
        
        if (inList) {
            processedLines.push(`</${listType}>`);
        }
        
        html = processedLines.join('\n');
        
        // 10. 处理段落
        const paragraphs = html.split(/\n\s*\n/);
        html = paragraphs.map(p => {
            p = p.trim();
            if (!p) return '';
            
            // 如果已经是HTML标签，不要包装
            if (p.match(/^<(h[1-6]|ul|ol|pre|blockquote|hr|div|img)/i)) {
                return p;
            }
            
            // 处理单个换行
            p = p.replace(/\n/g, '<br>');
            return `<p>${p}</p>`;
        }).join('\n\n');
        
        // 11. 清理多余的空行
        html = html.replace(/\n{3,}/g, '\n\n');
        
        // 12. 最后清理 - 修复可能的嵌套问题
        html = html.replace(/<p><img/g, '<img')
                  .replace(/><\/p>/g, '>');
        
        console.log('增强Markdown转换完成');
        return html;
    }
}

// 页面加载完成后初始化
document.addEventListener('DOMContentLoaded', () => {
    new ArticleManager();
});