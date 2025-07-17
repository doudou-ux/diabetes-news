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
            },
            {
                "id": 1752666063385,
                "title": "图片显示测试文章",
                "url": "",
                "description": "测试各种图片格式的显示效果",
                "content": "# 图片显示测试文章\n\n这篇文章包含了各种格式的图片链接，用于测试图片修复功能。\n\n## 标准Markdown图片格式\n\n![示例图片](https://picsum.photos/600/400?random=1)\n\n## 中文图片链接格式\n\n[图片](https://picsum.photos/600/400?random=2)\n\n## 可能被识别为普通链接的图片\n\n[https://picsum.photos/600/400?random=3](https://picsum.photos/600/400?random=3)\n\n## 直接的图片URL\n\nhttps://picsum.photos/600/400?random=4\n\n## 其他可能的格式\n\n【图片】(https://picsum.photos/600/400?random=5)\n\n[image](https://picsum.photos/600/400?random=6)\n\n## 普通文本中的图片链接\n\n这里有一个图片链接：[图片](https://picsum.photos/600/400?random=7)，应该被正确显示。\n\n## 代码块中的图片（不应该被转换）\n\n```\n![这个不应该被转换](https://example.com/image.jpg)\n```\n\n## 结论\n\n通过点击\"修复图片显示\"按钮，所有上述格式的图片都应该被正确识别和显示。",
                "tags": ["测试", "图片", "技术"],
                "rating": 4,
                "date": "2025-07-17T10:00:00.000Z"
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

        // 修复图片按钮
        const fixImagesBtn = document.getElementById('fixImages');
        fixImagesBtn.addEventListener('click', () => {
            this.aggressiveImageFix();
            this.showMessage('图片修复完成！', 'success');
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
        
        // 应用激进的图片修复
        setTimeout(() => {
            this.aggressiveImageFix();
        }, 100);
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

    aggressiveImageFix() {
        console.log('开始激进图片修复...');
        const contentContainer = document.getElementById('articleContent');
        if (!contentContainer) return;

        // 1. 查找所有可能的图片链接文本
        const walker = document.createTreeWalker(
            contentContainer,
            NodeFilter.SHOW_TEXT,
            null,
            false
        );

        const textNodes = [];
        let node;
        while (node = walker.nextNode()) {
            textNodes.push(node);
        }

        // 2. 处理文本节点中的图片链接
        textNodes.forEach(textNode => {
            let text = textNode.textContent;
            let hasChanges = false;

            // 检测各种图片链接格式
            const imagePatterns = [
                // 标准格式
                /\[图片\]\((https?:\/\/[^\s)]+)\)/g,
                /\[image\]\((https?:\/\/[^\s)]+)\)/gi,
                /【图片】\((https?:\/\/[^\s)]+)\)/g,
                // 直接URL
                /(https?:\/\/[^\s]+\.(jpg|jpeg|png|gif|webp|svg|bmp)(\?[^\s]*)?)/gi,
                // 其他可能的格式
                /\[([^\]]*)\]\((https?:\/\/[^\s)]+\.(jpg|jpeg|png|gif|webp|svg|bmp)(\?[^\s]*)?)\)/gi
            ];

            imagePatterns.forEach(pattern => {
                if (pattern.test(text)) {
                    text = text.replace(pattern, (match, ...args) => {
                        hasChanges = true;
                        let url, alt = '图片';
                        
                        if (args.length >= 2 && args[1]) {
                            // 格式: [alt](url)
                            alt = args[0] || '图片';
                            url = args[1];
                        } else {
                            // 直接URL
                            url = args[0];
                        }
                        
                        return `<img src="${url}" alt="${alt}" loading="lazy" style="max-width: 100%; height: auto; border-radius: 8px; margin: 1rem 0;">`;
                    });
                    // 重置正则表达式
                    pattern.lastIndex = 0;
                }
            });

            if (hasChanges) {
                // 创建临时容器来解析HTML
                const tempDiv = document.createElement('div');
                tempDiv.innerHTML = text;
                
                // 替换文本节点
                const parent = textNode.parentNode;
                while (tempDiv.firstChild) {
                    parent.insertBefore(tempDiv.firstChild, textNode);
                }
                parent.removeChild(textNode);
            }
        });

        // 3. 查找所有现有的链接，检查是否应该是图片
        const links = contentContainer.querySelectorAll('a');
        links.forEach(link => {
            const href = link.href;
            const text = link.textContent.trim();
            
            // 检查链接是否指向图片
            if (href.match(/\.(jpg|jpeg|png|gif|webp|svg|bmp)(\?.*)?$/i) || 
                text.match(/(图片|image|图|img)/i)) {
                
                const img = document.createElement('img');
                img.src = href;
                img.alt = text || '图片';
                img.loading = 'lazy';
                img.style.cssText = 'max-width: 100%; height: auto; border-radius: 8px; margin: 1rem 0; display: block;';
                
                link.parentNode.replaceChild(img, link);
                console.log('转换链接为图片:', href);
            }
        });

        // 4. 处理可能的图片加载错误
        const images = contentContainer.querySelectorAll('img');
        images.forEach(img => {
            // 添加错误处理
            img.onerror = function() {
                console.log('图片加载失败:', this.src);
                this.style.display = 'block';
                this.style.padding = '2rem';
                this.style.background = '#f8f9fa';
                this.style.border = '2px dashed #dee2e6';
                this.style.borderRadius = '8px';
                this.style.color = '#6c757d';
                this.style.textAlign = 'center';
                this.style.fontStyle = 'italic';
                this.alt = `🖼️ 图片加载失败: ${this.alt || '图片'}`;
            };

            // 确保图片样式正确
            if (!img.style.maxWidth) {
                img.style.maxWidth = '100%';
                img.style.height = 'auto';
                img.style.borderRadius = '8px';
                img.style.margin = '1rem 0';
                img.style.display = 'block';
            }
        });

        // 5. 添加图片点击放大功能
        this.addImageClickHandler();

        console.log('激进图片修复完成，处理了', images.length, '张图片');
    }

    addImageClickHandler() {
        const images = document.querySelectorAll('#articleContent img');
        images.forEach(img => {
            img.style.cursor = 'pointer';
            img.title = '点击查看大图';
            
            img.addEventListener('click', (e) => {
                e.preventDefault();
                this.showImageModal(img.src, img.alt);
            });
        });
    }

    showImageModal(src, alt) {
        // 创建模态框
        const modal = document.createElement('div');
        modal.style.cssText = `
            position: fixed;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            background: rgba(0,0,0,0.9);
            display: flex;
            justify-content: center;
            align-items: center;
            z-index: 10000;
            cursor: pointer;
        `;

        const img = document.createElement('img');
        img.src = src;
        img.alt = alt;
        img.style.cssText = `
            max-width: 90%;
            max-height: 90%;
            object-fit: contain;
            border-radius: 8px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.5);
        `;

        modal.appendChild(img);
        document.body.appendChild(modal);

        // 点击关闭
        modal.addEventListener('click', () => {
            document.body.removeChild(modal);
        });

        // ESC键关闭
        const handleEsc = (e) => {
            if (e.key === 'Escape') {
                document.body.removeChild(modal);
                document.removeEventListener('keydown', handleEsc);
            }
        };
        document.addEventListener('keydown', handleEsc);
    }

    showMessage(message, type = 'info') {
        // 创建消息元素
        const messageEl = document.createElement('div');
        messageEl.className = `message message-${type}`;
        messageEl.textContent = message;
        
        // 添加到页面
        document.body.appendChild(messageEl);
        
        // 3秒后自动移除
        setTimeout(() => {
            if (document.body.contains(messageEl)) {
                document.body.removeChild(messageEl);
            }
        }, 3000);
        
        // 点击移除
        messageEl.addEventListener('click', () => {
            if (document.body.contains(messageEl)) {
                document.body.removeChild(messageEl);
            }
        });
    }
}

// 页面加载完成后初始化
document.addEventListener('DOMContentLoaded', () => {
    new ArticleManager();
});