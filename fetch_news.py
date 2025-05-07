# fetch_news.py
import datetime
import html
import os
import time # 用于日期转换
import requests
import feedparser # 用于解析RSS feeds
from bs4 import BeautifulSoup # 用于清理HTML标签

# --- 配置资讯分类 ---
# 关键字将用于构建 Google News RSS 查询
CATEGORIES_CONFIG = {
    "最新研究": {"keywords": "糖尿病 最新研究", "emoji": "🔬"},
    "治疗进展": {"keywords": "糖尿病 治疗进展", "emoji": "💊"},
    "预防与生活方式": {"keywords": "糖尿病 预防 生活方式", "emoji": "🏃‍♀️"},
    "并发症管理": {"keywords": "糖尿病 并发症管理", "emoji": "🩺"},
    "饮食与营养": {"keywords": "糖尿病 饮食营养", "emoji": "🥗"}
}

# --- 帮助函数：判断日期是否在本周 ---
def is_this_week_rss(time_struct, today_date_obj):
    """
    判断给定的 time_struct (来自 feedparser) 是否在本周 (周一到周日)。
    today_date_obj 是今天的 datetime.date 对象。
    """
    if not time_struct:
        return False
    try:
        # feedparser 返回的 time_struct 是 time.struct_time 对象
        article_date = datetime.date(time_struct.tm_year, time_struct.tm_mon, time_struct.tm_mday)
        
        # 计算本周的开始 (周一) 和结束 (周日)
        start_of_week = today_date_obj - datetime.timedelta(days=today_date_obj.weekday())
        end_of_week = start_of_week + datetime.timedelta(days=6)
        
        return start_of_week <= article_date <= end_of_week
    except Exception as e:
        print(f"    [is_this_week_rss] 日期转换错误: {e} - Time Struct: {time_struct}")
        return False # 如果日期无效或解析失败

# --- 帮助函数：清理 HTML ---
def clean_html(raw_html):
    """使用 BeautifulSoup 清理 HTML 标签"""
    if not raw_html:
        return ""
    try:
        soup = BeautifulSoup(raw_html, "html.parser")
        return soup.get_text()
    except Exception as e:
        print(f"    [clean_html] HTML 清理错误: {e}")
        return raw_html # 出错则返回原始文本

# --- 从 Google News RSS 获取真实新闻 ---
def fetch_real_news_from_google_rss(category_name, keywords_for_rss):
    """
    从 Google News RSS feed 获取指定分类的新闻。
    """
    print(f"  正在为分类 '{category_name}' (关键词: '{keywords_for_rss}') 获取真实新闻...")
    # 构建 Google News RSS URL (中文新闻, 中国地区)
    # 注意：Google News RSS 的确切行为和URL参数可能会变化
    base_url = "https://news.google.com/rss/search"
    query_params = {
        "q": keywords_for_rss,
        "hl": "zh-CN", # 语言：简体中文
        "gl": "CN",    # 地理位置：中国
        "ceid": "CN:zh-Hans" # 版本标识
    }
    
    articles_this_week = []
    today = datetime.date.today() # 获取今天的日期对象，用于 is_this_week_rss

    try:
        headers = { # 添加 User-Agent，模拟浏览器访问
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
        }
        response = requests.get(base_url, params=query_params, headers=headers, timeout=15) # 设置15秒超时
        response.raise_for_status() # 如果请求失败 (如 404, 500), 会抛出异常

        feed = feedparser.parse(response.content)
        
        if not feed.entries:
            print(f"    未找到分类 '{category_name}' 的新闻条目。")
            return []

        print(f"    分类 '{category_name}' 原始获取到 {len(feed.entries)} 条新闻，开始筛选本周新闻...")

        for entry in feed.entries:
            title = entry.get("title", "无标题")
            link = entry.get("link", "javascript:void(0);")
            
            # 获取发布日期，优先使用 published_parsed，其次 updated_parsed
            published_time_struct = entry.get("published_parsed")
            if not published_time_struct:
                published_time_struct = entry.get("updated_parsed")

            # 使用 is_this_week_rss 进行筛选
            if is_this_week_rss(published_time_struct, today):
                summary_html = entry.get("summary", "暂无摘要")
                snippet = clean_html(summary_html) # 清理HTML标签
                
                source_info = entry.get("source")
                source_name = source_info.get("title", "Google News") if source_info else "Google News"

                # 将解析的日期转换为更易读的格式 (可选, HTML中可能不需要这么精确)
                time_display_str = "未知时间"
                if published_time_struct:
                    try:
                        time_display_str = time.strftime("%Y-%m-%d", published_time_struct)
                    except:
                        pass # 保持 "未知时间"

                articles_this_week.append({
                    "title": title,
                    "url": link,
                    "snippet": snippet,
                    "source": source_name,
                    "time": time_display_str # 或者使用更友好的 "X天前" 等，但需要额外逻辑
                })
                if len(articles_this_week) >= 10: # 每个分类最多获取10条本周新闻
                    break 
        
        print(f"    分类 '{category_name}' 筛选后得到 {len(articles_this_week)} 条本周新闻。")

    except requests.exceptions.RequestException as e:
        print(f"    获取分类 '{category_name}' 新闻时发生网络错误: {e}")
    except Exception as e:
        print(f"    处理分类 '{category_name}' 新闻时发生未知错误: {e}")
        
    return articles_this_week


# --- HTML 生成逻辑 (与之前版本基本一致) ---
def generate_html_content(all_news_data):
    """根据新闻数据生成完整的HTML页面内容"""
    current_time_str = datetime.datetime.now().strftime('%Y年%m月%d日 %H:%M:%S')
    app_timezone = os.getenv('APP_TIMEZONE', 'UTC') # 从环境变量读取时区
    if app_timezone != 'UTC':
        # 尝试使用环境变量设置的时间，如果失败则回退
        try:
            # 注意：这部分依赖于运行环境对时区的支持。
            # 在GitHub Actions中，设置TZ环境变量通常有效。
             current_time_str = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).strftime('%Y年%m月%d日 %H:%M:%S %Z') # 假设东八区
             if app_timezone == "Asia/Shanghai": # 更精确的示例
                 pass # 实际应用中可能需要 pytz 或 dateutil
        except Exception as e:
            print(f"应用时区 ({app_timezone}) 时出错: {e}。将使用默认服务器时间。")
            current_time_str = datetime.datetime.now().strftime('%Y年%m月%d日 %H:%M:%S (服务器时间)')


    current_year = datetime.datetime.now().year
    github_repo_url = f"https://github.com/{os.getenv('GITHUB_REPOSITORY', 'YOUR_USERNAME/YOUR_REPOSITORY_NAME')}"
    if github_repo_url == "https://github.com/YOUR_USERNAME/YOUR_REPOSITORY_NAME" and os.getenv('GITHUB_ACTOR'):
        # 尝试从 GITHUB_ACTOR 和仓库名推断 (不一定准确)
        actor = os.getenv('GITHUB_ACTOR')
        repo_name = os.getenv('GITHUB_REPOSITORY_NAME', 'diabetes-news') # GITHUB_REPOSITORY_NAME 可能不存在
        if actor and repo_name:
            github_repo_url = f"https://github.com/{actor}/{repo_name}"
        elif os.getenv('GITHUB_REPOSITORY'): # GITHUB_REPOSITORY 格式是 owner/repo
             github_repo_url = f"https://github.com/{os.getenv('GITHUB_REPOSITORY')}"


    html_output = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>糖尿病前沿资讯</title>
    <script src="https://cdn.tailwindcss.com"></script>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&display=swap" rel="stylesheet">
    <style>
        body {{ font-family: 'Inter', sans-serif; background-color: #f0f4f8; }}
        .news-card {{ transition: transform 0.3s ease, box-shadow 0.3s ease; display: flex; flex-direction: column; justify-content: space-between; height: 100%; background-color: #ffffff; border-radius: 0.5rem; }}
        .news-card:hover {{ transform: translateY(-5px); box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05); }}
        .category-section {{ background-color: #ffffff; border-radius: 0.75rem; box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06); padding: 1.25rem; }}
        @media (min-width: 768px) {{ .category-section {{ padding: 2rem; }} }}
        .category-title-text {{ border-bottom: 3px solid #3b82f6; padding-bottom: 0.75rem; margin-bottom: 1.75rem; color: #1e40af; font-size: 1.5rem; }}
        @media (min-width: 768px) {{ .category-title-text {{ font-size: 1.875rem; }} }}
        .news-item-link {{ color: #2563eb; text-decoration: none; font-weight: 600; }}
        .news-item-link:hover {{ text-decoration: underline; color: #1d4ed8; }}
        .source-tag, .time-tag {{ padding: 0.25rem 0.75rem; border-radius: 0.375rem; font-size: 0.75rem; font-weight: 500; line-height: 1.25; }}
        .source-tag {{ background-color: #e0e7ff; color: #3730a3; }}
        .time-tag {{ background-color: #d1fae5; color: #065f46; }}
        .card-content {{ flex-grow: 1; }}
        .card-footer {{ margin-top: auto; }}
        .header-main-title {{ font-size: 1.875rem; }}
        @media (min-width: 640px) {{ .header-main-title {{ font-size: 2.25rem; }} }}
        @media (min-width: 768px) {{ .header-main-title {{ font-size: 3rem; }} }}
        .main-container {{ padding: 1rem; }}
        @media (min-width: 768px) {{ .main-container {{ padding: 2rem; }} }}
        .loader {{ border: 5px solid #f3f3f3; border-radius: 50%; border-top: 5px solid #3498db; width: 50px; height: 50px; -webkit-animation: spin 1s linear infinite; animation: spin 1s linear infinite; margin: 20px auto; }}
        @-webkit-keyframes spin {{ 0% {{ -webkit-transform: rotate(0deg); }} 100% {{ -webkit-transform: rotate(360deg); }} }}
        @keyframes spin {{ 0% {{ transform: rotate(0deg); }} 100% {{ transform: rotate(360deg); }} }}
    </style>
</head>
<body class="bg-gray-100 text-gray-800">
    <div class="container mx-auto main-container">
        <header class="text-center mb-10 md:mb-16">
            <h1 class="font-bold text-blue-700 header-main-title">糖尿病前沿资讯</h1>
            <p class="text-gray-600 mt-3 text-base md:text-lg">本周最新动态（自动更新于：<span id="updateTime">{current_time_str}</span>）</p>
            <p class="text-sm text-gray-500 mt-2">资讯来源：Google News RSS Feeds</p>
        </header>
        <div id="news-content" class="space-y-12">
            <div id="loading-indicator" class="text-center py-10">
                 <div class="loader"></div>
                 <p class="text-gray-600 mt-2">正在加载最新资讯...</p>
            </div>
        </div>
    </div>
    <footer class="text-center p-6 mt-12 text-gray-600 text-sm border-t border-gray-300">
        <p>&copy; {current_year} 糖尿病资讯聚合. <a href="{github_repo_url}" target="_blank" rel="noopener noreferrer" class="text-blue-600 hover:underline">项目源码</a></p>
        <p class="mt-1">本站内容仅供参考, 不构成医疗建议。</p>
    </footer>
    <script>
        // 简单的JS来在实际内容加载后移除加载指示器 (如果HTML是静态生成，这个可能不是必须的)
        // 但如果Python脚本生成内容后，我们想用JS动态插入，则有用
        // 在当前纯静态生成HTML的模式下，这个脚本块可以简化或移除
        document.addEventListener('DOMContentLoaded', function() {{
            const newsContent = document.getElementById('news-content');
            const loadingIndicator = document.getElementById('loading-indicator');
            // 检查 news-content 是否已经有除了 loading-indicator 之外的子元素
            let hasRealContent = false;
            for (let i = 0; i < newsContent.childNodes.length; i++) {{
                if (newsContent.childNodes[i].id !== 'loading-indicator' && newsContent.childNodes[i].nodeType === 1) {{
                    hasRealContent = true;
                    break;
                }}
            }}
            if (loadingIndicator && hasRealContent) {{
                // loadingIndicator.style.display = 'none'; // 或者直接移除
            }} else if (loadingIndicator && !hasRealContent && newsContent.textContent.includes("未能加载")) {
                // 如果只有错误信息，也移除加载动画
                loadingIndicator.style.display = 'none';
            }
            // 更新时间戳 (如果需要更动态的时间更新)
            // document.getElementById('updateTime').textContent = new Date().toLocaleString('zh-CN', {{ hour12: false }});
        }});
    </script>
</body>
</html>"""

    # 生成每个分类的新闻内容
    news_html_parts = []
    found_any_news = False

    if not all_news_data or all(not articles for articles in all_news_data.values()):
        news_html_parts.append('<p class="text-center text-gray-500 text-xl py-10">抱歉，目前未能加载到本周相关的糖尿病资讯。</p>')
    else:
        for category, articles in all_news_data.items():
            category_emoji = CATEGORIES_CONFIG.get(category, {}).get("emoji", "")
            category_html = f"""
            <section class="category-section">
                <h2 class="font-semibold category-title-text">{category_emoji} {html.escape(category)}</h2>
            """
            if not articles:
                category_html += '<p class="text-gray-500">本周暂无该分类下的资讯。</p>'
            else:
                found_any_news = True
                category_html += '<div class="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-6 md:gap-8 news-grid">'
                for article in articles: # 已经在获取时限制了数量
                    title = html.escape(article.get('title', '无标题'))
                    url = html.escape(article.get('url', 'javascript:void(0);'))
                    snippet_raw = article.get('snippet', '暂无摘要')
                    # 限制摘要长度，避免过长
                    snippet = html.escape(snippet_raw[:150] + ('...' if len(snippet_raw) > 150 else ''))
                    source = html.escape(article.get('source', '未知来源'))
                    time_display = html.escape(article.get('time', '未知时间'))

                    category_html += f"""
                    <div class="news-card shadow-md hover:shadow-lg p-5">
                        <div class="card-content">
                            <h3 class="text-lg mb-2">
                                <a href="{url}" {'target="_blank" rel="noopener noreferrer"' if url != 'javascript:void(0);' else ''} class="news-item-link">
                                    {title}
                                </a>
                            </h3>
                            <p class="text-gray-600 text-sm mb-4 leading-relaxed">{snippet}</p>
                        </div>
                        <div class="card-footer text-xs text-gray-500 flex flex-wrap gap-2 items-center pt-3 border-t border-gray-200">
                            <span class="source-tag">{source}</span>
                            <span class="time-tag">{time_display}</span>
                        </div>
                    </div>
                    """
                category_html += "</div>" # close grid
            category_html += "</section>"
            news_html_parts.append(category_html)
    
    # 替换HTML模板中的占位符
    # 移除加载指示器，如果已经有内容或错误信息
    if found_any_news or "<p class=\"text-center text-gray-500 text-xl py-10\">" in news_html_parts[0]:
        html_output = html_output.replace(
            """<div id="loading-indicator" class="text-center py-10">
                 <div class="loader"></div>
                 <p class="text-gray-600 mt-2">正在加载最新资讯...</p>
            </div>""",
            "".join(news_html_parts)
        )
    else: # 如果什么都没找到，但也不是明确的错误信息，保留加载（理论上不应发生）
        html_output = html_output.replace(
             """<div id="news-content" class="space-y-12">
            <div id="loading-indicator" class="text-center py-10">
                 <div class="loader"></div>
                 <p class="text-gray-600 mt-2">正在加载最新资讯...</p>
            </div>
        </div>""",
             f"""<div id="news-content" class="space-y-12">{"".join(news_html_parts)}</div>"""
        )


    return html_output

# --- 主执行逻辑 ---
if __name__ == "__main__":
    print("开始从 Google News RSS 生成糖尿病资讯网页...")
    all_news_data_for_html = {}

    for category_name_zh, config in CATEGORIES_CONFIG.items():
        # 使用真实新闻获取函数
        articles = fetch_real_news_from_google_rss(category_name_zh, config["keywords"])
        all_news_data_for_html[category_name_zh] = articles
        print(f"  分类 '{category_name_zh}' 处理完毕，获取到 {len(articles)} 条本周新闻。")
        time.sleep(1) # 友好访问，避免请求过于频繁，在分类间稍作停顿
    
    final_html = generate_html_content(all_news_data_for_html)
    
    output_filename = "index.html" 
    try:
        with open(output_filename, "w", encoding="utf-8") as f:
            f.write(final_html)
        print(f"成功生成网页：{output_filename}")
    except IOError as e:
        print(f"错误：无法写入文件 {output_filename}。错误信息: {e}")
    except Exception as e:
        print(f"生成过程中发生未知错误: {e}")

    print("资讯网页生成完毕。")

