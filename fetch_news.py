# fetch_news.py
import datetime
import html
import os
import time
import requests
import feedparser
from bs4 import BeautifulSoup

# --- (1) 配置权威 RSS 源 ---
# 请您将来替换 YOUR_PUBMED_RSS_URL 和 YOUR_ADA_JOURNAL_RSS_URL 为实际的链接
# 您可以添加更多类似的条目
AUTHORITATIVE_RSS_FEEDS = [
    {
        "url": "https://www.medscape.com/rss/public/diabetes.xml",
        "source_override": "Medscape Diabetes", # 显示在新闻卡片上的来源名称
        "target_categories": ["最新研究", "治疗进展"] # 这条RSS源的新闻会被尝试放入这些分类
    },
    {
        "url": "https://www.healio.com/news/endocrinology/rss",
        "source_override": "Healio Endocrinology",
        "target_categories": ["最新研究", "治疗进展"]
    },
    # --- 请您替换或添加以下占位符 ---
    # {
    #     "url": "YOUR_PUBMED_DIABETES_RESEARCH_RSS_URL", # 例如: https://pubmed.ncbi.nlm.nih.gov/rss/search/...
    #     "source_override": "PubMed",
    #     "target_categories": ["最新研究"]
    # },
    # {
    #     "url": "YOUR_ADA_DIABETES_CARE_RSS_URL", # 例如: https://diabetesjournals.org/care/rss/current.xml
    #     "source_override": "Diabetes Care (ADA)",
    #     "target_categories": ["最新研究", "治疗进展"]
    # },
]

# --- (2) 配置网站展示的分类及对应的 Google News 补充关键词 ---
CATEGORIES_CONFIG = {
    "最新研究": {
        "keywords": "糖尿病 最新论文 OR 糖尿病技术突破 OR 糖尿病机制研究 OR 医学会议糖尿病 OR GLP-1糖尿病 OR SGLT2糖尿病 OR 胰岛β细胞 OR 胰岛素敏感性",
        "emoji": "🔬"
    },
    "治疗进展": {
        "keywords": "糖尿病新药 OR 糖尿病适应症扩展 OR 糖尿病设备研发 OR AI辅助诊疗糖尿病 OR 达格列净 OR 司美格鲁肽 OR CGM OR 连续血糖监测 OR 胰岛素泵",
        "emoji": "💊"
    },
    "饮食与营养": {
        "keywords": "糖尿病饮食指南 OR 碳水交换表 OR 糖尿病食谱 OR 低GI饮食糖尿病 OR 高蛋白饮食糖尿病 OR 间歇性断食糖尿病 OR 膳食纤维糖尿病",
        "emoji": "🥗"
    },
    "预防与生活方式": {
        "keywords": "糖尿病运动建议 OR 糖尿病睡眠 OR 糖尿病减重 OR 控糖计划 OR 糖尿病早期筛查 OR 糖耐量异常 OR 体脂管理糖尿病 OR 糖尿病步数目标",
        "emoji": "🏃‍♀️"
    },
    "并发症管理": {
        "keywords": "糖尿病足 OR 糖尿病视网膜病变 OR 糖尿病肾病 OR 糖尿病神经病变 OR 糖网病 OR 微血管病变糖尿病 OR 尿白蛋白糖尿病",
        "emoji": "🩺"
    },
    "患者故事与心理支持": {
        "keywords": "糖尿病控糖经验 OR 糖尿病心理支持 OR 糖尿病家庭支持 OR 糖尿病患者故事 OR 糖尿病医生问答",
        "emoji": "😊"
    },
    "政策/医保信息": {
        "keywords": "糖尿病药品纳保 OR 糖尿病医保报销 OR 糖尿病社区慢病随访 OR 国家药监局糖尿病政策 OR 医保局糖尿病政策",
        "emoji": "📄"
    }
}

# --- 帮助函数：判断日期是否在最近一个月内 ---
def is_within_last_month_rss(time_struct, today_date_obj):
    if not time_struct:
        return False
    try:
        article_date = datetime.date(time_struct.tm_year, time_struct.tm_mon, time_struct.tm_mday)
        thirty_days_ago = today_date_obj - datetime.timedelta(days=30)
        return thirty_days_ago <= article_date <= today_date_obj
    except Exception as e:
        print(f"      [is_within_last_month_rss] 日期转换错误: {e} - Time Struct: {time_struct}")
        return False

# --- 帮助函数：清理 HTML ---
def clean_html(raw_html):
    if not raw_html: return ""
    try:
        return BeautifulSoup(raw_html, "html.parser").get_text()
    except Exception: return raw_html

# --- (3) 通用的从单个 RSS URL 获取文章的函数 ---
def fetch_articles_from_rss(rss_url, source_name_override=None):
    """
    从给定的 RSS URL 获取文章列表。
    返回解析后的文章列表，每篇文章是一个包含 title, link, snippet, source, time_struct 的字典。
    """
    print(f"    正在从源获取: {rss_url} ({source_name_override or '未知源'})")
    articles = []
    try:
        headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'}
        response = requests.get(rss_url, headers=headers, timeout=20)
        response.raise_for_status()
        feed = feedparser.parse(response.content)

        if not feed.entries:
            print(f"      此源未返回任何条目: {rss_url}")
            return []

        print(f"      从此源原始获取到 {len(feed.entries)} 条新闻。")
        for entry in feed.entries:
            title = entry.get("title", "无标题")
            link = entry.get("link", f"javascript:void(0);_{title}") # 添加标题确保链接唯一性（如果link为空）
            
            published_time_struct = entry.get("published_parsed") or entry.get("updated_parsed")
            
            summary_html = entry.get("summary", entry.get("description", "暂无摘要")) # 有些源用description
            snippet = clean_html(summary_html)
            
            # 来源处理：优先使用 override，其次是 feedparser 提供的，最后是 Google News (如果是 Google News 源)
            actual_source_name = source_name_override
            if not actual_source_name:
                source_info = entry.get("source")
                actual_source_name = source_info.get("title") if source_info else "未知来源"
                if "news.google.com" in rss_url and not source_name_override : # 特殊处理Google News的来源
                    # Google News RSS的条目title通常是 "新闻标题 - 来源"，尝试提取
                    if ' - ' in title:
                        parts = title.rsplit(' - ', 1)
                        # title = parts[0] # 取消修改标题，让用户看到完整信息
                        actual_source_name = parts[1]


            articles.append({
                "title": title,
                "url": link,
                "snippet": snippet,
                "source": actual_source_name,
                "time_struct": published_time_struct # 保留 time_struct 用于后续日期过滤
            })
    except requests.exceptions.Timeout:
        print(f"      获取源时发生超时错误: {rss_url}")
    except requests.exceptions.RequestException as e:
        print(f"      获取源时发生网络错误: {e} (URL: {rss_url})")
    except Exception as e:
        print(f"      处理源时发生未知错误: {e} (URL: {rss_url})")
    return articles

# --- HTML 生成逻辑 (与之前版本基本一致) ---
def generate_html_content(all_news_data):
    current_time_str = datetime.datetime.now().strftime('%Y年%m月%d日 %H:%M:%S')
    app_timezone = os.getenv('APP_TIMEZONE', 'UTC') 
    if app_timezone != 'UTC':
        try:
             current_time_str = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).strftime('%Y年%m月%d日 %H:%M:%S %Z')
        except Exception as e:
            print(f"应用时区 ({app_timezone}) 时出错: {e}。将使用默认服务器时间。")
            current_time_str = datetime.datetime.now().strftime('%Y年%m月%d日 %H:%M:%S (服务器时间)')

    current_year = datetime.datetime.now().year
    github_repo_url = f"https://github.com/{os.getenv('GITHUB_REPOSITORY', 'doudou-ux/diabetes-news')}"

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
            <p class="text-gray-600 mt-3 text-base md:text-lg">最近一个月动态（自动更新于：<span id="updateTime">{current_time_str}</span>）</p>
            <p class="text-sm text-gray-500 mt-2">资讯综合来源</p>
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
        document.addEventListener('DOMContentLoaded', function() {{
            const newsContent = document.getElementById('news-content');
            const loadingIndicator = document.getElementById('loading-indicator');
            let hasRealContent = false;
            if (newsContent && newsContent.childNodes) {{
                for (let i = 0; i < newsContent.childNodes.length; i++) {{
                    if (newsContent.childNodes[i].id !== 'loading-indicator' && newsContent.childNodes[i].nodeType === 1) {{
                        hasRealContent = true;
                        break;
                    }}
                }}
            }}
            if (loadingIndicator && hasRealContent) {{
                // loadingIndicator.style.display = 'none'; 
            }} else if (loadingIndicator && newsContent && newsContent.textContent && newsContent.textContent.includes("未能加载")) {{
                loadingIndicator.style.display = 'none';
            }}
        }});
    </script>
</body>
</html>"""

    news_html_parts = []
    found_any_news = False

    if not all_news_data or all(not articles for articles in all_news_data.values()):
        news_html_parts.append('<p class="text-center text-gray-500 text-xl py-10">抱歉，目前未能加载到最近一个月相关的糖尿病资讯。</p>')
    else:
        for category, articles in all_news_data.items():
            category_emoji = CATEGORIES_CONFIG.get(category, {}).get("emoji", "")
            category_html = f"""
            <section class="category-section">
                <h2 class="font-semibold category-title-text">{category_emoji} {html.escape(category)}</h2>
            """
            if not articles:
                category_html += '<p class="text-gray-500">最近一个月暂无该分类下的资讯。</p>'
            else:
                found_any_news = True
                category_html += '<div class="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-6 md:gap-8 news-grid">'
                for article in articles: 
                    title = html.escape(article.get('title', '无标题'))
                    url = html.escape(article.get('url', 'javascript:void(0);'))
                    snippet_raw = article.get('snippet', '暂无摘要')
                    snippet = html.escape(snippet_raw[:150] + ('...' if len(snippet_raw) > 150 else ''))
                    source_display = html.escape(article.get('source', '未知来源')) # 使用 article 中已处理好的 source
                    time_display = html.escape(article.get('time_display_str', '未知时间')) # 使用 article 中已处理好的 time_display_str

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
                            <span class="source-tag">{source_display}</span>
                            <span class="time-tag">{time_display}</span>
                        </div>
                    </div>
                    """
                category_html += "</div>" 
            category_html += "</section>"
            news_html_parts.append(category_html)
    
    loading_indicator_html = """<div id="loading-indicator" class="text-center py-10">
                 <div class="loader"></div>
                 <p class="text-gray-600 mt-2">正在加载最新资讯...</p>
            </div>"""
    if loading_indicator_html in html_output:
        if found_any_news or (news_html_parts and "<p class=\"text-center text-gray-500 text-xl py-10\">" in news_html_parts[0]):
            html_output = html_output.replace(loading_indicator_html, "".join(news_html_parts))
        elif not news_html_parts: 
             html_output = html_output.replace(loading_indicator_html, '<p class="text-center text-gray-500 text-xl py-10">资讯加载时出现问题或暂无内容。</p>')
    else: 
        html_output = html_output.replace(
            '<div id="news-content" class="space-y-12">\n            \n            ', 
            f'<div id="news-content" class="space-y-12">\n{"".join(news_html_parts)}'
        )
    return html_output

# --- (4) 主执行逻辑 ---
if __name__ == "__main__":
    print("开始从多个 RSS 源生成糖尿病资讯网页...")
    
    # 初始化存储所有分类文章的字典
    # key 是网站分类名, value 是该分类下的文章列表
    all_articles_by_site_category = {category_name: [] for category_name in CATEGORIES_CONFIG.keys()}
    seen_article_links = set() # 用于去重
    today = datetime.date.today()
    MAX_ARTICLES_PER_CATEGORY = 10

    # --- 步骤一：从权威 RSS 源获取新闻 ---
    print("\n--- 正在从权威 RSS 源获取新闻 ---")
    for feed_info in AUTHORITATIVE_RSS_FEEDS:
        raw_articles_from_feed = fetch_articles_from_rss(feed_info["url"], feed_info["source_override"])
        print(f"  处理来自 {feed_info['source_override']} 的 {len(raw_articles_from_feed)} 条原始新闻...")
        
        for article_data in raw_articles_from_feed:
            if article_data["url"] in seen_article_links:
                print(f"    跳过重复文章 (来自权威源): {article_data['title'][:30]}...")
                continue

            if is_within_last_month_rss(article_data["time_struct"], today):
                time_display_str = "未知时间"
                if article_data["time_struct"]:
                    try:
                        time_display_str = time.strftime("%Y-%m-%d", article_data["time_struct"])
                    except: pass
                
                processed_article = {
                    "title": article_data["title"],
                    "url": article_data["url"],
                    "snippet": article_data["snippet"],
                    "source": article_data["source"], # 使用 fetch_articles_from_rss 中处理好的 source
                    "time_display_str": time_display_str
                }

                seen_article_links.add(article_data["url"])
                
                # 将文章放入其对应的所有目标分类
                for target_cat in feed_info["target_categories"]:
                    if target_cat in all_articles_by_site_category: # 确保目标分类是我们网站定义的分类
                        if len(all_articles_by_site_category[target_cat]) < MAX_ARTICLES_PER_CATEGORY:
                            all_articles_by_site_category[target_cat].append(processed_article)
                            print(f"      已添加 '{processed_article['title'][:30]}...' 到分类 '{target_cat}' (来自 {feed_info['source_override']})")
                        else:
                            print(f"      分类 '{target_cat}' 已满 (10条)，无法添加来自 {feed_info['source_override']} 的更多文章。")
                    else:
                        print(f"      警告: 权威源指定的目标分类 '{target_cat}' 不在网站分类配置中，跳过。")
        time.sleep(1) # 友好访问

    # --- 步骤二：从 Google News RSS 获取补充新闻 ---
    print("\n--- 正在从 Google News RSS 获取补充新闻 ---")
    for site_category_name, config in CATEGORIES_CONFIG.items():
        if len(all_articles_by_site_category[site_category_name]) < MAX_ARTICLES_PER_CATEGORY:
            print(f"  分类 '{site_category_name}' 需要补充新闻 (当前 {len(all_articles_by_site_category[site_category_name])} 条)，将从 Google News 获取...")
            
            google_news_rss_url = f"https://news.google.com/rss/search?q={html.escape(config['keywords'])}&hl=zh-CN&gl=CN&ceid=CN:zh-Hans"
            # 注意：Google News 的 source_override 留空，让 fetch_articles_from_rss 尝试从标题提取
            raw_articles_from_google = fetch_articles_from_rss(google_news_rss_url, source_name_override=None) 
            
            print(f"    处理来自 Google News (分类: {site_category_name}) 的 {len(raw_articles_from_google)} 条原始新闻...")
            for article_data in raw_articles_from_google:
                if len(all_articles_by_site_category[site_category_name]) >= MAX_ARTICLES_PER_CATEGORY:
                    break # 这个分类已经满了

                if article_data["url"] in seen_article_links:
                    print(f"      跳过重复文章 (来自Google News): {article_data['title'][:30]}...")
                    continue

                if is_within_last_month_rss(article_data["time_struct"], today):
                    time_display_str = "未知时间"
                    if article_data["time_struct"]:
                        try:
                            time_display_str = time.strftime("%Y-%m-%d", article_data["time_struct"])
                        except: pass

                    processed_article = {
                        "title": article_data["title"],
                        "url": article_data["url"],
                        "snippet": article_data["snippet"],
                        "source": article_data["source"], # 使用 fetch_articles_from_rss 中处理好的 source
                        "time_display_str": time_display_str
                    }
                    all_articles_by_site_category[site_category_name].append(processed_article)
                    seen_article_links.add(article_data["url"])
                    print(f"        已添加 '{processed_article['title'][:30]}...' 到分类 '{site_category_name}' (来自 Google News)")
            time.sleep(1) # 友好访问
        else:
            print(f"  分类 '{site_category_name}' 已有足够新闻 ({len(all_articles_by_site_category[site_category_name])}条)，跳过 Google News 获取。")

    # --- (5) 生成最终的HTML ---
    final_html = generate_html_content(all_articles_by_site_category)
    
    output_filename = "index.html" 
    try:
        with open(output_filename, "w", encoding="utf-8") as f:
            f.write(final_html)
        print(f"\n成功生成网页：{output_filename}")
    except IOError as e:
        print(f"\n错误：无法写入文件 {output_filename}。错误信息: {e}")
    except Exception as e:
        print(f"\n生成过程中发生未知错误: {e}")

    print("资讯网页生成完毕。")

