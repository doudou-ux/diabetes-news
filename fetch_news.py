# fetch_news.py
import datetime
import html
import os
import time
import requests
import feedparser
from bs4 import BeautifulSoup

# --- (1) 配置权威 RSS 源 ---
AUTHORITATIVE_RSS_FEEDS = [
    {
        "url": "https://www.medscape.com/rss/public/diabetes.xml",
        "source_override": "Medscape Diabetes",
        "target_categories": ["最新研究", "治疗进展"]
    },
    {
        "url": "https://www.healio.com/news/endocrinology/rss",
        "source_override": "Healio Endocrinology",
        "target_categories": ["最新研究", "治疗进展"]
    },
    # {
    #     "url": "YOUR_PUBMED_DIABETES_RESEARCH_RSS_URL",
    #     "source_override": "PubMed",
    #     "target_categories": ["最新研究"]
    # },
    # {
    #     "url": "YOUR_ADA_DIABETES_CARE_RSS_URL",
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
            link = entry.get("link", f"javascript:void(0);_{html.escape(title)}")
            published_time_struct = entry.get("published_parsed") or entry.get("updated_parsed")
            summary_html = entry.get("summary", entry.get("description", "暂无摘要"))
            snippet = clean_html(summary_html)
            actual_source_name = source_name_override
            if not actual_source_name:
                source_info = entry.get("source")
                actual_source_name = source_info.get("title") if source_info else "未知来源"
                if "news.google.com" in rss_url and not source_name_override :
                    if ' - ' in title:
                        parts = title.rsplit(' - ', 1)
                        actual_source_name = parts[1]
            articles.append({
                "title": title,
                "url": link,
                "snippet": snippet,
                "source": actual_source_name,
                "time_struct": published_time_struct
            })
    except requests.exceptions.Timeout:
        print(f"      获取源时发生超时错误: {rss_url}")
    except requests.exceptions.RequestException as e:
        print(f"      获取源时发生网络错误: {e} (URL: {rss_url})")
    except Exception as e:
        print(f"      处理源时发生未知错误: {e} (URL: {rss_url})")
    return articles

# --- (4) HTML 生成逻辑 ---
def generate_html_content(all_news_data_sorted): # 接收已排序的数据
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

    # --- HTML 头部和通用样式 ---
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
        /* .category-section (现在是 .tab-content) 的样式可以保持或调整 */
        .tab-content {{ background-color: #ffffff; border-radius: 0.75rem; box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06); padding: 1.25rem; display: none; }}
        .tab-content.active {{ display: block; }} /* Active tab content is visible */
        @media (min-width: 768px) {{ .tab-content {{ padding: 2rem; }} }}
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
        /* Tab button styles */
        .tab-buttons-container {{ 
            display: flex; 
            flex-wrap: wrap; /* Allow tabs to wrap on smaller screens */
            margin-bottom: 1.5rem; 
            border-bottom: 2px solid #d1d5db; /* gray-300 */
        }}
        .tab-button {{
            padding: 0.75rem 1.25rem; /* py-3 px-5 */
            margin-right: 0.5rem; /* mr-2 */
            margin-bottom: -2px; /* Overlap with bottom border */
            border: 2px solid transparent;
            border-bottom: none;
            border-top-left-radius: 0.375rem; /* rounded-t-md */
            border-top-right-radius: 0.375rem; /* rounded-t-md */
            cursor: pointer;
            font-weight: 500; /* font-medium */
            color: #4b5563; /* gray-600 */
            background-color: #f9fafb; /* gray-50 */
            transition: all 0.3s ease;
        }}
        .tab-button:hover {{
            color: #1f2937; /* gray-800 */
            background-color: #f3f4f6; /* gray-100 */
        }}
        .tab-button.active {{
            color: #1d4ed8; /* blue-700 */
            background-color: #ffffff; /* white */
            border-color: #d1d5db; /* gray-300 */
            border-bottom-color: #ffffff; /* white, to make it look connected to content */
        }}
    </style>
</head>
<body class="bg-gray-100 text-gray-800">
    <div class="container mx-auto main-container">
        <header class="text-center mb-10 md:mb-16">
            <h1 class="font-bold text-blue-700 header-main-title">糖尿病前沿资讯</h1>
            <p class="text-gray-600 mt-3 text-base md:text-lg">最近一个月动态（自动更新于：<span id="updateTime">{current_time_str}</span>）</p>
            <p class="text-sm text-gray-500 mt-2">资讯综合来源</p>
        </header>

        <div class="tab-buttons-container" id="tabButtons">
"""
    # --- 生成标签页按钮 ---
    first_category = True
    for category_name_key in all_news_data_sorted.keys(): # 使用排序后字典的键顺序
        category_config = CATEGORIES_CONFIG.get(category_name_key, {})
        emoji = category_config.get("emoji", "")
        # 创建一个对JavaScript友好的ID (例如，替换空格和特殊字符)
        tab_id = "tab-" + html.escape(category_name_key.replace(" ", "-").replace("/", "-").lower())
        active_class = "active" if first_category else ""
        html_output += f"""
            <button class="tab-button {active_class}" data-tab-target="#{tab_id}">{emoji} {html.escape(category_name_key)}</button>
        """
        first_category = False
    html_output += """
        </div>

        <div id="news-content">
"""
    # --- 生成标签页内容 ---
    first_category = True # Reset for content panes
    if not any(all_news_data_sorted.values()): # Check if all categories are empty after sorting
        html_output += '<p class="text-center text-gray-500 text-xl py-10">抱歉，目前未能加载到最近一个月相关的糖尿病资讯。</p>'
    else:
        for category_name_key, articles in all_news_data_sorted.items():
            category_config = CATEGORIES_CONFIG.get(category_name_key, {})
            emoji = category_config.get("emoji", "")
            tab_id = "tab-" + html.escape(category_name_key.replace(" ", "-").replace("/", "-").lower())
            active_class = "active" if first_category else ""
            
            category_html_content = f"""
            <div id="{tab_id}" class="tab-content {active_class}">
                <h2 class="font-semibold category-title-text">{emoji} {html.escape(category_name_key)}</h2>
            """
            if not articles:
                category_html_content += '<p class="text-gray-500">最近一个月暂无该分类下的资讯。</p>'
            else:
                category_html_content += '<div class="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-6 md:gap-8 news-grid">'
                for article in articles: 
                    title = html.escape(article.get('title', '无标题'))
                    url = html.escape(article.get('url', 'javascript:void(0);'))
                    snippet_raw = article.get('snippet', '暂无摘要')
                    snippet = html.escape(snippet_raw[:150] + ('...' if len(snippet_raw) > 150 else ''))
                    source_display = html.escape(article.get('source', '未知来源'))
                    time_display = html.escape(article.get('time_display_str', '未知时间'))

                    category_html_content += f"""
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
                category_html_content += "</div>" 
            category_html_content += "</section></div>" # Close tab-content and section (section might be redundant here)
            html_output += category_html_content
            first_category = False

    # --- HTML 尾部和 JavaScript ---
    html_output += f"""
        </div> </div> <footer class="text-center p-6 mt-12 text-gray-600 text-sm border-t border-gray-300">
        <p>&copy; {current_year} 糖尿病资讯聚合. <a href="{github_repo_url}" target="_blank" rel="noopener noreferrer" class="text-blue-600 hover:underline">项目源码</a></p>
        <p class="mt-1">本站内容仅供参考, 不构成医疗建议。</p>
    </footer>
    <script>
        document.addEventListener('DOMContentLoaded', function () {{
            const tabButtons = document.querySelectorAll('.tab-button');
            const tabContents = document.querySelectorAll('.tab-content');

            if (tabButtons.length > 0 && tabContents.length > 0) {{
                // Activate the first tab by default if not already active via class
                // (The Python script already adds 'active' to the first tab button and content)
                // if (!document.querySelector('.tab-button.active')) {{
                //     tabButtons[0].classList.add('active');
                // }}
                // if (!document.querySelector('.tab-content.active')) {{
                //    tabContents[0].classList.add('active');
                // }}
                
                tabButtons.forEach(button => {{
                    button.addEventListener('click', () => {{
                        // Deactivate all tabs and contents
                        tabButtons.forEach(btn => btn.classList.remove('active'));
                        tabContents.forEach(content => content.classList.remove('active'));

                        // Activate clicked tab and corresponding content
                        button.classList.add('active');
                        const targetContentId = button.dataset.tabTarget; // e.g., "#tab-最新研究"
                        const targetContent = document.querySelector(targetContentId);
                        if (targetContent) {{
                            targetContent.classList.add('active');
                        }}
                    }});
                }});
            }}
        }});
    </script>
</body>
</html>"""
    return html_output

# --- (5) 主执行逻辑 ---
if __name__ == "__main__":
    print("开始从多个 RSS 源生成糖尿病资讯网页...")
    
    all_articles_by_site_category_temp = {category_name: [] for category_name in CATEGORIES_CONFIG.keys()}
    seen_article_links = set()
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
                    "source": article_data["source"],
                    "time_display_str": time_display_str,
                    "time_struct": article_data["time_struct"] # 保留用于排序
                }
                seen_article_links.add(article_data["url"])
                for target_cat in feed_info["target_categories"]:
                    if target_cat in all_articles_by_site_category_temp:
                        # 先收集，不立即限制数量，排序后再限制
                        all_articles_by_site_category_temp[target_cat].append(processed_article)
                        print(f"      已收集 '{processed_article['title'][:30]}...' 到分类 '{target_cat}' (来自 {feed_info['source_override']})")
        time.sleep(1)

    # --- 步骤二：从 Google News RSS 获取补充新闻 ---
    print("\n--- 正在从 Google News RSS 获取补充新闻 ---")
    for site_category_name, config in CATEGORIES_CONFIG.items():
        # Google News 作为补充，如果权威源已经有很多，可以少取或不取
        # 这里简单处理：总是尝试获取，排序和去重会在后面统一处理
        print(f"  将为分类 '{site_category_name}' 从 Google News 获取补充新闻...")
        google_news_rss_url = f"https://news.google.com/rss/search?q={html.escape(config['keywords'])}&hl=zh-CN&gl=CN&ceid=CN:zh-Hans"
        raw_articles_from_google = fetch_articles_from_rss(google_news_rss_url, source_name_override=None)
        
        print(f"    处理来自 Google News (分类: {site_category_name}) 的 {len(raw_articles_from_google)} 条原始新闻...")
        for article_data in raw_articles_from_google:
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
                    "source": article_data["source"],
                    "time_display_str": time_display_str,
                    "time_struct": article_data["time_struct"] # 保留用于排序
                }
                all_articles_by_site_category_temp[site_category_name].append(processed_article)
                seen_article_links.add(article_data["url"])
                print(f"        已收集 '{processed_article['title'][:30]}...' 到分类 '{site_category_name}' (来自 Google News)")
        time.sleep(1)

    # --- 步骤三：对每个分类的文章按日期排序并截取前 MAX_ARTICLES_PER_CATEGORY 条 ---
    print("\n--- 正在对各分类新闻进行排序和截取 ---")
    all_articles_by_site_category_final_sorted = {}
    for category_name, articles_list in all_articles_by_site_category_temp.items():
        # 排序：最新的在前 (time_struct 越大越新)
        # 处理 time_struct 可能为 None 的情况，将其排在最后
        articles_list.sort(key=lambda x: time.mktime(x["time_struct"]) if x["time_struct"] else -float('inf'), reverse=True)
        all_articles_by_site_category_final_sorted[category_name] = articles_list[:MAX_ARTICLES_PER_CATEGORY]
        print(f"  分类 '{category_name}' 排序并截取后有 {len(all_articles_by_site_category_final_sorted[category_name])} 条新闻。")

    # --- (6) 生成最终的HTML ---
    final_html = generate_html_content(all_articles_by_site_category_final_sorted)
    
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
