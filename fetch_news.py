# fetch_news.py
import datetime
import html
import os
import re # 用于规范化标题
import time
import requests
import feedparser
from bs4 import BeautifulSoup
from urllib.parse import urljoin # 用于爬虫构建完整URL
from deep_translator import GoogleTranslator # 用于翻译

# --- (1) 配置权威 RSS 源 ---
# 添加 priority 字段，数值越大越优先
AUTHORITATIVE_RSS_FEEDS = [
    {
        "url": "https://www.medscape.com/rss/public/diabetes.xml",
        "source_override": "Medscape Diabetes",
        "target_categories": ["最新研究", "治疗进展"],
        "priority": 10
    },
    {
        "url": "https://www.healio.com/news/endocrinology/rss",
        "source_override": "Healio Endocrinology",
        "target_categories": ["最新研究", "治疗进展"],
        "priority": 9
    },
    {
        "url": "https://www.diabettech.com/feed/",
        "source_override": "Diabettech",
        "target_categories": ["治疗进展", "最新研究"], # 糖尿病技术相关
        "priority": 8
    },
    {
        "url": "https://thesavvydiabetic.com/feed/",
        "source_override": "The Savvy Diabetic",
        "target_categories": ["患者故事与心理支持", "预防与生活方式"], # 博客，经验分享
        "priority": 7
    },
    {
        "url": "https://forum.diabetes.org.uk/boards/forums/-/index.rss",
        "source_override": "Diabetes UK 论坛",
        "target_categories": ["患者故事与心理支持"], # 论坛讨论
        "priority": 6
    },
    # --- 请您替换或添加以下占位符，并指定 priority ---
    # {
    #     "url": "YOUR_PUBMED_DIABETES_RESEARCH_RSS_URL",
    #     "source_override": "PubMed",
    #     "target_categories": ["最新研究"],
    #     "priority": 12
    # },
    # {
    #     "url": "YOUR_ADA_DIABETES_CARE_RSS_URL",
    #     "source_override": "Diabetes Care (ADA)",
    #     "target_categories": ["最新研究", "治疗进展"],
    #     "priority": 11
    # },
]

# --- (1b) 配置爬虫源 ---
# 注意：这些爬虫目前不提取发布日期，可能无法通过日期过滤
SCRAPED_SOURCES_CONFIG = [
    {
        "name": "Breakthrough T1D News",
        "fetch_function": "fetch_breakthrought1d_articles", # 对应的爬虫函数名
        "source_override": "Breakthrough T1D",
        "target_categories": ["最新研究", "治疗进展"], # 主要关于1型糖尿病研究
        "priority": 8 # 权威机构，但爬虫稳定性可能稍差
    },
    {
        "name": "MyGlu Articles",
        "fetch_function": "fetch_myglu_articles",
        "source_override": "MyGlu",
        "target_categories": ["患者故事与心理支持", "预防与生活方式"],
        "priority": 7
    },
    {
        "name": "DZD News (2025)", # 假设每年更新URL或爬虫逻辑
        "fetch_function": "fetch_dzd_articles",
        "source_override": "DZD News",
        "target_categories": ["最新研究"],
        "priority": 9 # 德国糖尿病研究中心
    },
]

GOOGLE_NEWS_PRIORITY = 1 # Google News 作为补充来源，优先级较低

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

# --- 帮助函数：规范化标题 ---
def normalize_title(title):
    if not title: return ""
    title = title.lower()
    title = re.sub(r'[^\w\s-]', '', title)
    title = re.sub(r'\s+', ' ', title).strip()
    return title

# --- 帮助函数：判断日期是否在最近一个月内 ---
def is_within_last_month_rss(time_struct, today_date_obj):
    if not time_struct: return False # 如果没有时间信息，则认为不在范围内
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
    try: return BeautifulSoup(raw_html, "html.parser").get_text(separator=' ', strip=True)
    except Exception: return raw_html

# --- (A) RSS 源获取函数 ---
def fetch_articles_from_rss(rss_url, source_name_override=None):
    # (此函数内容与 diabetes_news_fetch_dedup_priority 中的 fetch_articles_from_rss 基本相同)
    # ... (为简洁起见，假设此函数已正确定义，与上一版本相同)
    print(f"    正在从 RSS 源获取: {rss_url} ({source_name_override or '未知源'})")
    articles = []
    try:
        headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'}
        response = requests.get(rss_url, headers=headers, timeout=20)
        response.raise_for_status()
        feed = feedparser.parse(response.content)
        if not feed.entries:
            print(f"      此 RSS 源未返回任何条目: {rss_url}")
            return []
        print(f"      从此 RSS 源原始获取到 {len(feed.entries)} 条新闻。")
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
                if "news.google.com" in rss_url and not source_name_override:
                    if ' - ' in title:
                        parts = title.rsplit(' - ', 1)
                        actual_source_name = parts[1]
            articles.append({
                "title": title, "url": link, "snippet": snippet,
                "source": actual_source_name, "time_struct": published_time_struct
            })
    except requests.exceptions.Timeout: print(f"      获取 RSS 源时发生超时错误: {rss_url}")
    except requests.exceptions.RequestException as e: print(f"      获取 RSS 源时发生网络错误: {e} (URL: {rss_url})")
    except Exception as e: print(f"      处理 RSS 源时发生未知错误: {e} (URL: {rss_url})")
    return articles


# --- (B) 爬虫函数定义 ---
def translate_text(text, target_lang='zh-CN'): # zh-CN for Simplified Chinese
    if not text: return ""
    try:
        # GoogleTranslator(source='auto', target='zh') is deprecated, use 'zh-CN' or 'zh-TW'
        return GoogleTranslator(source='auto', target=target_lang).translate(text)
    except Exception as e:
        print(f"      翻译错误: {e} - 原文: {text[:50]}...")
        return text # 翻译失败则返回原文

def fetch_breakthrought1d_articles():
    BASE_URL = "https://www.breakthrought1d.org/news/"
    print(f"    正在爬取: {BASE_URL}")
    articles = []
    try:
        headers = {"User-Agent": "Mozilla/5.0"}
        response = requests.get(BASE_URL, headers=headers, timeout=20)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, "html.parser")

        for article_el in soup.select("article"): # Element selector might need adjustment
            a_tag = article_el.select_one("h2 a")
            if not a_tag: continue

            title_en = a_tag.get_text(strip=True)
            link = urljoin(BASE_URL, a_tag.get("href"))
            summary_tag = article_el.select_one("p") # Selector for summary
            summary_en = summary_tag.get_text(strip=True) if summary_tag else ""
            
            # 尝试提取日期 - 这部分需要根据实际网页结构添加
            date_str = None # Placeholder
            time_struct = None # Placeholder
            # Example: date_tag = article_el.select_one(".date-class")
            # if date_tag: date_str = date_tag.get_text(strip=True)
            # if date_str: try: time_struct = feedparser._parse_date(date_str) except: pass
            if not time_struct: print(f"      警告: 未能从 {link} 提取发布日期。")


            title_cn = translate_text(title_en)
            summary_cn = translate_text(summary_en)

            articles.append({
                "title": title_cn, "url": link, "snippet": summary_cn,
                "source": "Breakthrough T1D", # 在 SCRAPED_SOURCES_CONFIG 中已定义
                "time_struct": time_struct # 目前为 None
            })
    except Exception as e:
        print(f"      爬取 Breakthrough T1D 时出错: {e}")
    return articles

def fetch_myglu_articles():
    BASE_URL = "https://myglu.org/articles"
    print(f"    正在爬取: {BASE_URL}")
    articles = []
    try:
        headers = {"User-Agent": "Mozilla/5.0"}
        response = requests.get(BASE_URL, headers=headers, timeout=20)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, "html.parser")

        for item in soup.select("div.node-article"): # Element selector
            a_tag = item.select_one("h2 a")
            if not a_tag: continue

            title_en = a_tag.get_text(strip=True)
            link = urljoin(BASE_URL, a_tag.get("href"))
            summary_tag = item.select_one("div.field-name-body div.field-item") # More specific selector
            summary_en = summary_tag.get_text(strip=True) if summary_tag else ""

            date_str = None
            time_struct = None
            # Example: date_tag = item.select_one("span.date-display-single")
            # if date_tag: date_str = date_tag.get_text(strip=True)
            # if date_str: try: time_struct = feedparser._parse_date(date_str) except: pass
            if not time_struct: print(f"      警告: 未能从 {link} 提取发布日期。")


            title_cn = translate_text(title_en)
            summary_cn = translate_text(summary_en)
            
            articles.append({
                "title": title_cn, "url": link, "snippet": summary_cn,
                "source": "MyGlu", # 在 SCRAPED_SOURCES_CONFIG 中已定义
                "time_struct": time_struct # 目前为 None
            })
    except Exception as e:
        print(f"      爬取 MyGlu 时出错: {e}")
    return articles

def fetch_dzd_articles():
    BASE_URL = "https://www.dzd-ev.de/en/press/press-releases/press-releases-2025/index.html" # URL may need annual update
    print(f"    正在爬取: {BASE_URL}")
    articles = []
    try:
        headers = {"User-Agent": "Mozilla/5.0"}
        response = requests.get(BASE_URL, headers=headers, timeout=20)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, "html.parser")

        for item in soup.select("div.teaser-text"): # Element selector
            title_tag = item.select_one("a")
            if not title_tag: continue

            title_en = title_tag.get_text(strip=True)
            link = urljoin(BASE_URL, title_tag.get("href"))
            summary_tag = item.select_one("p")
            summary_en = summary_tag.get_text(strip=True) if summary_tag else ""

            date_str = None
            time_struct = None
            # Example: date_container = item.find_previous_sibling("div", class_="news-list-date")
            # if date_container: date_str = date_container.get_text(strip=True)
            # if date_str: try: time_struct = feedparser._parse_date(date_str) except: pass
            if not time_struct: print(f"      警告: 未能从 {link} 提取发布日期。")

            title_cn = translate_text(title_en)
            summary_cn = translate_text(summary_en)

            articles.append({
                "title": title_cn, "url": link, "snippet": summary_cn,
                "source": "DZD News", # 在 SCRAPED_SOURCES_CONFIG 中已定义
                "time_struct": time_struct # 目前为 None
            })
    except Exception as e:
        print(f"      爬取 DZD News 时出错: {e}")
    return articles

# 将函数名字符串映射到实际函数对象
SCRAPER_FUNCTIONS_MAP = {
    "fetch_breakthrought1d_articles": fetch_breakthrought1d_articles,
    "fetch_myglu_articles": fetch_myglu_articles,
    "fetch_dzd_articles": fetch_dzd_articles,
}


# --- HTML 生成逻辑 (与之前版本 diabetes_news_fetch_tabs_v1 一致) ---
def generate_html_content(all_news_data_sorted):
    # (此函数内容与 diabetes_news_fetch_tabs_v1 中的 generate_html_content 完全相同，此处省略以减少重复)
    # 请确保您使用的是 diabetes_news_fetch_tabs_v1 中的完整 generate_html_content 函数
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
        .tab-content {{ background-color: #ffffff; border-radius: 0.75rem; box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06); padding: 1.25rem; display: none; }}
        .tab-content.active {{ display: block; }}
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
        .tab-buttons-container {{ display: flex; flex-wrap: wrap; margin-bottom: 1.5rem; border-bottom: 2px solid #d1d5db; }}
        .tab-button {{ padding: 0.75rem 1.25rem; margin-right: 0.5rem; margin-bottom: -2px; border: 2px solid transparent; border-bottom: none; border-top-left-radius: 0.375rem; border-top-right-radius: 0.375rem; cursor: pointer; font-weight: 500; color: #4b5563; background-color: #f9fafb; transition: all 0.3s ease; }}
        .tab-button:hover {{ color: #1f2937; background-color: #f3f4f6; }}
        .tab-button.active {{ color: #1d4ed8; background-color: #ffffff; border-color: #d1d5db; border-bottom-color: #ffffff; }}
    </style>
</head>
<body class="bg-gray-100 text-gray-800">
    <div class="container mx-auto main-container">
        <header class="text-center mb-10 md:mb-16">
            <h1 class="font-bold text-blue-700 header-main-title">糖尿病前沿资讯</h1>
            <p class="text-gray-600 mt-3 text-base md:text-lg">最近一个月动态（自动更新于：<span id="updateTime">{current_time_str}</span>）</p>
            <p class="text-sm text-gray-500 mt-2">资讯综合来源</p>
        </header>
        <div class="tab-buttons-container" id="tabButtons">"""
    first_category = True
    for category_name_key in all_news_data_sorted.keys():
        category_config = CATEGORIES_CONFIG.get(category_name_key, {})
        emoji = category_config.get("emoji", "")
        tab_id = "tab-" + html.escape(category_name_key.replace(" ", "-").replace("/", "-").lower())
        active_class = "active" if first_category else ""
        html_output += f"""<button class="tab-button {active_class}" data-tab-target="#{tab_id}">{emoji} {html.escape(category_name_key)}</button>"""
        first_category = False
    html_output += """</div><div id="news-content">"""
    first_category = True
    if not any(all_news_data_sorted.values()):
        html_output += '<p class="text-center text-gray-500 text-xl py-10">抱歉，目前未能加载到最近一个月相关的糖尿病资讯。</p>'
    else:
        for category_name_key, articles in all_news_data_sorted.items():
            category_config = CATEGORIES_CONFIG.get(category_name_key, {})
            emoji = category_config.get("emoji", "")
            tab_id = "tab-" + html.escape(category_name_key.replace(" ", "-").replace("/", "-").lower())
            active_class = "active" if first_category else ""
            category_html_content = f"""<div id="{tab_id}" class="tab-content {active_class}"><h2 class="font-semibold category-title-text">{emoji} {html.escape(category_name_key)}</h2>"""
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
                    category_html_content += f"""<div class="news-card shadow-md hover:shadow-lg p-5"><div class="card-content"><h3 class="text-lg mb-2"><a href="{url}" {'target="_blank" rel="noopener noreferrer"' if url != 'javascript:void(0);' else ''} class="news-item-link">{title}</a></h3><p class="text-gray-600 text-sm mb-4 leading-relaxed">{snippet}</p></div><div class="card-footer text-xs text-gray-500 flex flex-wrap gap-2 items-center pt-3 border-t border-gray-200"><span class="source-tag">{source_display}</span><span class="time-tag">{time_display}</span></div></div>"""
                category_html_content += "</div>" 
            category_html_content += "</div>"
            html_output += category_html_content
            first_category = False
    html_output += f"""</div> </div> <footer class="text-center p-6 mt-12 text-gray-600 text-sm border-t border-gray-300"><p>&copy; {current_year} 糖尿病资讯聚合. <a href="{github_repo_url}" target="_blank" rel="noopener noreferrer" class="text-blue-600 hover:underline">项目源码</a></p><p class="mt-1">本站内容仅供参考, 不构成医疗建议。</p></footer>
    <script>
        document.addEventListener('DOMContentLoaded', function () {{
            const tabButtons = document.querySelectorAll('.tab-button');
            const tabContents = document.querySelectorAll('.tab-content');
            if (tabButtons.length > 0 && tabContents.length > 0) {{
                tabButtons.forEach(button => {{
                    button.addEventListener('click', () => {{
                        tabButtons.forEach(btn => btn.classList.remove('active'));
                        tabContents.forEach(content => content.classList.remove('active'));
                        button.classList.add('active');
                        const targetContentId = button.dataset.tabTarget;
                        const targetContent = document.querySelector(targetContentId);
                        if (targetContent) {{ targetContent.classList.add('active'); }}
                    }});
                }});
            }}
        }});
    </script>
</body></html>"""
    return html_output

# --- (5) 主执行逻辑 ---
if __name__ == "__main__":
    print("开始从多个 RSS 源和爬虫生成糖尿病资讯网页...")
    
    unique_articles_candidates = {}
    globally_seen_urls = set()
    today = datetime.date.today()
    MAX_ARTICLES_PER_CATEGORY = 10

    # --- 步骤一：从权威 RSS 源获取新闻 ---
    print("\n--- 正在从权威 RSS 源获取新闻 ---")
    for feed_info in AUTHORITATIVE_RSS_FEEDS:
        current_priority = feed_info.get("priority", 5) 
        raw_articles_from_feed = fetch_articles_from_rss(feed_info["url"], feed_info["source_override"])
        for article_data in raw_articles_from_feed:
            if article_data["url"] in globally_seen_urls: continue
            if is_within_last_month_rss(article_data["time_struct"], today):
                normalized_title = normalize_title(article_data["title"])
                time_display_str = "未知时间"
                if article_data["time_struct"]:
                    try: time_display_str = time.strftime("%Y-%m-%d", article_data["time_struct"])
                    except: pass
                article_obj_for_storage = {
                    "title": article_data["title"], "url": article_data["url"],
                    "snippet": article_data["snippet"], "source": article_data["source"],
                    "time_display_str": time_display_str, "time_struct": article_data["time_struct"],
                    "source_priority": current_priority
                }
                if normalized_title not in unique_articles_candidates or \
                   current_priority > unique_articles_candidates[normalized_title]["priority"]:
                    unique_articles_candidates[normalized_title] = {
                        "article_obj": article_obj_for_storage, "priority": current_priority,
                        "target_categories": feed_info["target_categories"], "url": article_data["url"]
                    }
                    print(f"      候选(权威RSS): '{article_obj_for_storage['title'][:30]}...' (Prio: {current_priority})")
        time.sleep(1)

    # --- 步骤二：从爬虫源获取新闻 ---
    print("\n--- 正在从爬虫源获取新闻 ---")
    for scraper_info in SCRAPED_SOURCES_CONFIG:
        fetch_function = SCRAPER_FUNCTIONS_MAP.get(scraper_info["fetch_function"])
        if not fetch_function:
            print(f"    错误: 未找到爬虫函数 {scraper_info['fetch_function']}")
            continue
        
        current_priority = scraper_info.get("priority", 3) # 爬虫源默认优先级为3
        raw_articles_from_scraper = fetch_function() # 调用对应的爬虫函数
        
        for article_data in raw_articles_from_scraper: # article_data 结构应与RSS获取的类似
            if article_data["url"] in globally_seen_urls: continue
            # 注意：爬虫获取的文章目前 time_struct 为 None，所以 is_within_last_month_rss 会返回 False
            # 如果希望它们通过，需要修改 is_within_last_month_rss 或确保爬虫能提取日期
            if is_within_last_month_rss(article_data["time_struct"], today): # 这行可能导致爬虫内容被过滤
                normalized_title = normalize_title(article_data["title"])
                # time_display_str 对于爬虫内容，如果time_struct为None，则会是"未知时间"
                time_display_str = "未知时间"
                if article_data.get("time_struct"): # 确保有 time_struct
                    try: time_display_str = time.strftime("%Y-%m-%d", article_data["time_struct"])
                    except: pass
                else: # 如果爬虫没有提供日期，可以尝试使用一个固定文本或当天日期
                    # time_display_str = f"抓取于 {today.strftime('%Y-%m-%d')}" # 示例
                    print(f"      注意: 文章 '{article_data['title'][:30]}...' 无有效发布日期，将显示'未知时间'")


                article_obj_for_storage = {
                    "title": article_data["title"], "url": article_data["url"],
                    "snippet": article_data["snippet"], "source": scraper_info["source_override"],
                    "time_display_str": time_display_str, "time_struct": article_data["time_struct"],
                    "source_priority": current_priority
                }
                if normalized_title not in unique_articles_candidates or \
                   current_priority > unique_articles_candidates[normalized_title]["priority"]:
                    unique_articles_candidates[normalized_title] = {
                        "article_obj": article_obj_for_storage, "priority": current_priority,
                        "target_categories": scraper_info["target_categories"], "url": article_data["url"]
                    }
                    print(f"      候选(爬虫): '{article_obj_for_storage['title'][:30]}...' (Prio: {current_priority})")
            else:
                 print(f"      跳过(爬虫，日期过滤): '{article_data['title'][:30]}...' (Time struct: {article_data.get('time_struct')})")

        time.sleep(1)


    # --- 步骤三：从 Google News RSS 获取补充新闻 ---
    print("\n--- 正在从 Google News RSS 获取补充新闻 ---")
    for site_category_name, config in CATEGORIES_CONFIG.items():
        print(f"  为分类 '{site_category_name}' 从 Google News 获取补充...")
        google_news_rss_url = f"https://news.google.com/rss/search?q={html.escape(config['keywords'])}&hl=zh-CN&gl=CN&ceid=CN:zh-Hans"
        raw_articles_from_google = fetch_articles_from_rss(google_news_rss_url, source_name_override=None)
        
        for article_data in raw_articles_from_google:
            if article_data["url"] in globally_seen_urls and \
               any(article_data["url"] == cand["url"] for cand in unique_articles_candidates.values()):
                continue
            if is_within_last_month_rss(article_data["time_struct"], today):
                normalized_title = normalize_title(article_data["title"])
                time_display_str = "未知时间"
                if article_data["time_struct"]:
                    try: time_display_str = time.strftime("%Y-%m-%d", article_data["time_struct"])
                    except: pass
                article_obj_for_storage = {
                    "title": article_data["title"], "url": article_data["url"],
                    "snippet": article_data["snippet"], "source": article_data["source"],
                    "time_display_str": time_display_str, "time_struct": article_data["time_struct"],
                    "source_priority": GOOGLE_NEWS_PRIORITY
                }
                if normalized_title not in unique_articles_candidates or \
                   GOOGLE_NEWS_PRIORITY > unique_articles_candidates[normalized_title]["priority"]:
                    unique_articles_candidates[normalized_title] = {
                        "article_obj": article_obj_for_storage, "priority": GOOGLE_NEWS_PRIORITY,
                        "target_categories": [site_category_name], "url": article_data["url"]
                    }
                    print(f"      候选(Google): '{article_obj_for_storage['title'][:30]}...' (Prio: {GOOGLE_NEWS_PRIORITY}) for category {site_category_name}")
        time.sleep(1)

    # --- 步骤四：将 unique_articles_candidates 分配到最终的分类字典中 ---
    print("\n--- 正在将去重和优先选择后的新闻分配到各分类 ---")
    all_articles_by_site_category_temp = {category_name: [] for category_name in CATEGORIES_CONFIG.keys()}
    for candidate_info in unique_articles_candidates.values():
        article_to_add = candidate_info["article_obj"]
        if article_to_add["url"] not in globally_seen_urls:
            for target_cat in candidate_info["target_categories"]:
                if target_cat in all_articles_by_site_category_temp:
                    all_articles_by_site_category_temp[target_cat].append(article_to_add)
            globally_seen_urls.add(article_to_add["url"])

    # --- 步骤五：对每个分类的文章按日期排序并截取 ---
    print("\n--- 正在对各分类新闻进行排序和截取 ---")
    all_articles_by_site_category_final_sorted = {}
    for category_name, articles_list in all_articles_by_site_category_temp.items():
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
