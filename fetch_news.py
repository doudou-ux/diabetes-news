# fetch_news.py
import datetime
import html
import os
import re
import time
import requests # 使用 requests 进行 HTTP 调用
import json
# 移除了不再需要的 hmac, base64, hashlib, urllib.parse (用于签名)
import feedparser
from bs4 import BeautifulSoup
from urllib.parse import urljoin
# 移除了 deep_translator 的导入
# from deep_translator import GoogleTranslator
try:
    from Bio import Entrez # For PubMed API
except ImportError:
    print("错误：未找到 Biopython 库。请通过 'pip install biopython' 安装。")
    Entrez = None # Set Entrez to None if import fails

# --- (0) 从环境变量读取讯飞星火 API Keys ---
SPARK_API_PASSWORD = os.getenv("SPARK_API_PASSWORD") # 读取 APIPassword
# Spark Lite HTTP Endpoint
SPARK_LITE_HTTP_URL = "https://spark-api-open.xf-yun.com/v1/chat/completions"

# --- (1) 配置权威 RSS 源 ---
# (与 diabetes_news_fetch_all_sources_v2 版本相同)
AUTHORITATIVE_RSS_FEEDS = [
    {"url": "https://www.medscape.com/cx/rss/professional.xml", "source_override": "Medscape Professional", "priority": 10, "needs_translation": True},
    {"url": "https://www.healio.com/sws/feed/news/endocrinology", "source_override": "Healio Endocrinology", "priority": 9, "needs_translation": True},
    {"url": "https://www.diabettech.com/feed/", "source_override": "Diabettech", "priority": 8, "needs_translation": True},
    # {"url": "https://thesavvydiabetic.com/feed/", "source_override": "The Savvy Diabetic", "priority": 7, "needs_translation": True}, # 403
    # {"url": "https://forum.diabetes.org.uk/boards/forums/-/index.rss", "source_override": "Diabetes UK 论坛", "priority": 6, "needs_translation": True}, # Removed by user
    # {"url": "https://www.gov.uk/government/latest.atom?organisations%5B%5D=medicines-and-healthcare-products-regulatory-agency", "source_override": "MHRA (UK)", "priority": 9, "needs_translation": True}, # Atom feed, needs verification
    {"url": "https://www.fda.gov/about-fda/contact-fda/stay-informed/rss-feeds/press-releases/rss.xml", "source_override": "FDA (US) Press", "priority": 10, "needs_translation": True},
    # { "url": "YOUR_PUBMED_RSS_URL", "source_override": "PubMed (RSS Search)", "priority": 12, "needs_translation": True },
    # { "url": "YOUR_ADA_JOURNAL_RSS_URL", "source_override": "Diabetes Care (ADA)", "priority": 11, "needs_translation": True },
]

# --- (1b) 配置爬虫源 ---
# (与 diabetes_news_fetch_all_sources_v2 版本相同)
SCRAPED_SOURCES_CONFIG = [
    {"name": "Breakthrough T1D News", "fetch_function": "fetch_breakthrought1d_articles", "source_override": "Breakthrough T1D", "priority": 8},
    # {"name": "MyGlu Articles", "fetch_function": "fetch_myglu_articles", "source_override": "MyGlu", "priority": 7}, # 404
    {"name": "DZD News (2025)", "fetch_function": "fetch_dzd_articles", "source_override": "DZD News", "priority": 9},
    # {"name": "ADCES News", "fetch_function": "fetch_adces_articles", "source_override": "ADCES News", "priority": 8}, # 404
    # {"name": "PANTHER Program News", "fetch_function": "fetch_panther_articles", "source_override": "PANTHER Program", "priority": 7}, # 404
    # {"name": "NMPA Policies", "fetch_function": "fetch_nmpa_articles", "source_override": "NMPA", "priority": 10}, # 412
    {"name": "PubMed API Search", "fetch_function": "fetch_pubmed_articles", "source_override": "PubMed", "priority": 12},
    {"name": "IDF News", "fetch_function": "fetch_idf_articles", "source_override": "IDF News", "priority": 9},
]

GOOGLE_NEWS_PRIORITY = 1
SOURCE_TYPE_ORDER = {'authoritative_rss': 0, 'scraper': 1, 'google_news': 2, 'unknown': 99}

# --- (2) 配置网站展示的分类 ---
CATEGORIES_CONFIG = {
    "最新研究": {"emoji": "🔬"}, "治疗进展": {"emoji": "💊"}, "饮食与营养": {"emoji": "🥗"},
    "预防与生活方式": {"emoji": "🏃‍♀️"}, "并发症管理": {"emoji": "🩺"}, "患者故事与心理支持": {"emoji": "😊"},
    "政策/医保信息": {"emoji": "📄"}, "综合资讯": {"emoji": "📰"}
}
VALID_CATEGORY_NAMES = list(CATEGORIES_CONFIG.keys())
if "综合资讯" in VALID_CATEGORY_NAMES: VALID_CATEGORY_NAMES.remove("综合资讯")

# --- 帮助函数：规范化标题 ---
def normalize_title(title):
    if not title: return ""
    title = title.lower(); title = re.sub(r'[^\w\s-]', '', title); title = re.sub(r'\s+', ' ', title).strip()
    return title

# --- 帮助函数：判断日期是否在最近一个月内 ---
def is_within_last_month_rss(time_struct, today_date_obj):
    if not time_struct: return False
    try:
        article_date = datetime.date(time_struct.tm_year, time_struct.tm_mon, time_struct.tm_mday)
        thirty_days_ago = today_date_obj - datetime.timedelta(days=30)
        return thirty_days_ago <= article_date <= today_date_obj
    except Exception as e: return False

# --- 帮助函数：清理 HTML ---
def clean_html(raw_html):
    if not raw_html: return ""
    try: return BeautifulSoup(raw_html, "html.parser").get_text(separator=' ', strip=True)
    except Exception: return raw_html

# --- 翻译函数 (使用讯飞星火 HTTP API) ---
def translate_text_with_llm(text, target_lang='Chinese'): # 使用 'Chinese' 作为目标语言提示
    """使用讯飞星火 HTTP API 进行翻译"""
    if not text or not isinstance(text, str) or not text.strip():
        return text # 如果文本为空或非字符串，直接返回
    if not SPARK_API_PASSWORD:
        print("      错误: 讯飞星火 APIPassword 未配置，无法进行 LLM 翻译。")
        return text # 返回原文

    # 限制输入长度，避免超长和过多 token 消耗
    max_input_length = 500
    text_to_translate = text[:max_input_length]

    prompt = f"Please translate the following text to {target_lang}. Only return the translated text, without any introduction or explanation.\n\nText to translate:\n{text_to_translate}"

    print(f"      正在调用讯飞星火 HTTP API 翻译: {text_to_translate[:30]}...")

    try:
        headers = {
            "Content-Type": "application/json",
            "Authorization": f"Bearer {SPARK_API_PASSWORD}"
        }
        payload = {
            "model": "lite",
            "messages": [{"role": "user", "content": prompt}],
            "temperature": 0.3, # 翻译任务通常需要较低的温度以保证准确性
            "max_tokens": int(len(text_to_translate) * 1.5) + 50 # 估算输出 token 上限
        }
        response = requests.post(SPARK_LITE_HTTP_URL, headers=headers, json=payload, timeout=20) # 翻译可能需要更长超时
        response.raise_for_status()
        response_data = response.json()

        llm_output = ""
        if 'choices' in response_data and response_data['choices']:
            message = response_data['choices'][0].get('message', {})
            llm_output = message.get('content', '').strip()
        else:
            print(f"      警告: 未知的讯飞星火 API 翻译响应结构: {response_data}")
            return text # 返回原文

        # 基本的清洗，去除可能的引号
        cleaned_output = llm_output.strip().strip('"').strip("'")

        if cleaned_output:
            print(f"      翻译成功: {text_to_translate[:30]}... -> {cleaned_output[:30]}...")
            return cleaned_output
        else:
            print(f"      警告: 讯飞星火 API 翻译返回为空。")
            return text # 返回原文

    except requests.exceptions.RequestException as req_e:
         print(f"      调用讯飞星火 API 翻译时发生网络或HTTP错误: {req_e}")
         if 'response' in locals() and response is not None:
              print(f"      响应状态码: {response.status_code}")
              try: print(f"      响应内容: {response.json()}")
              except json.JSONDecodeError: print(f"      响应内容 (非JSON): {response.text}")
         return text # 返回原文
    except Exception as e:
        print(f"      调用讯飞星火 API 翻译时出错: {e}")
        return text # 返回原文

# --- (A) RSS 源获取函数 ---
def fetch_articles_from_rss(rss_url, source_name_override=None):
    # (与 diabetes_news_fetch_all_sources_v1 版本相同)
    print(f"    正在从 RSS 源获取: {rss_url} ({source_name_override or '未知源'})")
    articles = []
    try:
        headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'}
        response = requests.get(rss_url, headers=headers, timeout=20, allow_redirects=True)
        if response.url != rss_url and "apology" in response.url:
             print(f"      请求被重定向到错误页面: {response.url}")
             response.raise_for_status() 
        response.raise_for_status()
        feed = feedparser.parse(response.content)
        if not feed.entries:
            if feed.bozo: print(f"      警告: feedparser 解析 RSS 源时遇到问题: {feed.bozo_exception} (URL: {rss_url})")
            else: print(f"      此 RSS 源未返回任何条目: {rss_url}")
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
# 现在爬虫函数内部不再调用 translate_text，翻译将在主逻辑中处理
def fetch_breakthrought1d_articles():
    BASE_URL = "https://www.breakthrought1d.org/news/"
    print(f"    正在爬取: {BASE_URL}")
    articles = []
    try:
        headers = {"User-Agent": "Mozilla/5.0"}
        response = requests.get(BASE_URL, headers=headers, timeout=20)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, "html.parser")
        for article_el in soup.select("article"):
            a_tag = article_el.select_one("h2 a")
            if not a_tag: continue
            title_en = a_tag.get_text(strip=True)
            link = urljoin(BASE_URL, a_tag.get("href"))
            summary_tag = article_el.select_one("p")
            summary_en = summary_tag.get_text(strip=True) if summary_tag else ""
            time_struct = None
            if not time_struct: print(f"      警告: 未能从 {link} 提取发布日期。")
            # 不再翻译，返回原文
            articles.append({
                "title": title_en, "url": link, "snippet": summary_en,
                "source": "Breakthrough T1D", "time_struct": time_struct,
                "needs_translation": True # 添加标记
            })
    except Exception as e: print(f"      爬取 Breakthrough T1D 时出错: {e}")
    return articles

def fetch_myglu_articles():
    print("    跳过 MyGlu 爬虫 (已注释掉)")
    return []

def fetch_dzd_articles():
    BASE_URL = "https://www.dzd-ev.de/en/press/press-releases/press-releases-2025/index.html"
    print(f"    正在爬取: {BASE_URL}")
    articles = []
    try:
        headers = {"User-Agent": "Mozilla/5.0"}
        response = requests.get(BASE_URL, headers=headers, timeout=20)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, "html.parser")
        for item in soup.select("div.teaser-text"):
            title_tag = item.select_one("a")
            if not title_tag: continue
            title_en = title_tag.get_text(strip=True)
            link = urljoin(BASE_URL, title_tag.get("href"))
            summary_tag = item.select_one("p")
            summary_en = summary_tag.get_text(strip=True) if summary_tag else ""
            time_struct = None
            if not time_struct: print(f"      警告: 未能从 {link} 提取发布日期。")
            # 不再翻译，返回原文
            articles.append({
                "title": title_en, "url": link, "snippet": summary_en,
                "source": "DZD News", "time_struct": time_struct,
                "needs_translation": True # 添加标记
            })
    except Exception as e: print(f"      爬取 DZD News 时出错: {e}")
    return articles

def fetch_adces_articles():
    print("    跳过 ADCES News 爬虫 (已注释掉)")
    return []

def fetch_panther_articles():
    print("    跳过 PANTHER Program 爬虫 (已注释掉)")
    return []

def fetch_nmpa_articles():
    print("    跳过 NMPA 爬虫 (已注释掉)")
    return []

def fetch_pubmed_articles():
    if not Entrez: return []
    Entrez.email = os.getenv("PUBMED_API_EMAIL", "default_email@example.com")
    search_term = "(diabetes[Title/Abstract]) AND (treatment[Title/Abstract] OR research[Title/Abstract] OR prevention[Title/Abstract])"
    print(f"    正在通过 PubMed API 搜索: {search_term}")
    articles = []
    MAX_PUBMED_RESULTS = 20
    try:
        handle_search = Entrez.esearch(db="pubmed", term=search_term, retmax=MAX_PUBMED_RESULTS, sort="pub date")
        record_search = Entrez.read(handle_search)
        handle_search.close()
        id_list = record_search["IdList"]
        if not id_list: return []
        print(f"      PubMed API 返回 {len(id_list)} 个文章 ID，正在获取摘要...")
        handle_summary = Entrez.esummary(db="pubmed", id=",".join(id_list))
        summary_record = Entrez.read(handle_summary)
        handle_summary.close()
        for docsum in summary_record:
            title_en = docsum.get("Title", "No Title Available")
            pmid = docsum.get("Id", "")
            link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "javascript:void(0);"
            snippet_en = docsum.get("Source", "No Summary Available")
            authors = docsum.get("AuthorList", [])
            pubdate_str = docsum.get("PubDate", "")
            time_struct = None
            if pubdate_str:
                try:
                    dt_obj = None
                    try: dt_obj = datetime.datetime.strptime(pubdate_str, "%Y %b %d")
                    except ValueError: pass
                    if not dt_obj:
                        try: dt_obj = datetime.datetime.strptime(pubdate_str, "%Y %b")
                        except ValueError: pass
                    if not dt_obj:
                        try: dt_obj = datetime.datetime.strptime(pubdate_str, "%Y")
                        except ValueError: pass
                    if dt_obj: time_struct = dt_obj.timetuple()
                except Exception as date_e: print(f"      解析 PubMed 日期时出错: {date_e} - {pubdate_str}")
            # 不再翻译，返回原文
            articles.append({
                "title": title_en, "url": link, "snippet": f"{snippet_en} (作者: {', '.join(authors)[:50]}...)",
                "source": "PubMed", "time_struct": time_struct,
                "needs_translation": True # 添加标记
            })
            time.sleep(0.4)
    except Exception as e: print(f"      处理 PubMed API 时出错: {e}")
    return articles

def fetch_idf_articles():
    url = "https://www.idf.org/news"
    print(f"    正在爬取: {url}")
    articles = []
    try:
        headers = {"User-Agent": "Mozilla/5.0"}
        response = requests.get(url, headers=headers, timeout=20)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, 'html.parser')
        for item in soup.select("article.news-item, div.news-item"):
            title_tag = item.select_one("h3 a, h2 a, .title a")
            if not title_tag: continue
            title_en = title_tag.get_text(strip=True)
            link = urljoin(url, title_tag.get("href"))
            summary_tag = item.select_one("p, .summary, .excerpt")
            summary_en = summary_tag.get_text(strip=True) if summary_tag else ""
            time_struct = None
            if not time_struct: print(f"      警告: 未能从 {link} 提取发布日期。")
            # 不再翻译，返回原文
            articles.append({
                "title": title_en, "url": link, "snippet": summary_en,
                "source": "IDF News", "time_struct": time_struct,
                "needs_translation": True # 添加标记
            })
    except Exception as e: print(f"      爬取 IDF News 时出错: {e}")
    return articles

SCRAPER_FUNCTIONS_MAP = {
    "fetch_breakthrought1d_articles": fetch_breakthrought1d_articles,
    "fetch_myglu_articles": fetch_myglu_articles,
    "fetch_dzd_articles": fetch_dzd_articles,
    "fetch_adces_articles": fetch_adces_articles,
    "fetch_panther_articles": fetch_panther_articles,
    "fetch_nmpa_articles": fetch_nmpa_articles,
    "fetch_pubmed_articles": fetch_pubmed_articles,
    "fetch_idf_articles": fetch_idf_articles,
}

# --- (C) 使用讯飞星火 HTTP API 进行动态分类 ---
def categorize_article_with_llm(article_obj):
    # (与 diabetes_news_fetch_llm_http_categorization_v2 版本相同)
    if not SPARK_API_PASSWORD:
        print("      错误: 讯飞星火 APIPassword 未配置，无法进行 LLM 分类。将归入'综合资讯'。")
        return "综合资讯"
    title = article_obj.get("title", "")
    snippet = article_obj.get("snippet", "")
    text_to_classify = f"标题：{title}\n摘要：{snippet[:300]}"
    prompt = f"""请根据以下文章内容，判断它最符合下列哪个分类？请严格从列表中选择一个，并只返回分类名称，不要添加任何其他解释或文字。

可选分类列表：{', '.join(VALID_CATEGORY_NAMES)}

文章内容：
{text_to_classify}

最合适的分类名称是："""
    print(f"      正在调用讯飞星火 HTTP API 对 '{title[:30]}...' 进行分类...")
    try:
        headers = {"Content-Type": "application/json", "Authorization": f"Bearer {SPARK_API_PASSWORD}"}
        payload = {"model": "lite", "messages": [{"role": "user", "content": prompt}], "temperature": 0.5, "max_tokens": 50}
        response = requests.post(SPARK_LITE_HTTP_URL, headers=headers, json=payload, timeout=30)
        response.raise_for_status()
        response_data = response.json()
        llm_output = ""
        if 'choices' in response_data and response_data['choices']:
            message = response_data['choices'][0].get('message', {})
            llm_output = message.get('content', '').strip()
        else:
            print(f"      警告: 未知的讯飞星火 API 翻译响应结构: {response_data}")
            llm_output = ""
        print(f"      讯飞星火 API 返回: '{llm_output}'")
        cleaned_output = llm_output.strip().replace("\"", "").replace("'", "").replace("：","").replace(":","")
        if cleaned_output in VALID_CATEGORY_NAMES:
            print(f"      文章 '{title[:30]}...' 成功分类到 '{cleaned_output}'")
            return cleaned_output
        else:
            print(f"      警告: 讯飞星火 API 返回的分类 '{cleaned_output}' (原始: '{llm_output}') 无效或不在列表中。将归入'综合资讯'。")
            return "综合资讯"
    except requests.exceptions.RequestException as req_e:
         print(f"      调用讯飞星火 API 时发生网络或HTTP错误: {req_e}")
         if 'response' in locals() and response is not None:
              print(f"      响应状态码: {response.status_code}")
              try: print(f"      响应内容: {response.json()}")
              except json.JSONDecodeError: print(f"      响应内容 (非JSON): {response.text}")
         return "综合资讯"
    except Exception as e:
        print(f"      调用讯飞星火 API 时出错: {e}")
        return "综合资讯"

# --- HTML 生成逻辑 ---
def generate_html_content(all_news_data_sorted):
    # (与 diabetes_news_fetch_tabs_v1 中的 generate_html_content 完全相同)
    # ... (省略 HTML 生成代码) ...
    current_time_str = datetime.datetime.now().strftime('%Y年%m月%d日 %H:%M:%S')
    app_timezone = os.getenv('APP_TIMEZONE', 'UTC')
    if app_timezone != 'UTC':
        try:
             current_time_str = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).strftime('%Y年%m月%d日 %H:%M:%S %Z')
        except Exception as e:
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
            <p class="text-sm text-gray-500 mt-2">资讯综合来源 (由 AI 智能分类)</p>
        </header>
        <div class="tab-buttons-container" id="tabButtons">"""
    first_category = True
    for category_name_key in CATEGORIES_CONFIG.keys():
        category_config = CATEGORIES_CONFIG.get(category_name_key, {})
        emoji = category_config.get("emoji", "")
        tab_id = "tab-" + html.escape(category_name_key.replace(" ", "-").replace("/", "-").lower())
        is_active = False
        if first_category:
             if category_name_key in all_news_data_sorted and all_news_data_sorted[category_name_key]:
                 is_active = True
                 first_category = False 
             elif list(CATEGORIES_CONFIG.keys())[0] == category_name_key and not any(all_news_data_sorted.values()):
                  is_active = True 
        active_class = "active" if is_active else ""
        html_output += f"""<button class="tab-button {active_class}" data-tab-target="#{tab_id}">{emoji} {html.escape(category_name_key)}</button>"""

    html_output += """</div><div id="news-content">"""
    first_category_content = True 
    if not any(all_news_data_sorted.values()):
        html_output += '<p class="text-center text-gray-500 text-xl py-10">抱歉，目前未能加载到最近一个月相关的糖尿病资讯。</p>'
    else:
        for category_name_key in CATEGORIES_CONFIG.keys():
            articles = all_news_data_sorted.get(category_name_key, []) 
            category_config = CATEGORIES_CONFIG.get(category_name_key, {})
            emoji = category_config.get("emoji", "")
            tab_id = "tab-" + html.escape(category_name_key.replace(" ", "-").replace("/", "-").lower())
            is_active_content = False
            if first_category_content:
                 if articles: 
                      is_active_content = True
                      first_category_content = False 
                 elif list(CATEGORIES_CONFIG.keys())[0] == category_name_key and not any(all_news_data_sorted.values()):
                      is_active_content = True 
            active_class = "active" if is_active_content else ""
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

    # --- 注意：f-string 中的 JavaScript 大括号需要转义 ---
    html_output += f"""
        </div> </div> <footer class="text-center p-6 mt-12 text-gray-600 text-sm border-t border-gray-300"><p>&copy; {current_year} 糖尿病资讯聚合. <a href="{github_repo_url}" target="_blank" rel="noopener noreferrer" class="text-blue-600 hover:underline">项目源码</a></p><p class="mt-1">本站内容仅供参考, 不构成医疗建议。</p></footer>
    <script>
        document.addEventListener('DOMContentLoaded', function () {{ // JS 大括号需要转义
            const tabButtons = document.querySelectorAll('.tab-button');
            const tabContents = document.querySelectorAll('.tab-content');
            if (tabButtons.length > 0 && tabContents.length > 0) {{
                let firstActiveButton = document.querySelector('.tab-button.active');
                let firstActiveContent = document.querySelector('.tab-content.active');

                if (!firstActiveButton && tabButtons.length > 0) {{
                    tabButtons[0].classList.add('active'); 
                    firstActiveButton = tabButtons[0];
                }}
                if (!firstActiveContent && firstActiveButton) {{
                    const targetId = firstActiveButton.dataset.tabTarget;
                    const targetContent = document.querySelector(targetId);
                    if (targetContent) {{
                        targetContent.classList.add('active'); 
                    }} else if (tabContents.length > 0) {{
                         tabContents[0].classList.add('active');
                    }}
                }} else if (firstActiveButton && firstActiveContent) {{
                     if (firstActiveButton.dataset.tabTarget !== '#' + firstActiveContent.id) {{
                          tabContents.forEach(content => content.classList.remove('active'));
                          const targetContent = document.querySelector(firstActiveButton.dataset.tabTarget);
                          if (targetContent) targetContent.classList.add('active');
                     }}
                }}

                tabButtons.forEach(button => {{
                    button.addEventListener('click', () => {{ // JS 大括号需要转义
                        tabButtons.forEach(btn => btn.classList.remove('active'));
                        tabContents.forEach(content => content.classList.remove('active'));
                        button.classList.add('active');
                        const targetContentId = button.dataset.tabTarget;
                        const targetContent = document.querySelector(targetContentId);
                        if (targetContent) {{ targetContent.classList.add('active'); }} // JS 大括号需要转义
                    }}); // JS 大括号需要转义
                }}); // JS 大括号需要转义
            }}
        }}); // JS 大括号需要转义
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
    MAX_LLM_CALLS = 100 # 适当增加 LLM 调用次数上限 (翻译+分类)
    llm_call_count = 0

    # --- 步骤一：从权威 RSS 源获取新闻 ---
    print("\n--- 正在从权威 RSS 源获取新闻 ---")
    for feed_info in AUTHORITATIVE_RSS_FEEDS:
        current_priority = feed_info.get("priority", 5) 
        needs_translation = feed_info.get("needs_translation", False)
        raw_articles_from_feed = fetch_articles_from_rss(feed_info["url"], feed_info["source_override"])
        for article_data in raw_articles_from_feed:
            if article_data["url"] in globally_seen_urls: continue
            
            title_to_process = article_data["title"]
            snippet_to_process = article_data["snippet"]
            
            # 调用 LLM 进行翻译 (如果需要)
            if needs_translation and SPARK_API_PASSWORD and llm_call_count < MAX_LLM_CALLS:
                print(f"    需要翻译来自 {feed_info['source_override']} 的文章: {title_to_process[:30]}...")
                translated_title = translate_text_with_llm(title_to_process)
                llm_call_count += 1 # 计入调用次数
                time.sleep(1.1) # 调用间隔
                translated_snippet = translate_text_with_llm(snippet_to_process)
                llm_call_count += 1
                time.sleep(1.1)
                # 只有当翻译成功时才替换原文 (translate_text_with_llm 失败时返回原文)
                title_to_process = translated_title
                snippet_to_process = translated_snippet
            elif needs_translation:
                 print(f"    警告: 无法翻译来自 {feed_info['source_override']} 的文章，API Key缺失或达到调用上限。")


            if is_within_last_month_rss(article_data["time_struct"], today):
                normalized_title = normalize_title(title_to_process) # 使用处理后的标题
                time_display_str = "未知时间"
                if article_data["time_struct"]:
                    try: time_display_str = time.strftime("%Y-%m-%d", article_data["time_struct"])
                    except: pass
                
                article_obj_for_storage = {
                    "title": title_to_process, "url": article_data["url"], "snippet": snippet_to_process, 
                    "source": article_data["source"], "time_display_str": time_display_str, 
                    "time_struct": article_data["time_struct"], 
                    "source_priority": current_priority, "source_type": "authoritative_rss"
                }
                if normalized_title not in unique_articles_candidates or \
                   current_priority > unique_articles_candidates[normalized_title]["priority"]:
                    unique_articles_candidates[normalized_title] = {
                        "article_obj": article_obj_for_storage, "priority": current_priority,
                        "url": article_data["url"]
                    }
        time.sleep(1)

    # --- 步骤二：从爬虫源获取新闻 ---
    print("\n--- 正在从爬虫源获取新闻 ---")
    for scraper_info in SCRAPED_SOURCES_CONFIG:
        if scraper_info["fetch_function"] not in SCRAPER_FUNCTIONS_MAP: continue
        fetch_function = SCRAPER_FUNCTIONS_MAP[scraper_info["fetch_function"]]
        current_priority = scraper_info.get("priority", 3)
        raw_articles_from_scraper = fetch_function() # 爬虫函数内部不再翻译
        
        for article_data in raw_articles_from_scraper: 
            if article_data["url"] in globally_seen_urls: continue

            title_to_process = article_data["title"]
            snippet_to_process = article_data["snippet"]
            needs_translation = article_data.get("needs_translation", False) # 检查爬虫是否标记需要翻译

            # 调用 LLM 进行翻译 (如果需要)
            if needs_translation and SPARK_API_PASSWORD and llm_call_count < MAX_LLM_CALLS:
                print(f"    需要翻译来自 {scraper_info['source_override']} 的文章: {title_to_process[:30]}...")
                translated_title = translate_text_with_llm(title_to_process)
                llm_call_count += 1
                time.sleep(1.1)
                translated_snippet = translate_text_with_llm(snippet_to_process)
                llm_call_count += 1
                time.sleep(1.1)
                title_to_process = translated_title
                snippet_to_process = translated_snippet
            elif needs_translation:
                print(f"    警告: 无法翻译来自 {scraper_info['source_override']} 的文章，API Key缺失或达到调用上限。")


            if is_within_last_month_rss(article_data["time_struct"], today):
                normalized_title = normalize_title(title_to_process)
                time_display_str = "未知时间"
                if article_data.get("time_struct"):
                    try: time_display_str = time.strftime("%Y-%m-%d", article_data["time_struct"])
                    except: pass
                else: pass 
                article_obj_for_storage = {
                    "title": title_to_process, "url": article_data["url"], "snippet": snippet_to_process, 
                    "source": scraper_info["source_override"], "time_display_str": time_display_str, 
                    "time_struct": article_data["time_struct"], 
                    "source_priority": current_priority, "source_type": "scraper"
                }
                if normalized_title not in unique_articles_candidates or \
                   current_priority > unique_articles_candidates[normalized_title]["priority"]:
                    unique_articles_candidates[normalized_title] = {
                        "article_obj": article_obj_for_storage, "priority": current_priority,
                        "url": article_data["url"]
                    }
        time.sleep(1)

    # --- 步骤三：从 Google News RSS 获取补充新闻 ---
    print("\n--- 正在从 Google News RSS 获取补充新闻 (用于全局候选池) ---")
    google_search_term = "糖尿病 新闻 OR diabetes news" 
    print(f"  使用 Google News 搜索词: {google_search_term}")
    google_news_rss_url = f"https://news.google.com/rss/search?q={html.escape(google_search_term)}&hl=zh-CN&gl=CN&ceid=CN:zh-Hans"
    raw_articles_from_google = fetch_articles_from_rss(google_news_rss_url, source_name_override=None)
    for article_data in raw_articles_from_google:
        # Google News 默认是中文，不需要翻译
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
                "title": article_data["title"], "url": article_data["url"], "snippet": article_data["snippet"], 
                "source": article_data["source"], "time_display_str": time_display_str, 
                "time_struct": article_data["time_struct"], 
                "source_priority": GOOGLE_NEWS_PRIORITY, "source_type": "google_news"
            }
            if normalized_title not in unique_articles_candidates or \
               GOOGLE_NEWS_PRIORITY > unique_articles_candidates[normalized_title]["priority"]:
                unique_articles_candidates[normalized_title] = {
                    "article_obj": article_obj_for_storage, "priority": GOOGLE_NEWS_PRIORITY,
                    "url": article_data["url"]
                }
    time.sleep(1)

    # --- 步骤四：使用 LLM 动态分类所有候选文章 ---
    print("\n--- 正在对所有候选文章进行动态分类 (使用讯飞星火 HTTP API) ---")
    all_articles_by_site_category_temp = {category_name: [] for category_name in CATEGORIES_CONFIG.keys()}
    categorized_urls = set() 
    spark_api_ready = bool(SPARK_API_PASSWORD) 
    if not spark_api_ready:
        print("警告: 讯飞星火 APIPassword 未配置在环境变量中，将跳过 LLM 分类，所有文章归入'综合资讯'。")
    
    for candidate_info in unique_articles_candidates.values():
        article_to_categorize = candidate_info["article_obj"]
        article_url = article_to_categorize["url"]
        if article_url in categorized_urls: continue

        best_category = "综合资讯" 
        if spark_api_ready and llm_call_count < MAX_LLM_CALLS:
            try:
                best_category = categorize_article_with_llm(article_to_categorize)
                llm_call_count += 1
                time.sleep(1.1) 
            except Exception as llm_e:
                print(f"    LLM 分类时发生意外错误: {llm_e}，文章将归入'综合资讯'。")
                best_category = "综合资讯"
        elif llm_call_count >= MAX_LLM_CALLS:
             print(f"    已达到 LLM 调用次数上限 ({MAX_LLM_CALLS})，剩余文章将归入'综合资讯'。")
             best_category = "综合资讯"
        elif not spark_api_ready:
             best_category = "综合资讯" 

        if best_category in all_articles_by_site_category_temp:
            all_articles_by_site_category_temp[best_category].append(article_to_categorize)
            categorized_urls.add(article_url) 
        else:
            print(f"    警告: LLM 返回的分类 '{best_category}' 不在预设分类中，归入 '综合资讯'")
            all_articles_by_site_category_temp["综合资讯"].append(article_to_categorize)
            categorized_urls.add(article_url)

    # --- 步骤五：对每个分类的文章按来源类型和日期排序并截取 ---
    print("\n--- 正在对各分类新闻进行排序和截取 ---")
    all_articles_by_site_category_final_sorted = {}
    for category_name, articles_list in all_articles_by_site_category_temp.items():
        articles_list.sort(key=lambda x: (
            SOURCE_TYPE_ORDER.get(x.get("source_type", "unknown"), 99),
            -(time.mktime(x["time_struct"]) if x.get("time_struct") else -float('inf'))
        ))
        all_articles_by_site_category_final_sorted[category_name] = articles_list[:MAX_ARTICLES_PER_CATEGORY]
        print(f"  分类 '{category_name}' 排序并截取后有 {len(all_articles_by_site_category_final_sorted[category_name])} 条新闻。")
        if all_articles_by_site_category_final_sorted[category_name]:
            print(f"    排序后前几条来源类型: {[a.get('source_type', 'unknown') for a in all_articles_by_site_category_final_sorted[category_name][:5]]}")

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
