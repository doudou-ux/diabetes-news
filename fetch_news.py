import datetime
import html
import os
import re
import time
import requests # 使用 requests 进行 HTTP 调用
import json
import feedparser
from bs4 import BeautifulSoup
from urllib.parse import urljoin, quote_plus
from dateutil.parser import parse as dateutil_parse # 使用 dateutil 解析日期更灵活

try:
    from Bio import Entrez # For PubMed API
except ImportError:
    print("错误：未找到 Biopython 库。请通过 'pip install biopython' 安装。")
    Entrez = None # Set Entrez to None if import fails

# --- (0) 从环境变量读取讯飞星火 API Keys ---
SPARK_API_PASSWORD = os.getenv("SPARK_API_PASSWORD")
SPARK_LITE_HTTP_URL = "https://spark-api-open.xf-yun.com/v1/chat/completions"
SPARK_APPID = os.getenv("SPARK_APPID")
SPARK_API_KEY = os.getenv("SPARK_API_KEY")
SPARK_API_SECRET = os.getenv("SPARK_API_SECRET")

# --- (1) 配置权威 RSS 源 ---
AUTHORITATIVE_RSS_FEEDS = [
    # {"url": "https://www.medscape.com/cx/rss/professional.xml", "source_override": "Medscape Professional", "priority": 10, "needs_translation": True}, # Original - 404
    {"url": "http://rss.medscape.com/medscapetoday.rss", "source_override": "Medscape Today", "priority": 10, "needs_translation": True}, # 尝试新的 Medscape 链接 (可能不稳定，易出现DNS解析错误)
    {"url": "https://www.healio.com/sws/feed/news/endocrinology", "source_override": "Healio Endocrinology", "priority": 9, "needs_translation": True},
    {"url": "https://www.diabettech.com/feed/", "source_override": "Diabettech", "priority": 8, "needs_translation": True},
    {"url": "https://www.fda.gov/about-fda/contact-fda/stay-informed/rss-feeds/press-releases/rss.xml", "source_override": "FDA (US) Press", "priority": 10, "needs_translation": True}, # FDA 链接可能因服务器阻止而失败
]

# --- (1b) 配置爬虫源 ---
SCRAPED_SOURCES_CONFIG = [
    {"name": "Breakthrough T1D News", "fetch_function": "fetch_breakthrought1d_articles", "source_override": "Breakthrough T1D", "priority": 8}, # 爬虫可能需要更新选择器
    {"name": "DZD News (2025)", "fetch_function": "fetch_dzd_articles", "source_override": "DZD News", "priority": 9}, # 爬虫可能需要更新选择器和年份
    {"name": "PubMed API Search", "fetch_function": "fetch_pubmed_articles", "source_override": "PubMed", "priority": 12},
    {"name": "IDF News", "fetch_function": "fetch_idf_articles", "source_override": "IDF News", "priority": 9}, # 爬虫可能需要更新选择器
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

# --- 帮助函数：尝试解析日期字符串 ---
def parse_date_flexible(date_str):
    """尝试用多种方式解析日期字符串，返回 time.struct_time 或 None"""
    if not date_str:
        return None
    try:
        # dateutil.parser 可以处理多种格式
        dt_obj = dateutil_parse(date_str)
        return dt_obj.timetuple()
    except Exception as e:
        # print(f"      使用 dateutil 解析日期失败: {date_str} - {e}") # 可以取消注释以调试日期解析
        return None # 如果所有尝试都失败

# --- 帮助函数：判断日期是否在最近一个月内 (使用 time.struct_time) ---
# 这个函数现在不再用于过滤，但保留以备将来可能需要
def is_within_last_month(time_struct, today_date_obj):
    """检查 time.struct_time 对象是否在最近30天内"""
    if not time_struct:
        return False
    try:
        if not all(hasattr(time_struct, attr) for attr in ['tm_year', 'tm_mon', 'tm_mday']):
             print(f"      日期结构不完整: {time_struct}，无法判断是否在最近一月内。")
             return False

        article_date = datetime.date(time_struct.tm_year, time_struct.tm_mon, time_struct.tm_mday)
        thirty_days_ago = today_date_obj - datetime.timedelta(days=30)
        is_recent = thirty_days_ago <= article_date <= today_date_obj
        return is_recent
    except ValueError as ve:
        print(f"      解析日期时发生值错误: {time_struct} - {ve}")
        return False
    except Exception as e:
        print(f"      判断日期时发生未知错误: {time_struct} - {e}")
        return False

# --- 帮助函数：清理 HTML ---
def clean_html(raw_html):
    if not raw_html: return ""
    try: return BeautifulSoup(raw_html, "html.parser").get_text(separator=' ', strip=True)
    except Exception: return raw_html

# --- 翻译函数 (使用讯飞星火 HTTP API - Bearer Token Auth) ---
def translate_text_with_llm(text, target_lang='Chinese'):
    global llm_call_count, MAX_LLM_CALLS # 引用全局变量
    if not text or not isinstance(text, str) or not text.strip(): return text
    if not SPARK_API_PASSWORD:
        print("      错误: 讯飞星火 APIPassword 未配置，无法进行 LLM 翻译。")
        return text
    if llm_call_count >= MAX_LLM_CALLS:
        return text # 返回原文

    max_input_length = 500
    text_to_translate = text[:max_input_length]
    prompt = f"Please translate the following text to {target_lang}. Only return the translated text, without any introduction or explanation.\n\nText to translate:\n{text_to_translate}"
    llm_call_count += 1

    try:
        headers = {"Content-Type": "application/json", "Authorization": f"Bearer {SPARK_API_PASSWORD}"}
        payload = {"model": "lite", "messages": [{"role": "user", "content": prompt}], "temperature": 0.3, "max_tokens": int(len(text_to_translate) * 1.5) + 50}
        response = requests.post(SPARK_LITE_HTTP_URL, headers=headers, json=payload, timeout=20)
        response.raise_for_status()
        response_data = response.json()
        llm_output = ""
        if 'choices' in response_data and response_data['choices']:
            message = response_data['choices'][0].get('message', {})
            llm_output = message.get('content', '').strip()
        elif 'payload' in response_data and 'choices' in response_data['payload'] and \
             'text' in response_data['payload']['choices'] and response_data['payload']['choices']['text']:
            llm_output = response_data['payload']['choices']['text'][0].get('content', '').strip()
        else:
            print(f"      警告: 未知的讯飞星火 API 翻译响应结构: {response_data}")
            return text
        
        cleaned_output = llm_output.strip().strip('"').strip("'")
        if cleaned_output and cleaned_output.lower() != text_to_translate.lower():
            return cleaned_output
        else:
            return text # 返回原文
    except requests.exceptions.RequestException as req_e:
        print(f"      调用讯飞星火 API 翻译时发生网络或HTTP错误: {req_e}")
        if 'response' in locals() and response is not None:
            print(f"      响应状态码: {response.status_code}")
            try: print(f"      响应内容: {response.json()}")
            except json.JSONDecodeError: print(f"      响应内容 (非JSON): {response.text}")
        return text
    except Exception as e:
        print(f"      调用讯飞星火 API 翻译时出错: {e}")
        return text

# --- (A) RSS 源获取函数 ---
def fetch_articles_from_rss(rss_url, source_name_override=None):
    print(f"    正在从 RSS 源获取: {rss_url} ({source_name_override or '未知源'})")
    articles = []
    try:
        headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/100.0.4896.127 Safari/537.36'}
        response = requests.get(rss_url, headers=headers, timeout=25, allow_redirects=True)

        if response.url != rss_url and ("apology" in response.url or "abuse" in response.url or "error" in response.url):
            print(f"      请求被重定向到疑似错误/验证页面: {response.url}")
            if response.status_code >= 400:
                 print(f"      重定向后页面状态码为 {response.status_code}，判定为失败。")
                 response.raise_for_status()
            else:
                 print(f"      重定向后页面状态码为 {response.status_code}，尝试继续解析。")
        else:
             response.raise_for_status()

        feed = feedparser.parse(response.content)
        if feed.bozo:
             print(f"      警告: feedparser 解析 RSS 源时遇到问题 (可能不严重): {feed.bozo_exception} (URL: {rss_url})")
        if not feed.entries:
            print(f"      此 RSS 源未返回任何条目: {rss_url}")
            return []

        print(f"      从此 RSS 源原始获取到 {len(feed.entries)} 条新闻。")
        for entry in feed.entries:
            title = entry.get("title", "无标题")
            link = entry.get("link", f"javascript:void(0);_{html.escape(title)}")

            published_time_struct = entry.get("published_parsed") or entry.get("updated_parsed")
            if not published_time_struct:
                date_str = entry.get("published") or entry.get("updated") or entry.get("dc_date")
                if date_str:
                    published_time_struct = parse_date_flexible(date_str)

            summary_html = entry.get("summary", entry.get("description", "暂无摘要"))
            snippet = clean_html(summary_html)
            actual_source_name = source_name_override
            if not actual_source_name:
                source_info = entry.get("source")
                actual_source_name = source_info.get("title") if source_info else "未知来源"
                if "news.google.com" in rss_url and not source_name_override:
                    if ' - ' in title:
                        parts = title.rsplit(' - ', 1)
                        if len(parts) == 2 and parts[1]:
                             actual_source_name = parts[1]
                             title = parts[0]
            articles.append({
                "title": title.strip(),
                "url": link,
                "snippet": snippet,
                "source": actual_source_name,
                "time_struct": published_time_struct
            })
    except requests.exceptions.Timeout: print(f"      获取 RSS 源时发生超时错误: {rss_url}")
    except requests.exceptions.SSLError as ssl_e: print(f"      获取 RSS 源时发生 SSL 错误: {ssl_e} (URL: {rss_url})")
    except requests.exceptions.RequestException as e: print(f"      获取 RSS 源时发生网络错误: {e} (URL: {rss_url})")
    except Exception as e: print(f"      处理 RSS 源时发生未知错误: {e} (URL: {rss_url})")
    return articles

# --- (B) 爬虫函数定义 ---

def fetch_breakthrought1d_articles():
    BASE_URL = "https://www.breakthrought1d.org/news/"
    print(f"    正在爬取: {BASE_URL}")
    articles = []
    try:
        headers = {"User-Agent": "Mozilla/5.0"}
        response = requests.get(BASE_URL, headers=headers, timeout=20)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, "html.parser")
        for article_el in soup.select("article.post, div.news-item"):
            a_tag = article_el.select_one("h2 a, h3 a, .entry-title a, .news-title a")
            if not a_tag: continue
            title_en = a_tag.get_text(strip=True)
            link = urljoin(BASE_URL, a_tag.get("href"))
            summary_tag = article_el.select_one("div.entry-summary p, div.post-excerpt p, .news-summary p")
            summary_en = summary_tag.get_text(strip=True) if summary_tag else ""

            time_struct = None
            date_tag = article_el.select_one("time.published, span.posted-on time, .post-date, .news-date")
            if date_tag:
                date_str = date_tag.get_text(strip=True) or date_tag.get('datetime')
                if date_str:
                    time_struct = parse_date_flexible(date_str)

            articles.append({
                "title": title_en, "url": link, "snippet": summary_en,
                "source": "Breakthrough T1D", "time_struct": time_struct,
                "needs_translation": True
            })
    except Exception as e: print(f"      爬取 Breakthrough T1D 时出错: {e}")
    return articles

def fetch_dzd_articles():
    BASE_URL = "https://www.dzd-ev.de/en/press/press-releases/press-releases-2025/index.html"
    print(f"    正在爬取: {BASE_URL}")
    articles = []
    try:
        headers = {"User-Agent": "Mozilla/5.0"}
        response = requests.get(BASE_URL, headers=headers, timeout=20)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, "html.parser")

        for item in soup.select("div.news-list-item, div.teaser"):
            title_tag = item.select_one("h3 a, h2 a, .news-title a")
            if not title_tag: continue
            title_en = title_tag.get_text(strip=True)
            link = urljoin(BASE_URL, title_tag.get("href"))

            summary_tag = item.select_one("p, .teaser-text p, .news-teaser p")
            summary_en = summary_tag.get_text(strip=True) if summary_tag else ""

            time_struct = None
            date_span = item.select_one("span.date, div.date, .news-date")
            if date_span:
                date_str = date_span.get_text(strip=True)
                time_struct = parse_date_flexible(date_str)

            articles.append({
                "title": title_en, "url": link, "snippet": summary_en,
                "source": "DZD News", "time_struct": time_struct,
                "needs_translation": True
            })
    except Exception as e: print(f"      爬取 DZD News 时出错: {e}")
    return articles

def fetch_pubmed_articles():
    if not Entrez:
        print("    错误: Biopython (Entrez) 未加载，跳过 PubMed API 调用。")
        return []
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
        if not id_list:
            print("      PubMed API 未返回任何文章 ID。")
            return []
        
        print(f"      PubMed API 返回 {len(id_list)} 个文章 ID，正在获取摘要...")
        handle_summary = Entrez.esummary(db="pubmed", id=",".join(id_list))
        summary_record = Entrez.read(handle_summary)
        handle_summary.close()

        for docsum in summary_record:
            title_en = docsum.get("Title", "No Title Available")
            pmid = docsum.get("Id", "")
            link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "javascript:void(0);"
            journal_source = docsum.get("Source", "")
            authors = docsum.get("AuthorList", [])
            authors_str = ", ".join(authors[:3]) + (" et al." if len(authors) > 3 else "") if authors else "N/A"
            snippet_en = f"Journal: {journal_source}. Authors: {authors_str}." if journal_source or authors else "No Summary Available"

            pubdate_str = docsum.get("PubDate", "")
            time_struct = None
            if pubdate_str:
                time_struct = parse_date_flexible(pubdate_str)

            articles.append({
                "title": title_en, "url": link, "snippet": snippet_en,
                "source": "PubMed", "time_struct": time_struct,
                "needs_translation": True
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

        for item in soup.select("article.news-item, div.news-item, li.news-list-item, .card"):
            title_tag = item.select_one("h3 a, h2 a, .title a, .news-title a, .card-title a")
            if not title_tag: continue
            
            title_en = title_tag.get_text(strip=True)
            link = urljoin(url, title_tag.get("href"))

            summary_tag = item.select_one("p, .summary, .excerpt, .news-excerpt, .card-text")
            summary_en = summary_tag.get_text(strip=True) if summary_tag else ""

            time_struct = None
            date_tag = item.select_one("time, .date, .post-date, .news-date, .card-date")
            if date_tag:
                date_str = date_tag.get_text(strip=True) or date_tag.get('datetime')
                if date_str:
                     time_struct = parse_date_flexible(date_str)
            
            articles.append({
                "title": title_en, "url": link, "snippet": summary_en,
                "source": "IDF News", "time_struct": time_struct,
                "needs_translation": True
            })
    except Exception as e: print(f"      爬取 IDF News 时出错: {e}")
    return articles


SCRAPER_FUNCTIONS_MAP = {
    "fetch_breakthrought1d_articles": fetch_breakthrought1d_articles,
    "fetch_dzd_articles": fetch_dzd_articles,
    "fetch_pubmed_articles": fetch_pubmed_articles,
    "fetch_idf_articles": fetch_idf_articles,
}

# --- (C) 使用讯飞星火 HTTP API 进行动态分类 ---
def categorize_article_with_llm(article_obj):
    """使用讯飞星火 HTTP API 对文章进行分类 (using Bearer Token for SPARK_LITE_HTTP_URL)"""
    global llm_call_count, MAX_LLM_CALLS
    if not SPARK_API_PASSWORD:
        print("      错误: 讯飞星火 APIPassword 未配置，无法进行 LLM 分类。将归入'综合资讯'。")
        return "综合资讯"
    if llm_call_count >= MAX_LLM_CALLS:
        return "综合资讯"

    title = article_obj.get("title", "")
    snippet = article_obj.get("snippet", "")
    text_to_classify = f"标题：{title}\n摘要：{snippet[:300]}"

    prompt = f"""请根据以下文章内容，判断它最符合下列哪个分类？请严格从列表中选择一个，并只返回分类名称，不要添加任何其他解释或文字。

可选分类列表：{', '.join(VALID_CATEGORY_NAMES)}

文章内容：
{text_to_classify}

最合适的分类名称是："""

    llm_call_count += 1

    try:
        headers = {"Content-Type": "application/json", "Authorization": f"Bearer {SPARK_API_PASSWORD}"}
        payload = {
            "model": "lite",
            "messages": [{"role": "user", "content": prompt}],
            "temperature": 0.3,
            "max_tokens": 50
        }
        response = requests.post(SPARK_LITE_HTTP_URL, headers=headers, json=payload, timeout=30)
        response.raise_for_status()
        response_data = response.json()

        llm_output = ""
        if 'choices' in response_data and response_data['choices']:
            message = response_data['choices'][0].get('message', {})
            llm_output = message.get('content', '').strip()
        elif 'header' in response_data and response_data['header'].get('code') != 0:
            print(f"      讯飞星火 API 返回错误: code={response_data['header'].get('code')}, message={response_data['header'].get('message')}")
            return "综合资讯"
        else:
            print(f"      警告: 未知的讯飞星火 API 分类响应结构: {response_data}")
            return "综合资讯"

        cleaned_output = llm_output.strip().strip('"').strip("'").replace("：","").replace(":","")
        if cleaned_output in VALID_CATEGORY_NAMES:
            return cleaned_output
        else:
            for valid_cat in VALID_CATEGORY_NAMES:
                if valid_cat in cleaned_output:
                    print(f"      警告: LLM 返回的分类 '{cleaned_output}' (原始: '{llm_output}') 包含有效分类 '{valid_cat}'. 使用 '{valid_cat}'.")
                    return valid_cat
            print(f"      警告: 讯飞星火 API 返回的分类 '{cleaned_output}' (原始: '{llm_output}') 无效或不在列表中。将归入'综合资讯'。")
            return "综合资讯"

    except requests.exceptions.RequestException as req_e:
        print(f"      调用讯飞星火 API 分类时发生网络或HTTP错误: {req_e}")
        if 'response' in locals() and response is not None:
            print(f"      响应状态码: {response.status_code}")
            try: print(f"      响应内容: {response.json()}")
            except json.JSONDecodeError: print(f"      响应内容 (非JSON): {response.text}")
        return "综合资讯"
    except Exception as e:
        print(f"      调用讯飞星火 API 分类时发生错误: {e}")
        return "综合资讯"

# --- HTML 生成逻辑 ---
def generate_html_content(all_news_data_sorted):
    current_time_str = datetime.datetime.now().strftime('%Y年%m月%d日 %H:%M:%S')
    app_timezone = os.getenv('APP_TIMEZONE', 'UTC')
    if app_timezone != 'UTC':
        try:
            tz_offset = int(os.getenv('APP_TIMEZONE_OFFSET_HOURS', '8'))
            current_time_obj = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=tz_offset)))
            current_time_str = current_time_obj.strftime('%Y年%m月%d日 %H:%M:%S %Z')
        except Exception as e:
            print(f"Error applying timezone: {e}. Falling back to server time.")
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
        .card-footer {{ margin-top: auto; }} /* Pushes footer to bottom */
        .header-main-title {{ font-size: 1.875rem; }} /* 30px */
        @media (min-width: 640px) {{ .header-main-title {{ font-size: 2.25rem; }} }} /* 36px */
        @media (min-width: 768px) {{ .header-main-title {{ font-size: 3rem; }} }} /* 48px */
        .main-container {{ padding: 1rem; }}
        @media (min-width: 768px) {{ .main-container {{ padding: 2rem; }} }}
        .loader {{ border: 5px solid #f3f3f3; border-radius: 50%; border-top: 5px solid #3498db; width: 50px; height: 50px; -webkit-animation: spin 1s linear infinite; animation: spin 1s linear infinite; margin: 20px auto; }}
        @-webkit-keyframes spin {{ 0% {{ -webkit-transform: rotate(0deg); }} 100% {{ -webkit-transform: rotate(360deg); }} }}
        @keyframes spin {{ 0% {{ transform: rotate(0deg); }} 100% {{ transform: rotate(360deg); }} }}
        .tab-buttons-container {{ display: flex; flex-wrap: wrap; margin-bottom: 1.5rem; border-bottom: 2px solid #d1d5db; }}
        .tab-button {{ padding: 0.75rem 1.25rem; margin-right: 0.5rem; margin-bottom: -2px; /* Overlap border */ border: 2px solid transparent; border-bottom: none; border-top-left-radius: 0.375rem; border-top-right-radius: 0.375rem; cursor: pointer; font-weight: 500; color: #4b5563; background-color: #f9fafb; transition: all 0.3s ease; }}
        .tab-button:hover {{ color: #1f2937; background-color: #f3f4f6; }}
        .tab-button.active {{ color: #1d4ed8; background-color: #ffffff; border-color: #d1d5db; border-bottom-color: #ffffff; /* Make bottom border of active tab white to blend with content */ }}
    </style>
</head>
<body class="bg-gray-100 text-gray-800">
    <div class="container mx-auto main-container max-w-screen-xl"> 
        <header class="text-center mb-10 md:mb-16">
            <h1 class="font-bold text-blue-700 header-main-title">糖尿病前沿资讯</h1>
            <p class="text-gray-600 mt-3 text-base md:text-lg">最近一个月动态（自动更新于：<span id="updateTime">{current_time_str}</span>）</p>
            <p class="text-sm text-gray-500 mt-2">资讯综合来源 (由 AI 智能分类)</p>
        </header>
        <div class="tab-buttons-container" id="tabButtons">"""
    first_active_category_key = None
    if any(all_news_data_sorted.values()):
        for cat_key in CATEGORIES_CONFIG.keys():
            if all_news_data_sorted.get(cat_key):
                first_active_category_key = cat_key
                break
    if not first_active_category_key:
         first_active_category_key = list(CATEGORIES_CONFIG.keys())[0] if CATEGORIES_CONFIG else None

    for category_name_key in CATEGORIES_CONFIG.keys():
        category_config = CATEGORIES_CONFIG.get(category_name_key, {})
        emoji = category_config.get("emoji", "")
        tab_id = "tab-" + html.escape(category_name_key.replace(" ", "-").replace("/", "-").lower())
        is_active = (category_name_key == first_active_category_key)
        active_class = "active" if is_active else ""
        html_output += f"""<button class="tab-button {active_class}" data-tab-target="#{tab_id}">{emoji} {html.escape(category_name_key)}</button>"""

    html_output += """</div><div id="news-content">"""
    if not any(all_news_data_sorted.values()):
        html_output += '<p class="text-center text-gray-500 text-xl py-10">抱歉，目前未能加载到最近一个月相关的糖尿病资讯。</p>'
    else:
        for category_name_key in CATEGORIES_CONFIG.keys():
            articles = all_news_data_sorted.get(category_name_key, []) 
            category_config = CATEGORIES_CONFIG.get(category_name_key, {})
            emoji = category_config.get("emoji", "")
            tab_id = "tab-" + html.escape(category_name_key.replace(" ", "-").replace("/", "-").lower())
            is_active_content = (category_name_key == first_active_category_key)
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
                    
                    category_html_content += f"""
                    <div class="news-card shadow-md hover:shadow-lg p-5">
                        <div class="card-content">
                            <h3 class="text-lg mb-2 font-semibold">
                                <a href="{url}" {'target="_blank" rel="noopener noreferrer"' if url != 'javascript:void(0);' and url.startswith('http') else ''} class="news-item-link">{title}</a>
                            </h3>
                            <p class="text-gray-600 text-sm mb-4 leading-relaxed">{snippet}</p>
                        </div>
                        <div class="card-footer text-xs text-gray-500 flex flex-wrap gap-2 items-center pt-3 border-t border-gray-200">
                            <span class="source-tag">{source_display}</span>
                            <span class="time-tag">{time_display}</span>
                        </div>
                    </div>"""
                category_html_content += "</div>" 
            category_html_content += "</div>"
            html_output += category_html_content

    html_output += f"""
        </div> </div> <footer class="text-center p-6 mt-12 text-gray-600 text-sm border-t border-gray-300"><p>&copy; {current_year} 糖尿病资讯聚合. <a href="{github_repo_url}" target="_blank" rel="noopener noreferrer" class="text-blue-600 hover:underline">项目源码</a></p><p class="mt-1">本站内容仅供参考, 不构成医疗建议。</p></footer>
    <script>
        document.addEventListener('DOMContentLoaded', function () {{
            const tabButtons = document.querySelectorAll('.tab-button');
            const tabContents = document.querySelectorAll('.tab-content');
            
            let activeTabButton = document.querySelector('.tab-button.active');
            let activeTabContent = document.querySelector('.tab-content.active');

            if (!activeTabButton && tabButtons.length > 0) {{
                tabButtons[0].classList.add('active');
                activeTabButton = tabButtons[0];
            }}
            if (activeTabButton && !activeTabContent) {{
                const targetId = activeTabButton.dataset.tabTarget;
                const targetContent = document.querySelector(targetId);
                if (targetContent) {{
                    targetContent.classList.add('active');
                }} else if (tabContents.length > 0) {{
                    tabContents[0].classList.add('active');
                }}
            }} else if (activeTabButton && activeTabContent) {{
                if (activeTabButton.dataset.tabTarget !== '#' + activeTabContent.id) {{
                    tabContents.forEach(content => content.classList.remove('active'));
                    const targetContent = document.querySelector(activeTabButton.dataset.tabTarget);
                    if (targetContent) targetContent.classList.add('active');
                }}
            }}

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
    MAX_LLM_CALLS = int(os.getenv("MAX_LLM_CALLS", "500"))
    llm_call_count = 0

    # --- 步骤一：从权威 RSS 源获取新闻 ---
    print("\n--- 正在从权威 RSS 源获取新闻 ---")
    for feed_info in AUTHORITATIVE_RSS_FEEDS:
        print(f"  处理 RSS 源: {feed_info['source_override']}")
        current_priority = feed_info.get("priority", 5) 
        needs_translation_feed = feed_info.get("needs_translation", False)
        raw_articles_from_feed = fetch_articles_from_rss(feed_info["url"], feed_info["source_override"])
        
        processed_in_feed = 0
        # skipped_due_date = 0 # 不再需要这个计数器
        added_to_candidates = 0

        for article_data in raw_articles_from_feed:
            if article_data["url"] in globally_seen_urls: continue
            
            title_to_process = article_data["title"]
            snippet_to_process = article_data["snippet"]
            time_struct = article_data["time_struct"] # 获取 time_struct

            if needs_translation_feed and llm_call_count < MAX_LLM_CALLS:
                original_title_for_norm = title_to_process
                translated_title = translate_text_with_llm(title_to_process)
                if translated_title != title_to_process:
                    title_to_process = translated_title
                
                translated_snippet = translate_text_with_llm(snippet_to_process)
                if translated_snippet != snippet_to_process:
                     snippet_to_process = translated_snippet
                
                if llm_call_count % 2 == 0:
                    time.sleep(1.1)

            # *** 移除日期过滤 ***
            # if time_struct and is_within_last_month(time_struct, today): # <- 这行被移除
            
            # 只要文章没见过，就处理并尝试加入候选池
            normalized_title_key = normalize_title(article_data["title"])
            
            time_display_str = "未知时间"
            if time_struct: # 仅在 time_struct 有效时尝试格式化
                try:
                    time_display_str = time.strftime("%Y-%m-%d", time_struct)
                except Exception:
                    pass # 保持 "未知时间"
            
            article_obj_for_storage = {
                "title": title_to_process, "url": article_data["url"], "snippet": snippet_to_process, 
                "source": article_data["source"], "time_display_str": time_display_str, 
                "time_struct": time_struct, # 存储解析后的 time_struct，用于排序
                "source_priority": current_priority, "source_type": "authoritative_rss"
            }
            
            if normalized_title_key not in unique_articles_candidates or \
               current_priority > unique_articles_candidates[normalized_title_key]["priority"]:
                unique_articles_candidates[normalized_title_key] = {
                    "article_obj": article_obj_for_storage, "priority": current_priority,
                    "url": article_data["url"]
                }
                globally_seen_urls.add(article_data["url"])
                added_to_candidates += 1
            # else: # 不需要 else 分支了，因为日期过滤已移除
            #     skipped_due_date += 1
            processed_in_feed += 1

        print(f"    来自 {feed_info['source_override']} 处理了 {processed_in_feed} 条, 新增 {added_to_candidates} 条到候选池。") # 更新日志信息
        time.sleep(1)

    # --- 步骤二：从爬虫源获取新闻 ---
    print("\n--- 正在从爬虫源获取新闻 ---")
    for scraper_info in SCRAPED_SOURCES_CONFIG:
        print(f"  处理爬虫源: {scraper_info['source_override']}")
        if scraper_info["fetch_function"] not in SCRAPER_FUNCTIONS_MAP:
            print(f"    错误: 未找到爬虫函数 {scraper_info['fetch_function']}")
            continue
            
        fetch_function = SCRAPER_FUNCTIONS_MAP[scraper_info["fetch_function"]]
        current_priority = scraper_info.get("priority", 3)
        raw_articles_from_scraper = fetch_function()
        processed_in_scraper = 0
        # skipped_due_date_scraper = 0 # 不再需要
        added_to_candidates_scraper = 0

        for article_data in raw_articles_from_scraper: 
            if article_data["url"] in globally_seen_urls: continue

            title_to_process = article_data["title"]
            snippet_to_process = article_data["snippet"]
            needs_translation_scraper = article_data.get("needs_translation", False) 
            time_struct = article_data.get("time_struct")

            if needs_translation_scraper and llm_call_count < MAX_LLM_CALLS:
                original_title_for_norm = title_to_process
                translated_title = translate_text_with_llm(title_to_process)
                if translated_title != title_to_process:
                    title_to_process = translated_title

                translated_snippet = translate_text_with_llm(snippet_to_process)
                if translated_snippet != snippet_to_process:
                    snippet_to_process = translated_snippet
                
                if llm_call_count % 2 == 0:
                     time.sleep(1.1)
            
            # *** 移除日期过滤 ***
            # if time_struct and is_within_last_month(time_struct, today): # <- 这行被移除
            
            normalized_title_key = normalize_title(article_data["title"])
            
            time_display_str = "未知时间"
            if time_struct:
                try:
                    time_display_str = time.strftime("%Y-%m-%d", time_struct)
                except Exception: pass
            
            article_obj_for_storage = {
                "title": title_to_process, "url": article_data["url"], "snippet": snippet_to_process, 
                "source": scraper_info["source_override"], "time_display_str": time_display_str, 
                "time_struct": time_struct, 
                "source_priority": current_priority, "source_type": "scraper"
            }
            if normalized_title_key not in unique_articles_candidates or \
               current_priority > unique_articles_candidates[normalized_title_key]["priority"]:
                unique_articles_candidates[normalized_title_key] = {
                    "article_obj": article_obj_for_storage, "priority": current_priority,
                    "url": article_data["url"]
                }
                globally_seen_urls.add(article_data["url"])
                added_to_candidates_scraper += 1
            # else: # 不需要 else 分支了
            #     skipped_due_date_scraper += 1
            processed_in_scraper += 1

        print(f"    来自 {scraper_info['source_override']} 处理了 {processed_in_scraper} 条, 新增 {added_to_candidates_scraper} 条到候选池。") # 更新日志信息
        time.sleep(1)

    # --- 步骤三：从 Google News RSS 获取补充新闻 ---
    print("\n--- 正在从 Google News RSS 获取补充新闻 (用于全局候选池) ---")
    google_search_term = "糖尿病 新闻 OR diabetes news" 
    print(f"  使用 Google News 搜索词: {google_search_term}")
    google_news_rss_url = f"https://news.google.com/rss/search?q={quote_plus(google_search_term)}&hl=zh-CN&gl=CN&ceid=CN:zh-Hans"
    
    raw_articles_from_google = fetch_articles_from_rss(google_news_rss_url, source_name_override=None)
    processed_in_google = 0
    # skipped_due_date_google = 0 # 不需要
    added_to_candidates_google = 0

    for article_data in raw_articles_from_google:
        if article_data["url"] in globally_seen_urls: continue
            
        time_struct = article_data["time_struct"]

        # *** 移除日期过滤 ***
        # if time_struct and is_within_last_month(time_struct, today): # <- 这行被移除
        
        normalized_title_key = normalize_title(article_data["title"])
        
        time_display_str = "未知时间"
        if time_struct:
            try:
                time_display_str = time.strftime("%Y-%m-%d", time_struct)
            except Exception: pass
        
        title_to_process_g = article_data["title"]
        snippet_to_process_g = article_data["snippet"]

        article_obj_for_storage = {
            "title": title_to_process_g, "url": article_data["url"], "snippet": snippet_to_process_g, 
            "source": article_data["source"], "time_display_str": time_display_str, 
            "time_struct": time_struct, 
            "source_priority": GOOGLE_NEWS_PRIORITY, "source_type": "google_news"
        }

        if normalized_title_key not in unique_articles_candidates or \
           GOOGLE_NEWS_PRIORITY > unique_articles_candidates[normalized_title_key]["priority"]:
            unique_articles_candidates[normalized_title_key] = {
                "article_obj": article_obj_for_storage, "priority": GOOGLE_NEWS_PRIORITY,
                "url": article_data["url"]
            }
            globally_seen_urls.add(article_data["url"])
            added_to_candidates_google += 1
        # else: # 不需要 else
        #     skipped_due_date_google += 1
        processed_in_google += 1

    print(f"    来自 Google News 处理了 {processed_in_google} 条, 新增 {added_to_candidates_google} 条到候选池 (因重复跳过 {processed_in_google - added_to_candidates_google} 条)。") # 更新日志
    time.sleep(1)

    # --- 步骤四：使用 LLM 动态分类所有候选文章 ---
    print(f"\n--- 正在对 {len(unique_articles_candidates)} 条候选文章进行动态分类 (LLM 调用上限: {MAX_LLM_CALLS}) ---")
    all_articles_by_site_category_temp = {category_name: [] for category_name in CATEGORIES_CONFIG.keys()}
    
    spark_categorization_ready = bool(SPARK_API_PASSWORD) 
    if not spark_categorization_ready:
        print("警告: 讯飞星火 APIPassword 未配置，将跳过 LLM 分类，所有文章归入'综合资讯'。")
    
    categorized_count = 0
    llm_limit_hit_printed = False

    for candidate_key, candidate_info in unique_articles_candidates.items():
        article_to_categorize = candidate_info["article_obj"]

        best_category = "综合资讯" 
        if spark_categorization_ready and llm_call_count < MAX_LLM_CALLS:
            try:
                best_category = categorize_article_with_llm(article_to_categorize)
                time.sleep(1.1) 
            except Exception as llm_e:
                print(f"    LLM 分类时发生意外错误 for '{article_to_categorize['title'][:30]}...': {llm_e}，文章将归入'综合资讯'。")
                best_category = "综合资讯"
        elif llm_call_count >= MAX_LLM_CALLS and spark_categorization_ready and not llm_limit_hit_printed:
            print(f"    已达到 LLM 调用次数上限 ({MAX_LLM_CALLS})，剩余文章将归入'综合资讯'。")
            llm_limit_hit_printed = True
            best_category = "综合资讯"
        elif not spark_categorization_ready:
             best_category = "综合资讯"
        
        if best_category not in all_articles_by_site_category_temp:
            print(f"    警告: LLM 返回的分类 '{best_category}' 不在预设分类中 for article '{article_to_categorize['title'][:30]}...'. 归入 '综合资讯'")
            best_category = "综合资讯"
            
        all_articles_by_site_category_temp[best_category].append(article_to_categorize)
        categorized_count += 1

    print(f"--- 分类完成，共处理 {categorized_count} 篇文章 ---")


    # --- 步骤五：对每个分类的文章按来源类型和日期排序并截取 ---
    print("\n--- 正在对各分类新闻进行排序和截取 ---")
    all_articles_by_site_category_final_sorted = {}
    total_final_articles = 0
    for category_name, articles_list in all_articles_by_site_category_temp.items():
        # 确保排序时优先处理有日期的文章
        articles_list.sort(key=lambda x: (
            SOURCE_TYPE_ORDER.get(x.get("source_type", "unknown"), 99),
            # 对日期进行处理：有日期的排在前面，日期越新越靠前；无日期的排在后面
            0 if x.get("time_struct") else 1, # 有日期的 key[1] 为 0, 无日期的为 1
            -(time.mktime(x["time_struct"]) if x.get("time_struct") else -float('inf')), # 日期越新，负值越小
            -x.get("source_priority", 0) # 优先级越高，负值越小
        ))
        
        all_articles_by_site_category_final_sorted[category_name] = articles_list[:MAX_ARTICLES_PER_CATEGORY]
        count_in_category = len(all_articles_by_site_category_final_sorted[category_name])
        total_final_articles += count_in_category
        print(f"  分类 '{category_name}' 排序并截取后有 {count_in_category} 条新闻。")

    print(f"--- 排序截取完成，最终展示 {total_final_articles} 篇文章 ---")


    # --- (6) 生成最终的HTML ---
    print("\n--- 正在生成最终的 HTML 文件 ---")
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

    print(f"总共调用 LLM API {llm_call_count} 次。")
    print("资讯网页生成完毕。")

