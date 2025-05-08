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
try:
    from Bio import Entrez # For PubMed API
except ImportError:
    print("错误：未找到 Biopython 库。请通过 'pip install biopython' 安装。")
    Entrez = None # Set Entrez to None if import fails

# --- (0) 从环境变量读取讯飞星火 API Keys ---
# Spark Lite HTTP Endpoint for translation uses SPARK_API_PASSWORD
SPARK_API_PASSWORD = os.getenv("SPARK_API_PASSWORD") # 读取 APIPassword (用于翻译)
SPARK_LITE_HTTP_URL = "https://spark-api-open.xf-yun.com/v1/chat/completions" # Spark Lite HTTP Endpoint (用于翻译和分类)

# Spark API credentials for categorization (using HMAC)
SPARK_APPID = os.getenv("SPARK_APPID") # 读取 AppID (用于分类)
SPARK_API_KEY = os.getenv("SPARK_API_KEY") # 读取 APIKey (用于分类)
SPARK_API_SECRET = os.getenv("SPARK_API_SECRET") # 读取 APISecret (用于分类)

# --- (1) 配置权威 RSS 源 ---
AUTHORITATIVE_RSS_FEEDS = [
    {"url": "https://www.medscape.com/cx/rss/professional.xml", "source_override": "Medscape Professional", "priority": 10, "needs_translation": True},
    {"url": "https://www.healio.com/sws/feed/news/endocrinology", "source_override": "Healio Endocrinology", "priority": 9, "needs_translation": True},
    {"url": "https://www.diabettech.com/feed/", "source_override": "Diabettech", "priority": 8, "needs_translation": True},
    # {"url": "https://thesavvydiabetic.com/feed/", "source_override": "The Savvy Diabetic", "priority": 7, "needs_translation": True}, # 403
    # {"url": "https://forum.diabetes.org.uk/boards/forums/-/index.rss", "source_override": "Diabetes UK 论坛", "priority": 6, "needs_translation": True}, # Removed by user
    # {"url": "https://www.gov.uk/government/latest.atom?organisations%5B%5D=medicines-and-healthcare-products-regulatory-agency", "source_override": "MHRA (UK)", "priority": 9, "needs_translation": True}, # Atom feed
    {"url": "https://www.fda.gov/about-fda/contact-fda/stay-informed/rss-feeds/press-releases/rss.xml", "source_override": "FDA (US) Press", "priority": 10, "needs_translation": True},
    # { "url": "YOUR_PUBMED_RSS_URL", "source_override": "PubMed (RSS Search)", "priority": 12, "needs_translation": True },
    # { "url": "YOUR_ADA_JOURNAL_RSS_URL", "source_override": "Diabetes Care (ADA)", "priority": 11, "needs_translation": True },
]

# --- (1b) 配置爬虫源 ---
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

# --- 翻译函数 (使用讯飞星火 HTTP API - Bearer Token Auth) ---
def translate_text_with_llm(text, target_lang='Chinese'):
    global llm_call_count # 声明使用全局计数器
    if not text or not isinstance(text, str) or not text.strip(): return text
    if not SPARK_API_PASSWORD: # SPARK_API_PASSWORD is used for translation
        print("      错误: 讯飞星火 APIPassword 未配置，无法进行 LLM 翻译。")
        return text
    if llm_call_count >= MAX_LLM_CALLS: # 在调用前检查上限
        print(f"      警告: 已达到 LLM 调用次数上限 ({MAX_LLM_CALLS})，跳过翻译。")
        return text # 返回原文

    max_input_length = 500 # 讯飞星火 Lite 通常限制在 8192 tokens, 简单截断以避免超长
    text_to_translate = text[:max_input_length]
    prompt = f"Please translate the following text to {target_lang}. Only return the translated text, without any introduction or explanation.\n\nText to translate:\n{text_to_translate}"
    print(f"      正在调用讯飞星火 HTTP API 翻译: {text_to_translate[:30]}...")
    llm_call_count += 1 # 增加计数器

    try:
        headers = {"Content-Type": "application/json", "Authorization": f"Bearer {SPARK_API_PASSWORD}"}
        payload = {"model": "lite", "messages": [{"role": "user", "content": prompt}], "temperature": 0.3, "max_tokens": int(len(text_to_translate) * 1.5) + 50} # 根据输入长度调整 max_tokens
        response = requests.post(SPARK_LITE_HTTP_URL, headers=headers, json=payload, timeout=20) # 使用统一的 URL
        response.raise_for_status()
        response_data = response.json()
        llm_output = ""
        if 'choices' in response_data and response_data['choices']:
            message = response_data['choices'][0].get('message', {})
            llm_output = message.get('content', '').strip()
        elif 'payload' in response_data and 'choices' in response_data['payload'] and \
             'text' in response_data['payload']['choices'] and response_data['payload']['choices']['text']: # 兼容旧版可能的响应结构
            llm_output = response_data['payload']['choices']['text'][0].get('content', '').strip()
        else:
            print(f"      警告: 未知的讯飞星火 API 翻译响应结构: {response_data}")
            return text # 返回原文
        
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
        if response.url != rss_url and "apology" in response.url: # Check for redirects to error/apology pages
            print(f"      请求被重定向到错误页面: {response.url}")
            response.raise_for_status() # Will likely raise an error if it's an apology page
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
            if not actual_source_name: # If no override, try to get from feed or title
                source_info = entry.get("source")
                actual_source_name = source_info.get("title") if source_info else "未知来源"
                if "news.google.com" in rss_url and not source_name_override: # Special handling for Google News
                    if ' - ' in title:
                        parts = title.rsplit(' - ', 1)
                        actual_source_name = parts[1] # Often the source is at the end of the title
            articles.append({
                "title": title, "url": link, "snippet": snippet,
                "source": actual_source_name, "time_struct": published_time_struct
            })
    except requests.exceptions.Timeout: print(f"      获取 RSS 源时发生超时错误: {rss_url}")
    except requests.exceptions.RequestException as e: print(f"      获取 RSS 源时发生网络错误: {e} (URL: {rss_url})")
    except Exception as e: print(f"      处理 RSS 源时发生未知错误: {e} (URL: {rss_url})")
    return articles

# --- (B) 爬虫函数定义 ---
# (所有爬虫函数与 diabetes_news_fetch_all_sources_v1 版本相同)
# ... (为简洁起见，此处省略爬虫函数定义) ...
def fetch_breakthrought1d_articles():
    BASE_URL = "https://www.breakthrought1d.org/news/"
    print(f"    正在爬取: {BASE_URL}")
    articles = []
    try:
        headers = {"User-Agent": "Mozilla/5.0"}
        response = requests.get(BASE_URL, headers=headers, timeout=20)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, "html.parser")
        for article_el in soup.select("article"): # Common selector for articles
            a_tag = article_el.select_one("h2 a") # Link usually in a heading
            if not a_tag: continue
            title_en = a_tag.get_text(strip=True)
            link = urljoin(BASE_URL, a_tag.get("href"))
            summary_tag = article_el.select_one("p") # Summary often in a paragraph
            summary_en = summary_tag.get_text(strip=True) if summary_tag else ""
            
            # Attempt to find a date - this is highly site-specific
            time_struct = None
            date_tag = article_el.select_one("time, .date, .post-date") # Common date selectors
            if date_tag:
                date_str = date_tag.get_text(strip=True) or date_tag.get('datetime')
                if date_str:
                    try:
                        # Add more date parsing formats if needed
                        dt_obj = datetime.datetime.strptime(date_str, "%B %d, %Y") # Example: August 28, 2023
                        time_struct = dt_obj.timetuple()
                    except ValueError:
                        try:
                            dt_obj = datetime.datetime.strptime(date_str, "%Y-%m-%d") # Example: 2023-08-28
                            time_struct = dt_obj.timetuple()
                        except ValueError:
                             print(f"      警告: 未能解析 Breakthrough T1D 日期字符串: {date_str} for {link}")


            if not time_struct: print(f"      警告: 未能从 {link} 提取发布日期。")
            # 返回原文，标记需要翻译
            articles.append({
                "title": title_en, "url": link, "snippet": summary_en,
                "source": "Breakthrough T1D", "time_struct": time_struct,
                "needs_translation": True # Assume English content needs translation
            })
    except Exception as e: print(f"      爬取 Breakthrough T1D 时出错: {e}")
    return articles

def fetch_myglu_articles(): # 已注释掉
    print("    跳过 MyGlu 爬虫 (已注释掉)")
    return []

def fetch_dzd_articles():
    BASE_URL = "https://www.dzd-ev.de/en/press/press-releases/press-releases-2025/index.html" # Update year if needed
    print(f"    正在爬取: {BASE_URL}")
    articles = []
    try:
        headers = {"User-Agent": "Mozilla/5.0"}
        response = requests.get(BASE_URL, headers=headers, timeout=20)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, "html.parser")

        for item in soup.select("div.teaser-text"): # Selector based on previous script
            title_tag = item.select_one("a")
            if not title_tag: continue
            title_en = title_tag.get_text(strip=True)
            link = urljoin(BASE_URL, title_tag.get("href"))
            
            summary_tag = item.select_one("p")
            summary_en = summary_tag.get_text(strip=True) if summary_tag else ""

            time_struct = None
            # DZD specific date extraction - look for a date string within the teaser or related elements
            date_span = item.select_one("span.date, div.date") # Example selectors
            if date_span:
                date_str = date_span.get_text(strip=True)
                # Example: "24.07.2023" or "July 24, 2023" - adjust strptime format accordingly
                try:
                    dt_obj = datetime.datetime.strptime(date_str, "%d.%m.%Y")
                    time_struct = dt_obj.timetuple()
                except ValueError:
                    try:
                        dt_obj = datetime.datetime.strptime(date_str, "%B %d, %Y")
                        time_struct = dt_obj.timetuple()
                    except ValueError:
                        print(f"      警告: 未能解析 DZD 日期字符串: {date_str} for {link}")
            
            if not time_struct: print(f"      警告: 未能从 {link} 提取发布日期。")
            # 返回原文，标记需要翻译
            articles.append({
                "title": title_en, "url": link, "snippet": summary_en,
                "source": "DZD News", "time_struct": time_struct,
                "needs_translation": True
            })
    except Exception as e: print(f"      爬取 DZD News 时出错: {e}")
    return articles

def fetch_adces_articles(): # 已注释掉
    print("    跳过 ADCES News 爬虫 (已注释掉)")
    return []

def fetch_panther_articles(): # 已注释掉
    print("    跳过 PANTHER Program 爬虫 (已注释掉)")
    return []

def fetch_nmpa_articles(): # 已注释掉
    print("    跳过 NMPA 爬虫 (已注释掉)")
    return []

def fetch_pubmed_articles():
    if not Entrez:
        print("    错误: Biopython (Entrez) 未加载，跳过 PubMed API 调用。")
        return []
    Entrez.email = os.getenv("PUBMED_API_EMAIL", "default_email@example.com") # Set your email
    search_term = "(diabetes[Title/Abstract]) AND (treatment[Title/Abstract] OR research[Title/Abstract] OR prevention[Title/Abstract])"
    print(f"    正在通过 PubMed API 搜索: {search_term}")
    articles = []
    MAX_PUBMED_RESULTS = 20 # Limit results to avoid excessive API calls
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
            # PubMed esummary doesn't provide a full abstract, 'Source' is often the journal name.
            # For a better snippet, you might need efetch, but it's more complex.
            # Here, we'll use 'Source' as a placeholder or combine with authors if available.
            journal_source = docsum.get("Source", "")
            authors = docsum.get("AuthorList", [])
            authors_str = ", ".join(authors[:3]) + (" et al." if len(authors) > 3 else "") if authors else "N/A"
            
            snippet_en = f"Journal: {journal_source}. Authors: {authors_str}." if journal_source or authors else "No Summary Available"

            pubdate_str = docsum.get("PubDate", "") # e.g., "2023 Aug 21" or "2023 Aug" or "2023"
            time_struct = None
            if pubdate_str:
                try:
                    dt_obj = None
                    # Try parsing different common PubMed date formats
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
            
            articles.append({
                "title": title_en, "url": link, "snippet": snippet_en,
                "source": "PubMed", "time_struct": time_struct,
                "needs_translation": True # PubMed articles are typically in English
            })
            time.sleep(0.4) # NCBI E-utilities rate limit: max 3 requests/second without API key, 10/second with.
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

        # Selectors might need adjustment based on current IDF website structure
        for item in soup.select("article.news-item, div.news-item, li.news-list-item"): # Added a common list item selector
            title_tag = item.select_one("h3 a, h2 a, .title a, .news-title a") # More robust title selection
            if not title_tag: continue
            
            title_en = title_tag.get_text(strip=True)
            link = urljoin(url, title_tag.get("href"))

            summary_tag = item.select_one("p, .summary, .excerpt, .news-excerpt") # More robust summary selection
            summary_en = summary_tag.get_text(strip=True) if summary_tag else ""

            time_struct = None
            # Try to find a date element, e.g., <time>, <span class="date">, etc.
            date_tag = item.select_one("time, .date, .post-date, .news-date") # More robust date selection
            if date_tag:
                date_str = date_tag.get_text(strip=True) or date_tag.get('datetime')
                if date_str:
                    # Example date formats: "25 August 2023", "2023-08-25"
                    try:
                        dt_obj = datetime.datetime.strptime(date_str, "%d %B %Y")
                        time_struct = dt_obj.timetuple()
                    except ValueError:
                        try:
                            dt_obj = datetime.datetime.strptime(date_str, "%Y-%m-%d")
                            time_struct = dt_obj.timetuple()
                        except ValueError:
                            print(f"      警告: 未能解析 IDF 日期字符串: {date_str} for {link}")
            
            if not time_struct: print(f"      警告: 未能从 {link} 提取发布日期。")
            
            articles.append({
                "title": title_en, "url": link, "snippet": summary_en,
                "source": "IDF News", "time_struct": time_struct,
                "needs_translation": True # Assume English content
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

# --- (C) 使用讯飞星火 HTTP API 进行动态分类 (使用 AppID, APIKey, APISecret for some models/endpoints) ---
# Note: The original categorize_article_with_llm was using HMAC auth which is typically for WebSocket.
# The translate_text_with_llm uses Bearer token (just APIPassword).
# For HTTP chat completions (like translation), Bearer token is common.
# If categorization *requires* HMAC, it implies a different endpoint or the WebSocket API.
# Assuming categorization can also use the same SPARK_LITE_HTTP_URL with Bearer token for simplicity.
# If HMAC is strictly needed for a *different* categorization endpoint, that needs to be specified.

# For now, let's assume categorization can use the same Bearer token method as translation,
# as the error was about SPARK_APPID not being defined, not an auth failure itself.
# However, the original code for categorization had a different payload structure (header, parameter, payload).
# This structure is more typical for the WebSocket API or older/different HTTP endpoints.
# The /v1/chat/completions endpoint (SPARK_LITE_HTTP_URL) uses the simpler structure seen in translate_text_with_llm.

# REVISING categorize_article_with_llm to use the simpler Bearer token auth and payload,
# consistent with SPARK_LITE_HTTP_URL and translate_text_with_llm.
# If a different Spark API (e.g., v3.5 with HMAC) is intended for categorization,
# the get_spark_authorization_headers function and the correct endpoint URL would be needed.
# The user's original code for categorization was trying to call get_spark_authorization_headers,
# but that function was not provided in the snippet.

def categorize_article_with_llm(article_obj):
    """使用讯飞星火 HTTP API 对文章进行分类 (using Bearer Token for SPARK_LITE_HTTP_URL)"""
    global llm_call_count
    # For this simplified HTTP call, only SPARK_API_PASSWORD (as Bearer token) is directly used.
    # SPARK_APPID might be implicitly part of the API_PASSWORD or not needed for this specific endpoint.
    # If SPARK_APPID, API_KEY, API_SECRET are for a *different* Spark service/auth, this function needs adjustment.
    if not SPARK_API_PASSWORD: # Using SPARK_API_PASSWORD for the Bearer token
        print("      错误: 讯飞星火 APIPassword 未配置，无法进行 LLM 分类。将归入'综合资讯'。")
        return "综合资讯"
    # The check `if not all([SPARK_APPID, SPARK_API_KEY, SPARK_API_SECRET]):` was causing the error.
    # If these are truly needed for this *specific* categorization call (even with SPARK_LITE_HTTP_URL),
    # then the API call structure (headers, payload) might need to reflect that,
    # or a different endpoint/auth method (like HMAC) is required.
    # For now, assuming SPARK_LITE_HTTP_URL with Bearer token is sufficient.

    if llm_call_count >= MAX_LLM_CALLS:
        print(f"      警告: 已达到 LLM 调用次数上限 ({MAX_LLM_CALLS})，跳过分类。")
        return "综合资讯"

    title = article_obj.get("title", "")
    snippet = article_obj.get("snippet", "")
    text_to_classify = f"标题：{title}\n摘要：{snippet[:300]}" # Limit snippet length

    # Ensure VALID_CATEGORY_NAMES is defined and accessible here
    # It's defined globally, so it should be fine.
    prompt = f"""请根据以下文章内容，判断它最符合下列哪个分类？请严格从列表中选择一个，并只返回分类名称，不要添加任何其他解释或文字。

可选分类列表：{', '.join(VALID_CATEGORY_NAMES)}

文章内容：
{text_to_classify}

最合适的分类名称是："""

    print(f"      正在调用讯飞星火 HTTP API 对 '{title[:30]}...' 进行分类...")
    llm_call_count += 1

    try:
        # Using the same Bearer token auth and payload structure as translate_text_with_llm
        headers = {"Content-Type": "application/json", "Authorization": f"Bearer {SPARK_API_PASSWORD}"}
        # The payload structure for /v1/chat/completions
        payload = {
            "model": "lite", # Or another suitable model like "general" if available for this task
            "messages": [{"role": "user", "content": prompt}],
            "temperature": 0.3, # Lower temperature for more deterministic classification
            "max_tokens": 50 # Category name should be short
        }
        
        response = requests.post(SPARK_LITE_HTTP_URL, headers=headers, json=payload, timeout=30)
        response.raise_for_status()
        response_data = response.json()

        llm_output = ""
        if 'choices' in response_data and response_data['choices']:
            message = response_data['choices'][0].get('message', {})
            llm_output = message.get('content', '').strip()
        # Check for error in response header (some Spark APIs include error details here)
        elif 'header' in response_data and response_data['header'].get('code') != 0:
            print(f"      讯飞星火 API 返回错误: code={response_data['header'].get('code')}, message={response_data['header'].get('message')}")
            return "综合资讯"
        else:
            print(f"      警告: 未知的讯飞星火 API 分类响应结构: {response_data}")
            return "综合资讯"

        print(f"      讯飞星火 API 返回: '{llm_output}'")

        cleaned_output = llm_output.strip().strip('"').strip("'").replace("：","").replace(":","")
        if cleaned_output in VALID_CATEGORY_NAMES:
            print(f"      文章 '{title[:30]}...' 成功分类到 '{cleaned_output}'")
            return cleaned_output
        else:
            # Attempt to find a partial match if the LLM adds extra text
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
        # This is where the original "name 'SPARK_APPID' is not defined" would be caught if the old check was present
        print(f"      调用讯飞星火 API 分类时发生错误: {e}")
        return "综合资讯"


# --- HTML 生成逻辑 ---
def generate_html_content(all_news_data_sorted):
    # (此函数内容与 diabetes_news_fetch_tabs_v1 中的 generate_html_content 完全相同)
    # ... (省略 HTML 生成代码) ...
    current_time_str = datetime.datetime.now().strftime('%Y年%m月%d日 %H:%M:%S')
    app_timezone = os.getenv('APP_TIMEZONE', 'UTC') # Default to UTC if not set
    if app_timezone != 'UTC':
        try:
            # Attempt to use a fixed offset, e.g., UTC+8 for Beijing time
            # For more robust timezone handling, consider pytz library
            tz_offset = int(os.getenv('APP_TIMEZONE_OFFSET_HOURS', '8')) # Default to +8 if APP_TIMEZONE is not UTC
            current_time_obj = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=tz_offset)))
            current_time_str = current_time_obj.strftime('%Y年%m月%d日 %H:%M:%S %Z')
        except Exception as e:
            print(f"Error applying timezone: {e}. Falling back to server time.")
            current_time_str = datetime.datetime.now().strftime('%Y年%m月%d日 %H:%M:%S (服务器时间)')


    current_year = datetime.datetime.now().year
    github_repo_url = f"https://github.com/{os.getenv('GITHUB_REPOSITORY', 'doudou-ux/diabetes-news')}" # Default repo
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
    <div class="container mx-auto main-container max-w-screen-xl"> {{/* Added max-w-screen-xl for better layout on large screens */}}
        <header class="text-center mb-10 md:mb-16">
            <h1 class="font-bold text-blue-700 header-main-title">糖尿病前沿资讯</h1>
            <p class="text-gray-600 mt-3 text-base md:text-lg">最近一个月动态（自动更新于：<span id="updateTime">{current_time_str}</span>）</p>
            <p class="text-sm text-gray-500 mt-2">资讯综合来源 (由 AI 智能分类)</p>
        </header>
        <div class="tab-buttons-container" id="tabButtons">"""
    # Determine the first tab that should be active
    first_active_category_key = None
    if any(all_news_data_sorted.values()): # If there's any news
        for cat_key in CATEGORIES_CONFIG.keys():
            if all_news_data_sorted.get(cat_key):
                first_active_category_key = cat_key
                break
    if not first_active_category_key: # Fallback to the very first category if all are empty or no news
         first_active_category_key = list(CATEGORIES_CONFIG.keys())[0] if CATEGORIES_CONFIG else None

    for category_name_key in CATEGORIES_CONFIG.keys():
        category_config = CATEGORIES_CONFIG.get(category_name_key, {})
        emoji = category_config.get("emoji", "")
        tab_id = "tab-" + html.escape(category_name_key.replace(" ", "-").replace("/", "-").lower())
        is_active = (category_name_key == first_active_category_key)
        active_class = "active" if is_active else ""
        html_output += f"""<button class="tab-button {active_class}" data-tab-target="#{tab_id}">{emoji} {html.escape(category_name_key)}</button>"""

    html_output += """</div><div id="news-content">"""
    if not any(all_news_data_sorted.values()): # Check if any category has news
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
                    snippet = html.escape(snippet_raw[:150] + ('...' if len(snippet_raw) > 150 else '')) # Truncate snippet
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

    # --- 注意：f-string 中的 JavaScript 大括号需要转义 ---
    # JavaScript for tab functionality
    # Ensure that the first tab is correctly activated if no specific active class is set by Python
    html_output += f"""
        </div> </div> <footer class="text-center p-6 mt-12 text-gray-600 text-sm border-t border-gray-300"><p>&copy; {current_year} 糖尿病资讯聚合. <a href="{github_repo_url}" target="_blank" rel="noopener noreferrer" class="text-blue-600 hover:underline">项目源码</a></p><p class="mt-1">本站内容仅供参考, 不构成医疗建议。</p></footer>
    <script>
        document.addEventListener('DOMContentLoaded', function () {{ // JS 大括号需要转义: {{ becomes {{ and }} becomes }}
            const tabButtons = document.querySelectorAll('.tab-button');
            const tabContents = document.querySelectorAll('.tab-content');
            
            // Ensure consistent active state if Python logic didn't set one clearly
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
                     // Fallback if target not found, activate first content
                    tabContents[0].classList.add('active');
                }}
            }} else if (activeTabButton && activeTabContent) {{
                // If both are set, ensure they correspond
                if (activeTabButton.dataset.tabTarget !== '#' + activeTabContent.id) {{
                    tabContents.forEach(content => content.classList.remove('active')); // Deactivate all
                    const targetContent = document.querySelector(activeTabButton.dataset.tabTarget);
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
        }}); // JS 大括号需要转义
    </script>
</body></html>"""
    return html_output

# --- (5) 主执行逻辑 ---
if __name__ == "__main__":
    print("开始从多个 RSS 源和爬虫生成糖尿病资讯网页...")
    
    unique_articles_candidates = {} # Stores unique articles by normalized title, preferring higher priority
    globally_seen_urls = set() # To avoid processing the exact same URL multiple times if it appears in different feeds
    today = datetime.date.today()
    MAX_ARTICLES_PER_CATEGORY = 10 # Max articles to show per category in the final HTML
    MAX_LLM_CALLS = int(os.getenv("MAX_LLM_CALLS", "500")) # Max LLM calls for translation/categorization, configurable via env var, default to 500
    llm_call_count = 0

    # --- 步骤一：从权威 RSS 源获取新闻 ---
    print("\n--- 正在从权威 RSS 源获取新闻 ---")
    for feed_info in AUTHORITATIVE_RSS_FEEDS:
        print(f"  处理 RSS 源: {feed_info['source_override']}")
        current_priority = feed_info.get("priority", 5) 
        needs_translation_feed = feed_info.get("needs_translation", False)
        raw_articles_from_feed = fetch_articles_from_rss(feed_info["url"], feed_info["source_override"])
        
        processed_in_feed = 0
        for article_data in raw_articles_from_feed:
            if article_data["url"] in globally_seen_urls: continue # Skip if URL already processed from another source
            
            title_to_process = article_data["title"]
            snippet_to_process = article_data["snippet"]
            
            if needs_translation_feed and llm_call_count < MAX_LLM_CALLS:
                original_title_for_norm = title_to_process # Use original for normalization key if translation fails
                translated_title = translate_text_with_llm(title_to_process)
                if translated_title != title_to_process: # If translation happened and is different
                    title_to_process = translated_title
                else: # Translation failed or returned original
                    print(f"      翻译标题失败或返回原文: {title_to_process[:30]}...")

                translated_snippet = translate_text_with_llm(snippet_to_process)
                if translated_snippet != snippet_to_process:
                     snippet_to_process = translated_snippet
                else:
                    print(f"      翻译摘要失败或返回原文: {snippet_to_process[:30]}...")
                time.sleep(1.1) # API call delay

            if is_within_last_month_rss(article_data["time_struct"], today):
                normalized_title_key = normalize_title(article_data["title"]) # Normalize original title for uniqueness key
                
                time_display_str = "未知时间"
                if article_data["time_struct"]:
                    try: time_display_str = time.strftime("%Y-%m-%d", article_data["time_struct"])
                    except: pass # Keep "未知时间" if formatting fails
                
                article_obj_for_storage = {
                    "title": title_to_process, "url": article_data["url"], "snippet": snippet_to_process, 
                    "source": article_data["source"], "time_display_str": time_display_str, 
                    "time_struct": article_data["time_struct"], 
                    "source_priority": current_priority, "source_type": "authoritative_rss"
                }
                
                # Add or replace if this article is higher priority
                if normalized_title_key not in unique_articles_candidates or \
                   current_priority > unique_articles_candidates[normalized_title_key]["priority"]:
                    unique_articles_candidates[normalized_title_key] = {
                        "article_obj": article_obj_for_storage, "priority": current_priority,
                        "url": article_data["url"] # Keep URL for easy access
                    }
                    globally_seen_urls.add(article_data["url"]) # Mark URL as seen
                    processed_in_feed +=1
        print(f"    来自 {feed_info['source_override']} 的 {processed_in_feed} 条新文章加入候选池。")
        time.sleep(1) # Delay between different RSS feeds

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

        for article_data in raw_articles_from_scraper: 
            if article_data["url"] in globally_seen_urls: continue

            title_to_process = article_data["title"]
            snippet_to_process = article_data["snippet"]
            needs_translation_scraper = article_data.get("needs_translation", False) 
            
            if needs_translation_scraper and llm_call_count < MAX_LLM_CALLS:
                original_title_for_norm = title_to_process
                translated_title = translate_text_with_llm(title_to_process)
                if translated_title != title_to_process:
                    title_to_process = translated_title
                else:
                    print(f"      翻译标题失败或返回原文: {title_to_process[:30]}...")

                translated_snippet = translate_text_with_llm(snippet_to_process)
                if translated_snippet != snippet_to_process:
                    snippet_to_process = translated_snippet
                else:
                    print(f"      翻译摘要失败或返回原文: {snippet_to_process[:30]}...")
                time.sleep(1.1)
            
            # Date filtering for scrapers (assuming they also provide 'time_struct')
            if is_within_last_month_rss(article_data.get("time_struct"), today): # Use .get for safety
                normalized_title_key = normalize_title(article_data["title"]) # Original title for key
                
                time_display_str = "未知时间"
                if article_data.get("time_struct"):
                    try: time_display_str = time.strftime("%Y-%m-%d", article_data["time_struct"])
                    except: pass
                
                article_obj_for_storage = {
                    "title": title_to_process, "url": article_data["url"], "snippet": snippet_to_process, 
                    "source": scraper_info["source_override"], "time_display_str": time_display_str, 
                    "time_struct": article_data.get("time_struct"), 
                    "source_priority": current_priority, "source_type": "scraper"
                }
                if normalized_title_key not in unique_articles_candidates or \
                   current_priority > unique_articles_candidates[normalized_title_key]["priority"]:
                    unique_articles_candidates[normalized_title_key] = {
                        "article_obj": article_obj_for_storage, "priority": current_priority,
                        "url": article_data["url"]
                    }
                    globally_seen_urls.add(article_data["url"])
                    processed_in_scraper +=1
        print(f"    来自 {scraper_info['source_override']} 的 {processed_in_scraper} 条新文章加入候选池。")
        time.sleep(1)

    # --- 步骤三：从 Google News RSS 获取补充新闻 ---
    # This step can be merged or re-evaluated. If Google News is just another RSS,
    # it could be part of AUTHORITATIVE_RSS_FEEDS with its own priority.
    # The current logic adds it to the same candidate pool.
    print("\n--- 正在从 Google News RSS 获取补充新闻 (用于全局候选池) ---")
    google_search_term = "糖尿病 新闻 OR diabetes news" 
    print(f"  使用 Google News 搜索词: {google_search_term}")
    # Ensure URL encoding for the search term if it contains special characters, though html.escape might be for HTML context.
    # urllib.parse.quote_plus is better for URL parameters.
    from urllib.parse import quote_plus
    google_news_rss_url = f"https://news.google.com/rss/search?q={quote_plus(google_search_term)}&hl=zh-CN&gl=CN&ceid=CN:zh-Hans"
    
    # Google News items are often in Chinese if hl=zh-CN, but source might be English.
    # For simplicity, let's assume Google News items might need translation based on a flag or content check later.
    # Or, assume they are mostly in the target language (Chinese) due to hl=zh-CN.
    # The original code did not explicitly translate Google News items here.
    
    raw_articles_from_google = fetch_articles_from_rss(google_news_rss_url, source_name_override=None) # Let fetch_articles_from_rss determine source
    processed_in_google = 0
    for article_data in raw_articles_from_google:
        # Check if URL already exists from a higher-priority source
        if article_data["url"] in globally_seen_urls:
            # If it's already there, we only care if this Google News entry is somehow "better" (e.g. more complete)
            # but the current priority system would handle this if Google News has a low priority.
            # For now, if URL is seen, skip, as higher priority sources would have added it.
            continue
            
        if is_within_last_month_rss(article_data["time_struct"], today):
            normalized_title_key = normalize_title(article_data["title"])
            
            time_display_str = "未知时间"
            if article_data["time_struct"]:
                try: time_display_str = time.strftime("%Y-%m-%d", article_data["time_struct"])
                except: pass
            
            # Decide if Google News items need translation.
            # For now, assume they don't, or handle it like other feeds if a 'needs_translation' flag was set for Google News.
            # The original code did not translate them at this stage.
            title_to_process_g = article_data["title"]
            snippet_to_process_g = article_data["snippet"]
            # Example: if you wanted to translate all Google News items:
            # if llm_call_count < MAX_LLM_CALLS:
            #    title_to_process_g = translate_text_with_llm(article_data["title"])
            #    snippet_to_process_g = translate_text_with_llm(article_data["snippet"])
            #    time.sleep(1.1)


            article_obj_for_storage = {
                "title": title_to_process_g, "url": article_data["url"], "snippet": snippet_to_process_g, 
                "source": article_data["source"], "time_display_str": time_display_str, 
                "time_struct": article_data["time_struct"], 
                "source_priority": GOOGLE_NEWS_PRIORITY, "source_type": "google_news"
            }
            # Add or replace only if this Google News item is higher priority than what's there for the same title
            # (which is unlikely if GOOGLE_NEWS_PRIORITY is low, unless the title wasn't seen before)
            if normalized_title_key not in unique_articles_candidates or \
               GOOGLE_NEWS_PRIORITY > unique_articles_candidates[normalized_title_key]["priority"]:
                unique_articles_candidates[normalized_title_key] = {
                    "article_obj": article_obj_for_storage, "priority": GOOGLE_NEWS_PRIORITY,
                    "url": article_data["url"]
                }
                globally_seen_urls.add(article_data["url"]) # Mark as seen if added
                processed_in_google +=1
    print(f"    来自 Google News 的 {processed_in_google} 条新文章加入候选池。")
    time.sleep(1)

    # --- 步骤四：使用 LLM 动态分类所有候选文章 ---
    print(f"\n--- 正在对 {len(unique_articles_candidates)} 条候选文章进行动态分类 (LLM 调用上限: {MAX_LLM_CALLS}) ---")
    all_articles_by_site_category_temp = {category_name: [] for category_name in CATEGORIES_CONFIG.keys()}
    
    # spark_api_ready for categorization depends on SPARK_API_PASSWORD for the revised categorize_article_with_llm
    spark_categorization_ready = bool(SPARK_API_PASSWORD) 
    if not spark_categorization_ready:
        print("警告: 讯飞星火 APIPassword 未配置，将跳过 LLM 分类，所有文章归入'综合资讯'。")
    
    categorized_count = 0
    for candidate_key, candidate_info in unique_articles_candidates.items(): # Iterate through unique articles
        article_to_categorize = candidate_info["article_obj"]
        # URL has already been checked for global uniqueness when adding to candidates.

        best_category = "综合资讯" 
        if spark_categorization_ready and llm_call_count < MAX_LLM_CALLS:
            try:
                best_category = categorize_article_with_llm(article_to_categorize)
                # llm_call_count is incremented inside categorize_article_with_llm
                time.sleep(1.1) 
            except Exception as llm_e:
                print(f"    LLM 分类时发生意外错误 for '{article_to_categorize['title'][:30]}...': {llm_e}，文章将归入'综合资讯'。")
                best_category = "综合资讯" # Fallback category
        elif llm_call_count >= MAX_LLM_CALLS and spark_categorization_ready : # Only print this if API was ready but limit reached
            print(f"    已达到 LLM 调用次数上限 ({MAX_LLM_CALLS})，文章 '{article_to_categorize['title'][:30]}...' 将归入'综合资讯'。")
            best_category = "综合资讯"
        # If not spark_categorization_ready, it defaults to "综合资讯" anyway.
        
        # Ensure the category exists in our config
        if best_category not in all_articles_by_site_category_temp:
            print(f"    警告: LLM 返回的分类 '{best_category}' 不在预设分类中 for article '{article_to_categorize['title'][:30]}...'. 归入 '综合资讯'")
            best_category = "综合资讯"
            
        all_articles_by_site_category_temp[best_category].append(article_to_categorize)
        categorized_count += 1
        if categorized_count % 10 == 0: # Print progress every 10 articles
            print(f"    已分类 {categorized_count}/{len(unique_articles_candidates)} 文章...")


    # --- 步骤五：对每个分类的文章按来源类型和日期排序并截取 ---
    print("\n--- 正在对各分类新闻进行排序和截取 ---")
    all_articles_by_site_category_final_sorted = {}
    for category_name, articles_list in all_articles_by_site_category_temp.items():
        # Sort by: 1. Source Type (authoritative > scraper > google_news > unknown)
        #          2. Date (most recent first)
        #          3. Source Priority (higher number is better, so negate for ascending sort or use reverse)
        articles_list.sort(key=lambda x: (
            SOURCE_TYPE_ORDER.get(x.get("source_type", "unknown"), 99), # Lower is better
            -(time.mktime(x["time_struct"]) if x.get("time_struct") else -float('inf')), # Negative for descending date
            -x.get("source_priority", 0) # Negative for descending priority
        ))
        all_articles_by_site_category_final_sorted[category_name] = articles_list[:MAX_ARTICLES_PER_CATEGORY]
        print(f"  分类 '{category_name}' 排序并截取后有 {len(all_articles_by_site_category_final_sorted[category_name])} 条新闻。")
        if all_articles_by_site_category_final_sorted[category_name]:
            # Log first few items' source types and dates for verification
            # print(f"    排序后前几条: ")
            # for art_item in all_articles_by_site_category_final_sorted[category_name][:3]:
            #     print(f"      - [{art_item.get('source_type')}, {art_item.get('time_display_str')}, Prio:{art_item.get('source_priority')}] {art_item.get('title')[:30]}...")
            pass


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
