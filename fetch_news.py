# fetch_news.py
import datetime
import html
import os
import re # 用于规范化标题
import time
import requests
import feedparser
from bs4 import BeautifulSoup
from urllib.parse import urljoin
from deep_translator import GoogleTranslator
try:
    from Bio import Entrez
except ImportError:
    print("错误：未找到 Biopython 库。请通过 'pip install biopython' 安装。")
    Entrez = None

# --- (1) 配置权威 RSS 源 ---
# 添加 priority 和 needs_translation
AUTHORITATIVE_RSS_FEEDS = [
    {"url": "https://www.medscape.com/rss/public/diabetes.xml", "source_override": "Medscape Diabetes", "target_categories": ["最新研究", "治疗进展"], "priority": 10, "needs_translation": True},
    {"url": "https://www.healio.com/news/endocrinology/rss", "source_override": "Healio Endocrinology", "target_categories": ["最新研究", "治疗进展"], "priority": 9, "needs_translation": True},
    {"url": "https://www.diabettech.com/feed/", "source_override": "Diabettech", "target_categories": ["治疗进展", "最新研究"], "priority": 8, "needs_translation": True},
    {"url": "https://thesavvydiabetic.com/feed/", "source_override": "The Savvy Diabetic", "target_categories": ["患者故事与心理支持", "预防与生活方式"], "priority": 7, "needs_translation": True},
    {"url": "https://forum.diabetes.org.uk/boards/forums/-/index.rss", "source_override": "Diabetes UK 论坛", "target_categories": ["患者故事与心理支持"], "priority": 6, "needs_translation": True},
    {"url": "https://www.gov.uk/government/organisations/medicines-and-healthcare-products-regulatory-agency.rss", "source_override": "MHRA (UK)", "target_categories": ["政策/医保信息", "治疗进展"], "priority": 9, "needs_translation": True},
    {"url": "https://www.fda.gov/rss.xml", "source_override": "FDA (US)", "target_categories": ["政策/医保信息", "治疗进展", "最新研究"], "priority": 10, "needs_translation": True},
    # { "url": "YOUR_PUBMED_RSS_URL", "source_override": "PubMed (RSS Search)", "target_categories": ["最新研究"], "priority": 12, "needs_translation": True },
    # { "url": "YOUR_ADA_JOURNAL_RSS_URL", "source_override": "Diabetes Care (ADA)", "target_categories": ["最新研究", "治疗进展"], "priority": 11, "needs_translation": True },
]

# --- (1b) 配置爬虫源 ---
SCRAPED_SOURCES_CONFIG = [
    {"name": "Breakthrough T1D News", "fetch_function": "fetch_breakthrought1d_articles", "source_override": "Breakthrough T1D", "target_categories": ["最新研究", "治疗进展"], "priority": 8},
    {"name": "MyGlu Articles", "fetch_function": "fetch_myglu_articles", "source_override": "MyGlu", "target_categories": ["患者故事与心理支持", "预防与生活方式"], "priority": 7},
    {"name": "DZD News (2025)", "fetch_function": "fetch_dzd_articles", "source_override": "DZD News", "target_categories": ["最新研究"], "priority": 9},
    {"name": "ADCES News", "fetch_function": "fetch_adces_articles", "source_override": "ADCES News", "target_categories": ["治疗进展", "预防与生活方式", "政策/医保信息"], "priority": 8},
    {"name": "PANTHER Program News", "fetch_function": "fetch_panther_articles", "source_override": "PANTHER Program", "target_categories": ["最新研究", "治疗进展"], "priority": 7},
    {"name": "NMPA Policies", "fetch_function": "fetch_nmpa_articles", "source_override": "NMPA", "target_categories": ["政策/医保信息", "治疗进展"], "priority": 10},
    {"name": "PubMed API Search", "fetch_function": "fetch_pubmed_articles", "source_override": "PubMed", "target_categories": ["最新研究"], "priority": 12},
    {"name": "IDF News", "fetch_function": "fetch_idf_articles", "source_override": "IDF News", "target_categories": ["最新研究", "预防与生活方式", "政策/医保信息"], "priority": 9},
]

GOOGLE_NEWS_PRIORITY = 1
SOURCE_TYPE_ORDER = { # 定义来源类型的排序顺序，值越小越靠前
    'authoritative_rss': 0,
    'scraper': 1,
    'google_news': 2,
    'unknown': 99 # 未知类型排最后
}

# --- (2) 配置网站展示的分类及对应的 Google News 补充关键词 ---
CATEGORIES_CONFIG = {
    "最新研究": {"keywords": "糖尿病 最新论文 OR 糖尿病技术突破 OR 糖尿病机制研究 OR 医学会议糖尿病 OR GLP-1糖尿病 OR SGLT2糖尿病 OR 胰岛β细胞 OR 胰岛素敏感性", "emoji": "🔬"},
    "治疗进展": {"keywords": "糖尿病新药 OR 糖尿病适应症扩展 OR 糖尿病设备研发 OR AI辅助诊疗糖尿病 OR 达格列净 OR 司美格鲁肽 OR CGM OR 连续血糖监测 OR 胰岛素泵", "emoji": "💊"},
    "饮食与营养": {"keywords": "糖尿病饮食指南 OR 碳水交换表 OR 糖尿病食谱 OR 低GI饮食糖尿病 OR 高蛋白饮食糖尿病 OR 间歇性断食糖尿病 OR 膳食纤维糖尿病", "emoji": "🥗"},
    "预防与生活方式": {"keywords": "糖尿病运动建议 OR 糖尿病睡眠 OR 糖尿病减重 OR 控糖计划 OR 糖尿病早期筛查 OR 糖耐量异常 OR 体脂管理糖尿病 OR 糖尿病步数目标", "emoji": "🏃‍♀️"},
    "并发症管理": {"keywords": "糖尿病足 OR 糖尿病视网膜病变 OR 糖尿病肾病 OR 糖尿病神经病变 OR 糖网病 OR 微血管病变糖尿病 OR 尿白蛋白糖尿病", "emoji": "🩺"},
    "患者故事与心理支持": {"keywords": "糖尿病控糖经验 OR 糖尿病心理支持 OR 糖尿病家庭支持 OR 糖尿病患者故事 OR 糖尿病医生问答", "emoji": "😊"},
    "政策/医保信息": {"keywords": "糖尿病药品纳保 OR 糖尿病医保报销 OR 糖尿病社区慢病随访 OR 国家药监局糖尿病政策 OR 医保局糖尿病政策", "emoji": "📄"}
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
    if not time_struct: return False
    try:
        article_date = datetime.date(time_struct.tm_year, time_struct.tm_mon, time_struct.tm_mday)
        thirty_days_ago = today_date_obj - datetime.timedelta(days=30)
        return thirty_days_ago <= article_date <= today_date_obj
    except Exception as e:
        # print(f"      [is_within_last_month_rss] 日期转换错误: {e} - Time Struct: {time_struct}") # 日志可能过多，暂时注释
        return False

# --- 帮助函数：清理 HTML ---
def clean_html(raw_html):
    if not raw_html: return ""
    try: return BeautifulSoup(raw_html, "html.parser").get_text(separator=' ', strip=True)
    except Exception: return raw_html

# --- 翻译函数 ---
def translate_text(text, target_lang='zh-CN'):
    if not text: return ""
    try:
        translated_text = GoogleTranslator(source='auto', target=target_lang).translate(text=text, timeout=10)
        if not translated_text or translated_text == text:
             # print(f"      翻译可能未成功或无需翻译: 原文: {text[:50]}...") # 日志可能过多，暂时注释
             return text
        print(f"      翻译成功: {text[:30]}... -> {translated_text[:30]}...")
        return translated_text
    except Exception as e:
        print(f"      翻译错误: {e} - 原文: {text[:50]}...")
        return text

# --- (A) RSS 源获取函数 ---
def fetch_articles_from_rss(rss_url, source_name_override=None):
    # (与 diabetes_news_fetch_all_sources_v1 版本相同)
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
# (所有爬虫函数 fetch_breakthrought1d_articles, fetch_myglu_articles, fetch_dzd_articles,
# fetch_adces_articles, fetch_panther_articles, fetch_nmpa_articles,
# fetch_pubmed_articles, fetch_idf_articles 与 diabetes_news_fetch_all_sources_v1 版本相同)
# ... (为简洁起见，此处省略爬虫函数定义，假设它们已正确定义) ...
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
            time_struct = None # 日期提取逻辑缺失
            if not time_struct: print(f"      警告: 未能从 {link} 提取发布日期。")
            title_cn = translate_text(title_en)
            summary_cn = translate_text(summary_en)
            articles.append({
                "title": title_cn, "url": link, "snippet": summary_cn,
                "source": "Breakthrough T1D", "time_struct": time_struct
            })
    except Exception as e: print(f"      爬取 Breakthrough T1D 时出错: {e}")
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
        for item in soup.select("div.node-article"):
            a_tag = item.select_one("h2 a")
            if not a_tag: continue
            title_en = a_tag.get_text(strip=True)
            link = urljoin(BASE_URL, a_tag.get("href"))
            summary_tag = item.select_one("div.field-name-body div.field-item")
            summary_en = summary_tag.get_text(strip=True) if summary_tag else ""
            time_struct = None # 日期提取逻辑缺失
            if not time_struct: print(f"      警告: 未能从 {link} 提取发布日期。")
            title_cn = translate_text(title_en)
            summary_cn = translate_text(summary_en)
            articles.append({
                "title": title_cn, "url": link, "snippet": summary_cn,
                "source": "MyGlu", "time_struct": time_struct
            })
    except Exception as e: print(f"      爬取 MyGlu 时出错: {e}")
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
        for item in soup.select("div.teaser-text"):
            title_tag = item.select_one("a")
            if not title_tag: continue
            title_en = title_tag.get_text(strip=True)
            link = urljoin(BASE_URL, title_tag.get("href"))
            summary_tag = item.select_one("p")
            summary_en = summary_tag.get_text(strip=True) if summary_tag else ""
            time_struct = None # 日期提取逻辑缺失
            if not time_struct: print(f"      警告: 未能从 {link} 提取发布日期。")
            title_cn = translate_text(title_en)
            summary_cn = translate_text(summary_en)
            articles.append({
                "title": title_cn, "url": link, "snippet": summary_cn,
                "source": "DZD News", "time_struct": time_struct
            })
    except Exception as e: print(f"      爬取 DZD News 时出错: {e}")
    return articles

def fetch_adces_articles():
    BASE_URL = "https://www.adces.org/news"
    print(f"    正在爬取: {BASE_URL}")
    articles = []
    try:
        headers = {"User-Agent": "Mozilla/5.0"}
        response = requests.get(BASE_URL, headers=headers, timeout=20)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, "html.parser")
        for item in soup.select("div.blog-listing > div.card"):
            title_tag = item.select_one("h4.card-title a")
            if not title_tag: continue
            title_en = title_tag.get_text(strip=True)
            link = urljoin(BASE_URL, title_tag.get("href"))
            summary_tag = item.select_one("div.card-body p")
            summary_en = summary_tag.get_text(strip=True) if summary_tag else ""
            time_struct = None # 日期提取逻辑缺失
            if not time_struct: print(f"      警告: 未能从 {link} 提取发布日期。")
            title_cn = translate_text(title_en)
            summary_cn = translate_text(summary_en)
            articles.append({
                "title": title_cn, "url": link, "snippet": summary_cn,
                "source": "ADCES News", "time_struct": time_struct
            })
    except Exception as e: print(f"      爬取 ADCES News 时出错: {e}")
    return articles

def fetch_panther_articles():
    BASE_URL = "https://www.pantherprogram.org/news-events"
    print(f"    正在爬取: {BASE_URL}")
    articles = []
    try:
        headers = {"User-Agent": "Mozilla/5.0"}
        response = requests.get(BASE_URL, headers=headers, timeout=20)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, "html.parser")
        for item in soup.select("div.view-content > div.views-row"):
            title_tag = item.select_one("h3 a")
            if not title_tag: continue
            title_en = title_tag.get_text(strip=True)
            link = urljoin(BASE_URL, title_tag.get("href"))
            summary_tag = item.select_one("div.field-content p")
            summary_en = summary_tag.get_text(strip=True) if summary_tag else ""
            time_struct = None # 日期提取逻辑缺失
            if not time_struct: print(f"      警告: 未能从 {link} 提取发布日期。")
            title_cn = translate_text(title_en)
            summary_cn = translate_text(summary_en)
            articles.append({
                "title": title_cn, "url": link, "snippet": summary_cn,
                "source": "PANTHER Program", "time_struct": time_struct
            })
    except Exception as e: print(f"      爬取 PANTHER Program 时出错: {e}")
    return articles

def fetch_nmpa_articles():
    url = "https://www.nmpa.gov.cn/yaowu/zhengce/zhengcefabu.html"
    print(f"    正在爬取: {url}")
    articles = []
    try:
        headers = {"User-Agent": "Mozilla/5.0"}
        response = requests.get(url, headers=headers, timeout=20)
        response.raise_for_status()
        response.encoding = 'utf-8'
        soup = BeautifulSoup(response.text, 'html.parser')
        news_list_container = soup.select_one(".list")
        if not news_list_container:
             print(f"      错误: 未能在 NMPA 页面找到新闻列表容器。")
             return []
        for item in news_list_container.select("li"):
            a_tag = item.select_one("a")
            if not a_tag or not a_tag.get("href"): continue
            title = a_tag.get_text(strip=True)
            link = urljoin(url, a_tag.get("href"))
            snippet = title
            date_tag = item.select_one("span")
            date_str = date_tag.get_text(strip=True) if date_tag else None
            time_struct = None
            if date_str:
                try:
                    dt_obj = datetime.datetime.strptime(date_str, "%Y-%m-%d")
                    time_struct = dt_obj.timetuple()
                    print(f"      成功解析 NMPA 日期: {date_str}")
                except ValueError: print(f"      警告: 未能解析 NMPA 日期格式: {date_str}")
            else: print(f"      警告: 未能从 {link} 提取发布日期。")
            articles.append({
                "title": title, "url": link, "snippet": snippet,
                "source": "NMPA", "time_struct": time_struct
            })
    except Exception as e: print(f"      爬取 NMPA 时出错: {e}")
    return articles

def fetch_pubmed_articles():
    if not Entrez:
        print("    错误: Biopython (Entrez) 未加载，跳过 PubMed API 获取。")
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
                    if dt_obj:
                        time_struct = dt_obj.timetuple()
                        # print(f"      成功解析 PubMed 日期: {pubdate_str}") # 日志可能过多
                    # else: print(f"      警告: 未能解析 PubMed 日期格式: {pubdate_str}")
                except Exception as date_e: print(f"      解析 PubMed 日期时出错: {date_e} - {pubdate_str}")
            # else: print(f"      警告: 未能从 PubMed 获取发布日期 (PMID: {pmid})。")
            title_cn = translate_text(title_en)
            snippet_cn = translate_text(snippet_en)
            articles.append({
                "title": title_cn, "url": link, "snippet": f"{snippet_cn} (作者: {', '.join(authors)[:50]}...)",
                "source": "PubMed", "time_struct": time_struct
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
            time_struct = None # TODO: Add date extraction logic here
            if not time_struct: print(f"      警告: 未能从 {link} 提取发布日期。")
            title_cn = translate_text(title_en)
            summary_cn = translate_text(summary_en)
            articles.append({
                "title": title_cn, "url": link, "snippet": summary_cn,
                "source": "IDF News", "time_struct": time_struct
            })
    except Exception as e: print(f"      爬取 IDF News 时出错: {e}")
    return articles

# 更新 SCRAPER_FUNCTIONS_MAP
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

# --- HTML 生成逻辑 ---
def generate_html_content(all_news_data_sorted):
    # (此函数内容与 diabetes_news_fetch_tabs_v1 中的 generate_html_content 完全相同)
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
        needs_translation = feed_info.get("needs_translation", False)
        raw_articles_from_feed = fetch_articles_from_rss(feed_info["url"], feed_info["source_override"])
        
        for article_data in raw_articles_from_feed:
            if article_data["url"] in globally_seen_urls: continue
            title_to_process = article_data["title"]
            snippet_to_process = article_data["snippet"]
            if needs_translation:
                print(f"    需要翻译来自 {feed_info['source_override']} 的文章: {title_to_process[:30]}...")
                title_to_process = translate_text(title_to_process)
                snippet_to_process = translate_text(snippet_to_process)
                time.sleep(0.5) 

            if is_within_last_month_rss(article_data["time_struct"], today):
                normalized_title = normalize_title(title_to_process)
                time_display_str = "未知时间"
                if article_data["time_struct"]:
                    try: time_display_str = time.strftime("%Y-%m-%d", article_data["time_struct"])
                    except: pass
                
                article_obj_for_storage = {
                    "title": title_to_process, "url": article_data["url"], "snippet": snippet_to_process, 
                    "source": article_data["source"], "time_display_str": time_display_str, 
                    "time_struct": article_data["time_struct"], 
                    "source_priority": current_priority,
                    "source_type": "authoritative_rss" # 添加来源类型标记
                }

                if normalized_title not in unique_articles_candidates or \
                   current_priority > unique_articles_candidates[normalized_title]["priority"]:
                    unique_articles_candidates[normalized_title] = {
                        "article_obj": article_obj_for_storage, "priority": current_priority,
                        "target_categories": feed_info["target_categories"], "url": article_data["url"]
                    }
                    # print(f"      候选(权威RSS): '{article_obj_for_storage['title'][:30]}...' (Prio: {current_priority})")
        time.sleep(1)

    # --- 步骤二：从爬虫源获取新闻 ---
    print("\n--- 正在从爬虫源获取新闻 ---")
    for scraper_info in SCRAPED_SOURCES_CONFIG:
        if scraper_info["fetch_function"] not in SCRAPER_FUNCTIONS_MAP:
            print(f"    错误: 未在 SCRAPER_FUNCTIONS_MAP 中找到爬虫函数 {scraper_info['fetch_function']}")
            continue
        fetch_function = SCRAPER_FUNCTIONS_MAP[scraper_info["fetch_function"]]
        current_priority = scraper_info.get("priority", 3)
        raw_articles_from_scraper = fetch_function()
        
        for article_data in raw_articles_from_scraper: 
            if article_data["url"] in globally_seen_urls: continue
            if is_within_last_month_rss(article_data["time_struct"], today):
                normalized_title = normalize_title(article_data["title"])
                time_display_str = "未知时间"
                if article_data.get("time_struct"):
                    try: time_display_str = time.strftime("%Y-%m-%d", article_data["time_struct"])
                    except: pass
                else: 
                    print(f"      注意: 文章 '{article_data['title'][:30]}...' 无有效发布日期，将显示'未知时间'")
                article_obj_for_storage = {
                    "title": article_data["title"], "url": article_data["url"], "snippet": article_data["snippet"], 
                    "source": scraper_info["source_override"], "time_display_str": time_display_str, 
                    "time_struct": article_data["time_struct"], 
                    "source_priority": current_priority,
                    "source_type": "scraper" # 添加来源类型标记
                }
                if normalized_title not in unique_articles_candidates or \
                   current_priority > unique_articles_candidates[normalized_title]["priority"]:
                    unique_articles_candidates[normalized_title] = {
                        "article_obj": article_obj_for_storage, "priority": current_priority,
                        "target_categories": scraper_info["target_categories"], "url": article_data["url"]
                    }
                    # print(f"      候选(爬虫): '{article_obj_for_storage['title'][:30]}...' (Prio: {current_priority})")
            # else:
            #      print(f"      跳过(爬虫，日期过滤): '{article_data['title'][:30]}...' (Time struct: {article_data.get('time_struct')})")
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
                    "title": article_data["title"], "url": article_data["url"], "snippet": article_data["snippet"], 
                    "source": article_data["source"], "time_display_str": time_display_str, 
                    "time_struct": article_data["time_struct"], 
                    "source_priority": GOOGLE_NEWS_PRIORITY,
                    "source_type": "google_news" # 添加来源类型标记
                }
                if normalized_title not in unique_articles_candidates or \
                   GOOGLE_NEWS_PRIORITY > unique_articles_candidates[normalized_title]["priority"]:
                    unique_articles_candidates[normalized_title] = {
                        "article_obj": article_obj_for_storage, "priority": GOOGLE_NEWS_PRIORITY,
                        "target_categories": [site_category_name], "url": article_data["url"]
                    }
                    # print(f"      候选(Google): '{article_obj_for_storage['title'][:30]}...' (Prio: {GOOGLE_NEWS_PRIORITY}) for category {site_category_name}")
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

    # --- 步骤五：对每个分类的文章按来源类型和日期排序并截取 ---
    print("\n--- 正在对各分类新闻进行排序和截取 ---")
    all_articles_by_site_category_final_sorted = {}
    for category_name, articles_list in all_articles_by_site_category_temp.items():
        # 复合排序：首先按来源类型排序 (权威RSS -> 爬虫 -> Google News)，然后按日期降序
        articles_list.sort(key=lambda x: (
            SOURCE_TYPE_ORDER.get(x.get("source_type", "unknown"), 99), # 按来源类型排序值升序
            -(time.mktime(x["time_struct"]) if x.get("time_struct") else -float('inf')) # 按日期时间戳降序
        ))
        all_articles_by_site_category_final_sorted[category_name] = articles_list[:MAX_ARTICLES_PER_CATEGORY]
        print(f"  分类 '{category_name}' 排序并截取后有 {len(all_articles_by_site_category_final_sorted[category_name])} 条新闻。")
        # 打印排序后的前几条来源类型，用于验证
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
