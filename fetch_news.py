# fetch_news.py
import datetime
import json
import html # 用于转义HTML特殊字符
import os # 引入os模块用于读取环境变量 (为将来使用API Key做准备)

# --- 配置资讯分类 ---
CATEGORIES_CONFIG = {
    "最新研究": {"keywords": ["diabetes research", "糖尿病 研究"], "emoji": "🔬"},
    "治疗进展": {"keywords": ["diabetes treatment", "糖尿病 治疗"], "emoji": "💊"},
    "预防与生活方式": {"keywords": ["diabetes prevention lifestyle", "糖尿病 预防 生活方式"], "emoji": "🏃‍♀️"},
    "并发症管理": {"keywords": ["diabetes complications", "糖尿病 并发症"], "emoji": "🩺"},
    "饮食与营养": {"keywords": ["diabetes diet nutrition", "糖尿病 饮食 营养"], "emoji": "🥗"}
}

# --- 模拟新闻获取 ---
# 重要提示：在实际应用中，您需要替换下面的函数，
# 使用真实的新闻API（如 NewsAPI, Google News API 等，可能需要注册并获取API Key）
# 并将API Key存储在GitHub Secrets中 (例如 NEWS_API_KEY), 然后通过 os.getenv('NEWS_API_KEY') 读取。
def fetch_mock_news_for_category(category_name, keywords):
    """
    模拟获取新闻数据。实际应用中应替换为API调用。
    确保返回的文章包含 'title', 'url', 'snippet', 'source', 'time' 字段。
    'time' 字段应尽量模拟真实的时间字符串，如 "X 小时前", "X 天前", "YYYY-MM-DD"。
    """
    print(f"模拟获取 '{category_name}' 分类的新闻...")
    # 示例数据结构
    mock_articles = [
        {
            "title": f"{category_name}的示例新闻1：一项重大发现", # 示例标题
            "url": "javascript:void(0);", # 使用 javascript:void(0); 作为占位符链接
            "snippet": f"这是关于“{category_name}”的一条模拟新闻摘要。内容关注最新的动态和进展，帮助您了解行业前沿。", # 示例摘要
            "source": "模拟新闻来源A", # 示例来源
            "time": f"{datetime.datetime.now().day % 7 + 1} 天前" # 模拟一周内的时间
        },
        {
            "title": f"{category_name}的示例新闻2：专家深度解读", # 示例标题
            "url": "javascript:void(0);",
            "snippet": f"来自“{category_name}”领域的专家对近期热点进行了深入解读，提供了宝贵的见解和分析。", # 示例摘要
            "source": "模拟新闻来源B", # 示例来源
            "time": f"{datetime.datetime.now().hour % 24} 小时前" # 模拟小时前
        },
        {
            "title": f"{category_name}的示例新闻3：未来发展趋势探讨", # 示例标题
            "url": "javascript:void(0);",
            "snippet": f"探讨“{category_name}”方向的未来发展趋势，以及可能带来的影响和机遇。", # 示例摘要
            "source": "模拟新闻来源C", # 示例来源
            "time": "昨天" # 模拟昨天
        }
    ]
    # 为每个分类返回3-5条模拟新闻 (数量随机以模拟变化)
    return mock_articles[:(datetime.datetime.now().second % 3 + 2)] # 保证至少返回2条


# --- HTML 生成逻辑 ---
def generate_html_content(all_news_data):
    """根据新闻数据生成完整的HTML页面内容"""
    
    # 获取当前时间用于显示更新时间
    # 尝试从环境变量获取时区，如果未设置，则默认为UTC
    tz_info = None
    try:
        app_timezone = os.getenv('APP_TIMEZONE', 'UTC') # 从环境变量读取时区，默认为UTC
        if app_timezone != 'UTC': # 如果不是UTC，尝试创建时区对象
             # 注意：datetime.timezone 需要 Python 3.9+ 对于 ZoneInfo
             # 对于更广泛的兼容性，可以考虑使用 pytz 库，但这会增加依赖
             # 这里简化处理，仅支持 UTC 或依赖系统正确解析
            pass # 实际项目中可能需要 pytz 或 dateutil
    except Exception:
        print("Warning: Could not determine timezone from APP_TIMEZONE. Defaulting to system/UTC for timestamps.")

    current_time = datetime.datetime.now(tz_info)
    current_time_str = current_time.strftime('%Y年%m月%d日 %H:%M:%S')
    current_year = current_time.year

    # HTML头部和样式 (Tailwind CSS via CDN)
    html_output = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>糖尿病前沿资讯</title>
    <script src="https://cdn.tailwindcss.com"></script>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&display=swap" rel="stylesheet">
    <style>
        body {{
            font-family: 'Inter', sans-serif;
            background-color: #f0f4f8; /* 淡雅背景色 */
        }}
        .news-card {{
            transition: transform 0.3s ease, box-shadow 0.3s ease;
            display: flex;
            flex-direction: column;
            justify-content: space-between;
            height: 100%;
            background-color: #ffffff; /* 卡片背景白色 */
            border-radius: 0.5rem; /* md rounded corners */
        }}
        .news-card:hover {{
            transform: translateY(-5px);
            box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05);
        }}
        .category-section {{
            background-color: #ffffff; /* 分类区块背景白色 */
            border-radius: 0.75rem; /* xl rounded corners */
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06); /* subtle shadow */
            padding: 1.25rem; /* p-5 */
        }}
        @media (min-width: 768px) {{
            .category-section {{
                padding: 2rem; /* md:p-8 */
            }}
        }}
        .category-title-text {{
            border-bottom: 3px solid #3b82f6; /* 蓝色下划线 */
            padding-bottom: 0.75rem; /* pb-3 */
            margin-bottom: 1.75rem; /* mb-7 */
            color: #1e40af; /* 深蓝色标题 */
            font-size: 1.5rem; /* text-2xl */
        }}
        @media (min-width: 768px) {{
            .category-title-text {{
                font-size: 1.875rem; /* md:text-3xl */
            }}
        }}
        .news-item-link {{
            color: #2563eb; /* 链接颜色 */
            text-decoration: none;
            font-weight: 600; /* font-semibold */
        }}
        .news-item-link:hover {{
            text-decoration: underline;
            color: #1d4ed8; /* 悬停时更深的蓝色 */
        }}
        .source-tag, .time-tag {{
            padding: 0.25rem 0.75rem; /* px-3 py-1 */
            border-radius: 0.375rem; /* rounded-md */
            font-size: 0.75rem; /* text-xs */
            font-weight: 500; /* font-medium */
            line-height: 1.25;
        }}
        .source-tag {{ background-color: #e0e7ff; color: #3730a3; }}
        .time-tag {{ background-color: #d1fae5; color: #065f46; }}
        .card-content {{ flex-grow: 1; }}
        .card-footer {{ margin-top: auto; /* Pushes footer to bottom */ }}
        .header-main-title {{
            font-size: 1.875rem; /* text-3xl */
        }}
        @media (min-width: 640px) {{ /* sm breakpoint */
            .header-main-title {{ font-size: 2.25rem; }} /* sm:text-4xl */
        }}
        @media (min-width: 768px) {{ /* md breakpoint */
            .header-main-title {{ font-size: 3rem; }} /* md:text-5xl */
        }}
        .main-container {{ padding: 1rem; }} /* p-4 */
        @media (min-width: 768px) {{ /* md breakpoint */
            .main-container {{ padding: 2rem; }} /* md:p-8 */
        }}
    </style>
</head>
<body class="bg-gray-100 text-gray-800">
    <div class="container mx-auto main-container">
        <header class="text-center mb-10 md:mb-16">
            <h1 class="font-bold text-blue-700 header-main-title">糖尿病前沿资讯</h1>
            <p class="text-gray-600 mt-3 text-base md:text-lg">本周最新动态（自动更新于：{current_time_str}）</p>
            <p class="text-sm text-orange-600 mt-2 font-semibold">提示：当前显示为模拟数据。请在实际部署时替换为真实新闻源。</p>
        </header>
        <div id="news-content" class="space-y-12">
"""

    # 生成每个分类的新闻内容
    if not all_news_data or all(not articles for articles in all_news_data.values()):
        html_output += '<p class="text-center text-gray-500 text-xl py-10">抱歉，目前未能加载到相关的糖尿病资讯。</p>'
    else:
        for category, articles in all_news_data.items():
            category_emoji = CATEGORIES_CONFIG.get(category, {}).get("emoji", "")
            html_output += f"""
            <section class="category-section">
                <h2 class="font-semibold category-title-text">{category_emoji} {html.escape(category)}</h2>
            """
            if not articles:
                html_output += '<p class="text-gray-500">本周暂无该分类下的资讯。</p>'
            else:
                html_output += '<div class="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-6 md:gap-8 news-grid">'
                for article in articles[:6]: # 每个分类最多显示6条
                    title = html.escape(article.get('title', '无标题'))
                    url = html.escape(article.get('url', 'javascript:void(0);')) # Default to void if no URL
                    snippet_raw = article.get('snippet', '暂无摘要')
                    snippet = html.escape(snippet_raw[:120] + ('...' if len(snippet_raw) > 120 else ''))
                    source = html.escape(article.get('source', '未知来源'))
                    time_display = html.escape(article.get('time', '未知时间'))

                    html_output += f"""
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
                html_output += "</div>" # close grid
            html_output += "</section>"

    # HTML尾部
    # 在页脚添加指向您的 GitHub 仓库的链接 (请记得替换占位符)
    github_repo_url = "https://github.com/doudou-ux/diabetes-news" # 请替换这里
    html_output += f"""
        </div> </div> <footer class="text-center p-6 mt-12 text-gray-600 text-sm border-t border-gray-300">
        <p>&copy; {current_year} 糖尿病资讯聚合. <a href="{github_repo_url}" target="_blank" rel="noopener noreferrer" class="text-blue-600 hover:underline">项目源码</a></p>
        <p class="mt-1">本站内容仅供参考, 不构成医疗建议。</p>
    </footer>
</body>
</html>
"""
    return html_output

# --- 主执行逻辑 ---
if __name__ == "__main__":
    print("开始生成糖尿病资讯网页...")
    all_news_data_for_html = {}

    for category_name_zh, config in CATEGORIES_CONFIG.items():
        articles = fetch_mock_news_for_category(category_name_zh, config["keywords"])
        all_news_data_for_html[category_name_zh] = articles
    
    # 生成HTML内容
    final_html = generate_html_content(all_news_data_for_html)
    
    # 将HTML内容写入index.html文件
    # 在GitHub Actions中，这个文件会被提交到仓库的根目录
    output_filename = "index.html" 
    try:
        # 确保脚本写入到相对于脚本位置的 index.html，或者指定绝对路径
        # 对于 GitHub Actions，通常写入到工作区的根目录即可
        with open(output_filename, "w", encoding="utf-8") as f:
            f.write(final_html)
        print(f"成功生成网页：{output_filename}")
    except IOError as e:
        print(f"错误：无法写入文件 {output_filename}。错误信息: {e}")
    except Exception as e:
        print(f"生成过程中发生未知错误: {e}")


    print("资讯网页生成完毕。")
