import urllib3
import selenium
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.chrome.service import Service
import selenium.webdriver.common.driver_finder
import logging
import os

loggingHandlers = [] 
loggingHandlers.append(logging.StreamHandler())
formatter = logging.Formatter("%(asctime)s;%(levelname)s;%(message)s",
                              "%Y-%m-%d %H:%M:%S")
loggingHandlers[0].setFormatter(formatter)
logging.getLogger().setLevel(logging.DEBUG)

options = Options()
options.binary_location = "/Applications/Google\ Chrome.app/Contents/MacOS/Google\ Chrome"
#s = Service(os.path.expanduser("~/Downloads/chromedriver_mac_arm64/chromedriver"))
#driver = webdriver.Chrome(options=options, service=s)
#driver = webdriver.Chrome(options=options)
driver = webdriver.Chrome()

driver.get("https://naif.jpl.nasa.gov/pub/naif/pds/pds4/maven/maven_spice/spice_kernels/ck/")
driver.get("https://naif.jpl.nasa.gov/pub/naif/pds/pds4/maven/maven_spice/spice_kernels/ck/mvn_app_pred_141205_141209_v01.bc")

driver.get("https://www.baidu.com")

http = urllib3.PoolManager()
url_ = 'https://naif.jpl.nasa.gov/pub/naif/pds/pds4/maven/maven_spice/spice_kernels/ck/mvn_app_pred_141205_141209_v01.bc'
headers = {
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7",
        "Accept-Encoding": "gzip, deflate, br, zstd",
        "Accept-Language": "en,zh-CN;q=0.9,zh;q=0.8",
        "Connection": "keep-alive",
        "Host": "naif.jpl.nasa.gov",
        "Referer": "https://naif.jpl.nasa.gov/pub/naif/pds/pds4/maven/maven_spice/spice_kernels/ck/",
        "Sec-Fetch-Dest": "document",
        "Sec-Fetch-Mode": "navigate",
        "Sec-Fetch-Site": "same-origin",
        "Sec-Fetch-User": "?1",
        "Upgrade-Insecure-Requests": "1",
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/129.0.0.0 Safari/537.36",
        "sec-ch-ua": '''"Google Chrome";v="129", "Not=A?Brand";v="8", "Chromium";v="129"''',
        "sec-ch-ua-mobile": "?0",
        "sec-ch-ua-platform": '''"macOS"''',
        }
http = urllib3.PoolManager(headers=headers)
resp = http.request('GET', url_)
