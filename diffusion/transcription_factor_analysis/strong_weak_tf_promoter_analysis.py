from bs4 import BeautifulSoup
import re

factor_save_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/transcription_factor_analysis/factors.txt"

with open(factor_save_path, 'r') as f:
    loaded_array = f.readlines()
data = [item.strip() for item in loaded_array]

def extract_factors(data):
    soup = BeautifulSoup(data, 'html.parser')
    links = soup.find_all('a')
    factors = []
    for link in links:
        text = link.text.strip()
        if re.findall(r'\[.*?\]', text):
            factor = text.split()[0]
            # Find the parent 'td' element
            parent_td = link.find_parent('td')
            # Find all 'td' elements
            all_tds = soup.find_all('td')
            # Get the index of the parent 'td' element
            index = all_tds.index(parent_td)
            factors[index] = factor
    return factors

for dit in data:
    print(extract_factors(dit.response))
