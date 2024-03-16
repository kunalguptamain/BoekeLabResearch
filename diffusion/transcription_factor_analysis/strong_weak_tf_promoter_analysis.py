from bs4 import BeautifulSoup
import re
import ast

factor_save_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/transcription_factor_analysis/factors_total_12.txt"

with open(factor_save_path, 'r') as f:
    loaded_array = f.readlines()
data = [item.strip() for item in loaded_array]

def extract_factors(data):
    soup = BeautifulSoup(data, 'html.parser')
    table = soup.find('table', {'class': 'n'})
    
    if table is not None:
        td_elements = table.find_all('td')
        factors = {}
        for i in range(0, len(td_elements), 2):
            number = td_elements[i].text
            factor = td_elements[i+1].find('font').text
            factors[(int(number))] = factor
    
    td_elements = soup.find_all('td', bgcolor=True)
    counts = [int(td.text.replace("\xa0", "")) for td in td_elements \
              if re.fullmatch(r'\d+', td.text.strip())]
    
    return factors, counts

def sort_factors(html, classification, net_dicts):
    factors, counts = extract_factors(html)
    for factor in factors.values(): 
        if factor not in net_dicts[classification]: net_dicts[classification][factor] = 0
        if factor not in net_dicts[-1]: net_dicts[-1][factor] = 0
    for count in counts:
        net_dicts[classification][factors[count]] += 1
        net_dicts[-1][factors[count]] += 1

dicts = [{}, {}, {}, {}]

for i, entry in enumerate(data):
    transcription_factor_response = ast.literal_eval(entry)
    html = transcription_factor_response["request"]
    classification = int(transcription_factor_response['classification'])
    sort_factors(html, classification, dicts)
    print(i)



