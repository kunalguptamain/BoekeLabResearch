from bs4 import BeautifulSoup
import re
import ast
import csv

factor_save_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/transcription_factor_analysis/factors_7.txt"
factor_csv_save_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/transcription_factor_analysis/factors_7.csv"

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

dicts = [{}, {}, {}, {}, {}]
total = [0, 0, 0, 0, 0]

classification_dict = {
    4: "ALL",
    3: "GENERATED",
    2: "STRONG",
    1: "WEAK",
    0: "NON",
}

for i, entry in enumerate(data):
    print(i)
    transcription_factor_response = ast.literal_eval(entry)
    html = transcription_factor_response["request"]
    classification = int(transcription_factor_response['classification'])
    sort_factors(html, classification, dicts)
    total[classification] += 1
    total[-1] += 1

print(total)
factors = list(dicts[-1].keys())
counts = list(dicts[-1].values())

sorted_indices = sorted(range(len(counts)), key=lambda i: counts[i])

sorted_factors = [factors[i] for i in sorted_indices]
sorted_counts = [counts[i] for i in sorted_indices]

with open(factor_csv_save_path, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow([''] + [f'{classification_dict[i]}' for i in range(len(dicts))])
    for key in sorted_factors:
        row = [key] + [d.get(key, 0) / total[i] for i, d in enumerate(dicts)]
        writer.writerow(row)
