import requests

classification_dict = {
    "STRONG": 2,
    "WEAK": 1,
    "NON": 0,
}

promoter_abs_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/strong_weak_promoters/promoter.txt"
non_promoter_abs_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/strong_weak_promoters/non_promoter.txt"
factor_save_path = "C:/Users/kunal/Documents/BoekeLabResearch/diffusion/transcription_factor_analysis/factors.txt"
cycle = 0
save_threshold = 100

def promo_request(sequence, similarity = 8, id_con = 171055014000):
    url = 'https://alggen.lsi.upc.es/cgi-bin/promo_v3/promo/promo.cgi'
    data = {
        'Dissim': f'{similarity}',
        'MultiSearch': '0',
        'dirDB': 'TF_8.3',
        'idCon': f'{id_con}',
        'Option': 3,
        'String': sequence,
        'File': '(binary)',
        'B1': 'Submit'
    }
    response = requests.post(url, files=data)
    return response.text

def save(factor_save_path):
    with open(factor_save_path, 'w') as f:
        for item in responses:
            f.write("%s\n" % item)

responses = []

with open(promoter_abs_path, 'r') as file:
    lines = file.readlines()
    for i, line in enumerate(lines):
        if(i % 2 == 0):
            placed = True
            tokens = line.strip().split()
            forward = tokens[2] == "FORWARD"
            if(tokens[-1] not in classification_dict):
                placed = False
                continue
            classification = classification_dict[tokens[-1]]
        else:
            cycle += 1
            if(not placed): continue
            sequence = line.strip()
            response = promo_request(sequence)
            responses.append({
                "request": response,
                "classification": classification,
                "index": i
            })
            if(cycle % save_threshold == 0): save(factor_save_path)


with open(non_promoter_abs_path, 'r') as file:
    lines = file.readlines()
    for i, line in enumerate(lines):
        cycle += 1
        sequence = line.strip()
        response = promo_request(sequence)
        responses.append({
            "request": response,
            "classification": classification_dict["NON"],
            "index": -1
        })
        if(cycle % save_threshold == 0): save(factor_save_path)