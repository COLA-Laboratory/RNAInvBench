import requests
sequence_to_fold = "AUCUGUACUAGUUAGCUAACUAGAUCUGUAUCUGGCGGUUCCGUGGAAGAACUGACGUGUUCAUAUUCCCGACCGCAGCCCUGGGAGACGUCUCAGAGGC"
url = 'http://10.207.119.13:5000/ss_predict'
seq = sequence_to_fold
data = {'seq': [seq]}
response = requests.post(url, json=data)
print("response: ", response.json())
if response.status_code == 200:
    print(''.join(response.json()["predicted_sequence"]))