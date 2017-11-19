import csv
import json
 
def read_json(filename):
	return json.loads(open(filename).read())
 
def write_csv(data, filename):
	with open(filename, 'w') as outf:
		dw = csv.DictWriter(outf, data[0].keys())
		dw.writeheader()
		for row in data:
			dw.writerow(row)
 
write_csv(read_json('clinical.json'), 'clinical.csv')