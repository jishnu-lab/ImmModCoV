import sys 
import pandas as pd 
import argparse 
from urllib.request import Request,urlopen
from urllib.parse import urlencode
import json 

# command line arguments for the input file 

infile = str(sys.argv[1]) # the file with the snps from the hgi consortium 
outfile = str(sys.argv[2]) # the output file where all the information will be stored

data = pd.read_csv(infile,sep='\t',header=None)

info_ids = data[1] # the ids in the format chr:position:ref:allele

phenotypes = data[2] # the different phenotypes based on the HGI studies 

ids2phenotypes = list(zip(info_ids,phenotypes))



url = 'http://bisque.yulab.org/cgi-bin/run.cgi'
handle = open(outfile,'w')

for ids, phen in ids2phenotypes:
	print(ids) 
	names = ids.split(":")
	chromo,pos,mut = names[0],names[1],names[2]+names[3]
							
	params = {
	'id': chromo,
	'output': 'uniprot',
	'position': pos,
	'mutation': mut,
	'type': 'hg38' # the build type for the input id 
	}				


                            
	data = urlencode(params).encode("utf-8")
	request = Request(url, data)

	

	response = urlopen(request)
	page = response.read(2000000)


	try:
		output_data = json.loads(page)
	
	except:
		continue 
	if len(output_data)!=0:
		for item in output_data:
			out_pos = item['output_position']	
			prot = item['output_identifier']
			out_mut = item['output_mutation']
			try:	
				handle.write(ids+'\t'+prot+'\t'+str(out_pos)+'\t'+out_mut+'\t'+phen+'\n')
			except TypeError:
				continue 


	


handle.close()




