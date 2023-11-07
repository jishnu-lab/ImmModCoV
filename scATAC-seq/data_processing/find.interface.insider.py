import pickle as pkl 
from collections import defaultdict 
import pandas as pd 
import sys 


# command line arguments 

hgi_infile = str(sys.argv[1]) # the file with the HGI mutations mapped to the respective AA 

interface_file = str(sys.argv[2]) # the pkl file with the dictionary containing the proteins mapped to the interface residues

interface_file_covid = str(sys.argv[3]) # the interface annotation file from the SARS-COV-2 3D network paper 

outfile = str(sys.argv[4]) # the filename where the human proteins of interest with SNPs in the Hu-Hu interface 

covid_pkl_file = str(sys.argv[5]) # the PPI 2 residue information dictionary pickled file

viral_out = str(sys.argv[6]) # the file where the human proteins of interest with SNPs in the Hu-viral interface 

rsid_outfile = str(sys.argv[7])

data_hgi = pd.read_csv(hgi_infile,header=None,sep='\t')

prot2res_hgi = list(zip(data_hgi[0],data_hgi[1],data_hgi[2]))

# without collapsing the protein isoforms -  it is not a linear mapping at an AA sequence level
new_prot2res_hgi = [(prot,(int(res),rsid)) for rsid,prot,res in prot2res_hgi] # collapsing the isoforms to the main protein

dict_prot2res_hgi = defaultdict(list)

prot2rsids = defaultdict(list)

for p,r in new_prot2res_hgi:
	dict_prot2res_hgi[p].append(r)
	

handle = open(interface_file,'rb')

dict_insider = pkl.load(handle)
handle.close()
check = []
check1 = []
#int_rsids = []
handle0 = open(rsid_outfile,'w')
#print(dict_prot2res_hgi)
#sys.exit()
for prot in dict_prot2res_hgi:
	int_rsids = []

	try:
		#common_res = list(set(dict_prot2res_hgi[prot]) & set(dict_insider[prot]))
		common_res = []
		for res,ids in set(dict_prot2res_hgi[prot]):
			if res in set(dict_insider[prot]):
				common_res.append(res)
				int_rsids.append(ids)
		if common_res:
			check.append((prot,common_res))
			check1.append(prot)
			for ids in int_rsids:
				handle0.write(str(prot)+'\t'+str(ids)+'\n')
				
	except KeyError:
		continue 

handle1  = open(outfile,'w')

for prot,res in check:
	handle1.write(str(prot)+'\t'+str(res)+'\n')

handle1.close()

data_covid = pd.read_csv(interface_file_covid,sep='\t',header=0)

uniprot_human = set(data_covid['UniProt Human'])

test = set(check1) & set(uniprot_human)

print("First degree interactors: "+ str(test))

handle2 = open(covid_pkl_file,'rb')

prot2viral_interface = defaultdict(list) 

ppi2res = pkl.load(handle2)

#print(ppi2res[:10])

for ppi,resid in ppi2res.items():
	_,human_prot = ppi 
	_,human_resid = resid 
	prot2viral_interface[human_prot].extend(human_resid)
	
print (list(dict_prot2res_hgi.items())[0])
print (list(prot2viral_interface.items())[0])
print("check above..")
handle3 = open(viral_out,'w')
outlist = []
test = []
#print(prot2viral_interface)
#sys.exit()
for prot in prot2viral_interface:
		
	try:
		hgi_res = []
		print(dict_prot2res_hgi[prot])
		for res, ids in dict_prot2res_hgi[prot]:
			hgi_res.append(res)
		print(hgi_res)
		#print(prot2viral_interface[prot])
		common_res = list(set(hgi_res) & set(prot2viral_interface[prot]))
		all_res = dict_prot2res_hgi[prot] + prot2viral_interface[prot]
		if common_res:
			outlist.append((prot,common_res))	
			
		if all_res:
			test.append((prot,common_res)) 		
	except KeyError:
		continue 


print(len(outlist))	

print("test ...")	
#print(test)





