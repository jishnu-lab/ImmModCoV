import pandas as pd 
import sys 
import argparse 
import os 
import multiprocessing 
import numpy as np
from functools import partial 

def snp_finder(path,alpha):
	name = path.split('/')[-1]
	name = name.split('_')[2]
	data = pd.read_csv(path,header=0,sep='\t',compression='infer',usecols=['rsid','all_inv_var_meta_p','SNP'],dtype={'rsid':str,'all_inv_var_meta_p':np.float64,'SNP':str})
	rsids = list(data['rsid'])
	p_vals = list(data['all_inv_var_meta_p'])
	snps = list(data['SNP'])
	names = [name for i in range(len(rsids))]
	ids2vals = list(zip(snps,rsids,p_vals,names))
	ids2vals1 = [(s,ids,vals,n) for s,ids,vals,n  in ids2vals if vals<float(alpha)]

	
	return ids2vals1

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='get the snps which pass the threshold')
	parser.add_argument('-d','--dir_name',help='The directory with  the SNP file(s)')
	parser.add_argument('-o','--out',help='The output filename')
	#args = parser.parse_args()
	parser.add_argument('-a','--alpha',help = 'The alpha threshold for the SNPs')
	args = parser.parse_args()
	path = args.dir_name
	filenames = []

	for root, directories, files in os.walk(path,topdown=False):
		for name in files:
			filenames.append(os.path.join(root,name))

	print('start multiprocessing..')	
	pool = multiprocessing.Pool()
	snps_info = pool.map(partial(snp_finder,alpha = args.alpha),filenames)
	
	final_snp_infos = [item for sublist in snps_info for item in sublist]  

	print('write to the file...')
	with open(args.out,'w') as handle:
		for s,ids,p,n  in final_snp_infos:
			
			 
			handle.write(str(s)+'\t'+str(ids)+'\t'+str(p)+'\t'+str(n)+'\n')


