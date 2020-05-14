#! /usr/local/bin/python3.6

import os
import sys
import subprocess
import gzip
import re
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import multiprocessing
import traceback

								#######################################
													# Main #
								#######################################

def main():
	list_path=recup_path()
	### creation of required directories ###
	result_path=sys.argv[2]
	if(result_path[-1]!="/"):
		result_path=result_path+"/"
	os.makedirs(result_path+"fichierpoubelle", exist_ok=True)
	os.makedirs(result_path+"archive", exist_ok=True)
	os.makedirs(result_path+"results", exist_ok=True)
	### multiprocessing of run_softwares ###
	multithreading(list_path)
	shutil.rmtree(result_path+"fichierpoubelle")
	### creation of a result file gathering strain result files ###
	toutenun(result_path+"results")
	#insertion db


									######################################
													# MacSyFinder #
									######################################

def macsyfinder_wrapper(pathrep, result_path) :
	### recovery of the file name from the path ###
	path_input=pathrep+pathrep.split("/")[-2]+"_cds_from_genomic.fna"
	faafile=result_path+"fichierpoubelle/"+pathrep.split("/")[-2]+"/prot_from_cds.faa"
	os.makedirs(result_path+"fichierpoubelle/"+pathrep.split("/")[-2], exist_ok=True)
	outdir=result_path+"archive/"+pathrep.split("/")[-2]+"/MSF"
	logfile = open(result_path+"archive/"+pathrep.split("/")[-2]+"/MSF.log.txt", 'w')
	if not os.path.exists(outdir+".tar.gz"):
		if not os.path.exists(outdir):
			### translation from CDS is needed to take in account pseudogenes and duplication ###
			translationCDS(path_input, faafile)
			### execution of MacSyFinder ###
			msfcmd=["macsyfinder", "--sequence-db", faafile, "--db-type", "ordered_replicon", "-d", "DataGetter/My_macsyfinder_profils/DEF/", "-p", "DataGetter/My_macsyfinder_profils/Profiles/", "-o", outdir, "all"]
			print(' '.join(msfcmd))
			subprocess.run(msfcmd, stderr=logfile, stdout=subprocess.PIPE)	
			shutil.move(faafile, outdir)
		### recovery of data ###
		recup_macsyfinder(path_input, outdir, result_path)
	if os.path.exists(result_path+"fichierpoubelle/"+pathrep.split("/")[-2]):
		shutil.rmtree(result_path+"fichierpoubelle/"+pathrep.split("/")[-2])
	logfile.close()

	

def translationCDS(path_input, faafile):
	CDScounter=1
	translation=list()
	dico_id=dict()
	accession_research=re.compile('[YNPW]{2}_[0-9]+\.[0-9]')
	with gzip.open(path_input+".gz", "rt") as handle:
		for record in SeqIO.parse(handle, 'fasta') :
			record.IUPAC='unambiguous_dna'
			accession_num=accession_research.findall(record.id)		
			translation.append(record.translate())
			if (accession_num!=[] and accession_num[0] in dico_id.keys()):
				translation[CDScounter-1].id=accession_num[0]+"-"+str(dico_id[accession_num[0]])
				dico_id[accession_num[0]]+=1
			elif accession_num!=[]:
				translation[CDScounter-1].id=accession_num[0]
				dico_id[accession_num[0]]=1
			else:
				translation[CDScounter-1].id="pseudogene"+str(CDScounter)
			CDScounter+=1
	SeqIO.write(translation, faafile, "fasta")



def recup_macsyfinder(path_input, outdir, result_path):
	if os.path.exists(outdir+"/hmmer_results") :
		shutil.rmtree(outdir+"/hmmer_results")
		os.remove(outdir+"/macsyfinder.out")
		os.remove(outdir+"/macsyfinder.log")
		os.remove(outdir+"/macsyfinder.tab")
		os.remove(outdir+"/macsyfinder.conf")
	reportfile=open(outdir+"/macsyfinder.report", 'r')
	dict_acc_data=dict()
	for l in reportfile:
		if '#' not in l:
			list_data=l.rstrip('\n').split("\t")
			if "pseudo" in list_data[0] :
				key="cds_"+list_data[2]
			else:
				key=list_data[0].split("-")[0]+"_"+list_data[2]
			dict_acc_data[key]=[list_data[7], list_data[14], list_data[15]]
	location_research=re.compile("\[location=.+\]")
	coordiate_research=re.compile("[0-9]+") 
	acc_research=re.compile("[A-Z]{2}_[A-Z0-9]+\.[0-9]")
	with gzip.open(path_input+".gz", "rt") as handle:
		for record in SeqIO.parse(handle, 'fasta') :
			for key, val in dict_acc_data.items():
				accnum_prot=re.compile(key)
				accnum_replicon=accnum_prot.findall(record.id)
				if accnum_replicon!=[]:
					localisation=location_research.findall(record.description)
					posdep=coordiate_research.findall(localisation[0])
					genomeid=acc_research.findall(record.id)
					val[1]=int(posdep[0])+int(val[1])
					val[2]=int(posdep[0])+int(val[2])
					val.append(genomeid[0])
	reportfile.close()
	dict_without_duplication=dict()
	for key, val in dict_acc_data.items() :
		key2=val[0]
		if key2 in dict_without_duplication:
			val2=dict_without_duplication[key2]
			val2[2]=val[2]
			val2[4].append(key)
		else:
			dict_without_duplication[key2]=val
			val2=dict_without_duplication[key2]
			val2.append(0)
			val2[4]=[key]
	tabfile=open(result_path+"results/"+outdir.split("/")[-2]+".txt", 'a')
	for k, v in dict_without_duplication.items():
		tabfile.write("MSF\t"+str(v[3])+"\t"+str(v[1])+"\t"+str(v[2])+"\t"+'_'.join(k.split(_)[1:-1])+"\t"+str(v[4])+"\n")		# other info = liste des id des proteines qui constitue le systeme macromoleculaire
	tabfile.close()
	if (os.path.exists(outdir) and os.path.exists(outdir+".tar.gz")==False):
		shutil.make_archive(result_path+"archive/"+outdir.split("/")[-2]+"/MSF", 'gztar', root_dir=outdir.replace("/MSF", ""), base_dir='MSF')
		shutil.rmtree(outdir)

									######################################
												       # Alien_Hunter #
									######################################

def alienhunter_wrapper(pathrep, result_path) :
	### recovery of the file name from the path ###
	path_fna=pathrep+pathrep.split("/")[-2]+"_genomic.fna" 
	outdir=result_path+"archive/"+pathrep.split("/")[-2]+"/AH/"
	### creation of required directories ###
	os.makedirs(outdir, exist_ok=True)
	### execution of Alien_Hunter ###
	logfile = open(outdir+"/AH.log.txt", 'w')
	with gzip.open(path_fna+".gz", "rt") as handle:
		for record in SeqIO.parse(handle, 'fasta'):
			if not os.path.exists(outdir+record.id+".tar.gz"):
				fileout=outdir+record.id+"/result"+record.id+".embl"
				if not os.path.exists(fileout):
					SeqIO.write(record, result_path+"fichierpoubelle/"+record.id+".fasta", "fasta")
					os.makedirs(outdir+record.id, exist_ok=True)
					ahcmd = ["alien_hunter", result_path+"fichierpoubelle/"+record.id+".fasta", fileout]
					print(' '.join(ahcmd))
					subprocess.run(ahcmd, stderr=logfile, stdout=subprocess.PIPE)
					os.remove(result_path+"fichierpoubelle/"+record.id+".fasta")
				### recovery data ###
				if os.path.exists(fileout):
					resultAH=open(fileout, 'r')
					location_research=re.compile("[0-9]+")
					score_research=re.compile("[0-9]+\.[0-9]+")
					dict_numacc_data=dict()
					countPGI=1
					for l in resultAH:
						if "misc_feature" in l:
							key=record.id+"_PGI"+"_"+str(countPGI)
							dict_numacc_data[key]=list()	#potential genetic island (PGI)
							val=dict_numacc_data[key]
							val.append(location_research.findall(l)[0])				#start position in the genome
							val.append(location_research.findall(l)[1])				#end position in the genome
						elif "score" in l:
							key=record.id+"_PGI"+"_"+str(countPGI)
							val=dict_numacc_data[key]
							val.append(score_research.findall(l)[0])					#score
							countPGI+=1
					resultAH.close()
					### fill final file for this strain ###
					tabfile=open(result_path+"results/"+pathrep.split("/")[-2]+".txt", 'a')
					for key, val in dict_numacc_data.items():
						tabfile.write("AH\t"+key.split("_")[0]+"_"+key.split("_")[1]+"\t"+str(val[0])+"\t"+str(val[1])+"\t"+key.split("_")[-2]+"\t"+str(val[2])+"\n")				#other info = score
					tabfile.close()
				else:
					logfile.write("Pas de resultat pour "+record.id+"\n")
				shutil.make_archive(outdir+record.id, 'gztar', root_dir=outdir, base_dir=record.id)
				shutil.rmtree(outdir+record.id)
	logfile.close()

								#######################################
												# Integron_Finder #
								#######################################

def integronfinder_wrapper(pathrep, result_path) :
	file_fna=pathrep+pathrep.split("/")[-2]+"_genomic.fna" 
	os.makedirs(result_path+"archive/"+pathrep.split("/")[-2]+"/IF", exist_ok=True)
	logfile = open(result_path+"archive/"+pathrep.split("/")[-2]+"/IF/IF.log.txt", 'w')
	with gzip.open(file_fna+".gz", "rt") as handle:
		for record in SeqIO.parse(handle, 'fasta') :
			idrecord=record.id
			if not os.path.exists(result_path+"archive/"+pathrep.split("/")[-2]+"/IF/"+idrecord+".tar.gz "):
				SeqIO.write(record, result_path+"fichierpoubelle/"+idrecord+".fasta", "fasta")
				if not os.path.exists(result_path+"archive/"+pathrep.split("/")[-2]+"/IF/Results_Integron_Finder_"+idrecord+"/"+idrecord+".integrons"):
					ifcmd=["integron_finder",result_path+"fichierpoubelle/"+idrecord+".fasta"]
					print(' '.join(ifcmd))
					subprocess.run(ifcmd, stderr=logfile, stdout=subprocess.PIPE)
					os.remove(result_path+"fichierpoubelle/"+idrecord+".fasta")
					shutil.move("Results_Integron_Finder_"+idrecord, result_path+"archive/"+pathrep.split("/")[-2]+"/IF/")
				recup_integron(idrecord, pathrep, result_path)
				shutil.make_archive(result_path+"archive/"+pathrep.split("/")[-2]+"/IF/"+idrecord, 'gztar', root_dir=result_path+"archive/"+pathrep.split("/")[-2]+"/IF/", base_dir="Results_Integron_Finder_"+idrecord)
				shutil.rmtree(result_path+"archive/"+pathrep.split("/")[-2]+"/IF/"+"Results_Integron_Finder_"+idrecord)
	logfile.close()

def recup_integron(id, pathrep, result_path):
	resultIF=open(result_path+"archive/"+pathrep.split("/")[-2]+"/IF/Results_Integron_Finder_"+id+"/"+id+".integrons")
	dico_id_data=dict()
	for line in resultIF:
		if "# No Integron found" not in line:
			if "ID_integron" not in line:
				l=line.rstrip('\n').split("\t")
				key=l[1]+"_"+l[0]
				if key in dico_id_data.keys():
					dico_id_data[key][1]=l[4]
				else:
					dico_id_data[key]=[l[3], l[4], l[10]]
	tabfile=open(result_path+"results/"+pathrep.split("/")[-2]+".txt", 'a')
	for key, val in dico_id_data.items():
		if (str(val[2])=="CALIN" or str(val[2])=="In0"):
			val[2]="Partial_integron"
		else:
			val[2]="Complete_integron"
		tabfile.write("IF\t"+key.split("_")[0]+"_"+key.split("_")[1]+"\t"+str(val[0])+"\t"+str(val[1])+"\t"+str(val[2])+"\n")		#colonne other info vide
	tabfile.close()


									######################################
													# Other functions #
									######################################

def recup_path():
	path=sys.argv[1]
	if(path[-1]!="/"):
		path=path+"/"
	list_path=list()
	for l in os.listdir(path):
		list_path.append(path+l)
	return list_path

def multithreading(list):
	with multiprocessing.Pool(multiprocessing.cpu_count()) as p :	
		it=p.imap(run_softwares, list)
		for i, truc in enumerate(it):
			print(i)

def run_softwares(path):
	result_path=sys.argv[2]
	if result_path[-1]!="/":
		result_path=result_path+"/"
	if path[-1]!="/":
		path=path+"/"
	os.makedirs(result_path+"archive/"+path.split("/")[-2], exist_ok=True)
	if not os.path.exists(result_path+"results/"+path.split("/")[-2]+".txt"):
		tab=open(result_path+"results/"+path.split("/")[-2]+".txt", 'w')
		tab.write("#method\tid\tposition_start\tposition_end\ttype\tother_info\n")
		tab.close()
	try:
		macsyfinder_wrapper(path, result_path)
		alienhunter_wrapper(path, result_path)
		integronfinder_wrapper(path, result_path)
	except Exception as  e:
		traceback.print_exc(file=sys.stdout)
		sys.stdout.flush()
		raise e

def toutenun_pour_R(cheminversouestlefichierresults):
	if(cheminversouestlefichierresults[-1]!="/"):
		cheminversouestlefichierresults=cheminversouestlefichierresults+"/"
	fichierfinal=open(cheminversouestlefichierresults+"all_results_pour_R.tab", 'w')
	fichierfinal.write("code_assemblie\tgenre\tespece\ttype_replicon\ttaille_replicon\tmethod\tid\tposition_start\tposition_end\ttype\n")				#on a pas besoin de la colonne other_info dans R
	dict_ca_ge=crea_dict_codeassembly_genreespece()
	dict_nar_tetr=crea_dict_numassreplicon_typeettaillereplicon()
	for namefile in os.listdir(cheminversouestlefichierresults):
		if "all_results.tab" not in namefile:
			f=open(cheminversouestlefichierresults+namefile, 'r')
			for l in f:
				if "#" not in l :
					codeassembly=namefile.split('_')[0]+'_'+namefile.split('_')[1]
					numaccrep=l.split('\t')[1]
					if numaccrep in dict_nar_tetr:
						genre=dict_ca_ge[codeassembly].split(' ')[0]
						espece=dict_ca_ge[codeassembly].split(' ')[1]
						type_replicon=dict_nar_tetr[numaccrep][0]
						taille_replicon=dict_nar_tetr[numaccrep][1]
						print(l.split('\t')[0])
						if "IF" not in l.split('\t')[0]:
							fichierfinal.write(codeassembly+"\t"+genre+"\t"+espece+"\t"+type_replicon+"\t"+taille_replicon+"\t"+'\t'.join(l.split('\t')[0:-1])+"\n")
						else:
							fichierfinal.write(codeassembly+"\t"+genre+"\t"+espece+"\t"+type_replicon+"\t"+taille_replicon+"\t"+l)	
			f.close()
	fichierfinal.close()

def toutenun(cheminversouestlefichierresults):
	if(cheminversouestlefichierresults[-1]!="/"):
		cheminversouestlefichierresults=cheminversouestlefichierresults+"/"
	fichierfinal=open(cheminversouestlefichierresults+"all_results.tab", 'w')
	fichierfinal.write("method\tid\tposition_start\tposition_end\ttype\tother_info\n")
	for namefile in os.listdir(cheminversouestlefichierresults):
		if "all_results.tab" not in namefile:
			f=open(cheminversouestlefichierresults+namefile, 'r')
			for l in f:
				if "#" not in l :
					fichierfinal.write(l)
			f.close()
	fichierfinal.close()

##### creation dico 			###
def crea_dict_numassreplicon_typeettaillereplicon():
	f=open("/dipro/clemplacement/PanteroDB_v0.3/03.database/genome_replicons.tab","r")
	dico_path=dict()
	for ligne in f :
		val=ligne.rstrip('\n').split("\t")
		dico_path[val[1]]=[val[3], val[4]]
	return (dico_path)

def crea_dict_codeassembly_genreespece():
	f=open("/dipro/clemplacement/PanteroDB_v0.3/03.database/genome_assemblies.tab","r")
	dico_path=dict()
	for ligne in f :
		val=ligne.split("\t")
		dico_path[val[0]]=val[3]
	return (dico_path)

		###

main()