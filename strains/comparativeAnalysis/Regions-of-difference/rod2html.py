#!/usr/bin/env python

import sys


class Feature(object):
	pubmedbase = "http://www.ncbi.nlm.nih.gov/pubmed/"

	def __init__(self,sequence_id,start,end,direction,unigene,description):
		self.sequence_id = sequence_id
		self.start = start
		self.end = end
		self.direction = direction
		self.unigene = unigene
		self.description = description
		self.pmids = []

	def __str__(self):
		return "\t".join([str(i) for i in (self.sequence_id,self.start,self.end,self.direction,self.unigene,self.description,)])

	def pubmed_anchor(self):
		return ";".join([('<a href="%s%i">%i</a>' % (self.pubmedbase,pmid,pmid,)) for pmid in self.pmids])


class ROD(object):
	def __init__(self):
		self.features = []



def read_annotation_db(filename=None):
	if filename is None:
		filename = "../../../strains/TY2482/seqProject/BGI/annotations/era7bioinformatics/era7_BGI_V4_annotation/Era7_EHEC_BGI_V4_Annotation.txt"

	infile = file(filename)
	headers = infile.readline().split("\t")
	headers[-1] = headers[-1].rstrip()
	
	feature_dict = {}
	for line in infile:
		fields = line.split("\t")
		fields[-1] = fields[-1].rstrip()
		
		## sequence name is column A --> python index 0
		## start_pos is column E --> python index 4
		## end_pos is column G --> python index 6
		## unigene accession column I --> python index 8
		## "Protein names" column J --> python index 9
		## "PubMed ID" column U --> python index 19
		sequence_id = fields[0]
		start = int(fields[4])
		end = int(fields[6])
		direction = 1
		if end < start:
			start,end = end,start
			direction = -1
		unigene = fields[8]
		description = fields[9]
		feature = Feature(sequence_id,start,end,direction,unigene,description)
		feature.pmids = [int(i) for i in fields[20].split(";") if i.strip()]
		feature_dict[(sequence_id,start,end,)] = feature
	return feature_dict		


def get_features_in_rod(rod,feature_dict):
	## find features in ROD
	feature_list = []
	for key, feature in feature_dict.iteritems():
		if rod.sequence_id == feature.sequence_id:
			center = (feature.start + feature.end) / 2.0
			if (center > rod.start) and (center < rod.end):
				feature_list.append( (center,feature,) )

	feature_list.sort( reverse=(rod.direction < 0) )
	return feature_list


def get_rod_dict(feature_dict,filename=None):
	if filename is None:
		filename = 	"rod_list.csv"
	infile = file(filename)
	headers = infile.readline().split("\t")
	headers[-1] = headers[-1].rstrip()
	
	rod_dict = {}
	for line in infile:
		fields = line.split("\t")
		fields[-1] = fields[-1].rstrip()
		rod = rod_dict.setdefault(fields[0],ROD())
		rod.id = fields[0]
		rod.sequence_id = fields[1]
		rod.start =  int(fields[2])
		rod.end = int(fields[3])
		rod.direction = int(fields[4])
		rod.features.extend(get_features_in_rod(rod,feature_dict))
	return rod_dict


def dump_html(feature_list,outfile=sys.stdout):
	if not hasattr(outfile,"write"):
		outfile = file(outfile,"wb")
	outfile.write("<html><body><table>\n")	
	for center,feature in feature_list:
		if len(feature.unigene) == 6:
			url = "http://www.uniprot.org/uniref/?query=member:%s+identity:0.9&lucky=yes" % feature.unigene
			unigene_anchor = '<a href="%s">%s</a>' % (url,feature.unigene,)
		else:
			unigene_anchor = ''
		outfile.write("<tr><td>%s</td><td>%i</td><td>%i</td><td>%i</td><td>%s</td><td>%s</td><td>%s</td></tr>\n" %
		  (feature.sequence_id,feature.start,feature.end,feature.direction,unigene_anchor,
		   feature.description or feature.unigene,feature.pubmed_anchor(),) )	
	outfile.write("</table></body></html>\n")	


def main():

	feature_dict = read_annotation_db()
	rod_dict = get_rod_dict(feature_dict)

	for rod_id,rod in rod_dict.iteritems():
		print rod_id,len(rod.features)
		dump_html(rod.features,rod.id + ".html")


if __name__ == "__main__":
	main()
