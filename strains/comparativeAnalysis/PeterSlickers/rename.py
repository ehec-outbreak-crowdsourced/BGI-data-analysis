import fasta


def main():
	for feat in ("B7MQ09","wrbA",):
		infilename = "/home/piet/Ecol/Features/%s/Alignment-featDB.fna" % feat
		outfilename = "%s.fas" % feat
		fastatup_list = fasta.read_fasta_file_all(infilename)
		for i,fastatup in enumerate(fastatup_list):
			desc,seq = fastatup
			desc = desc.replace("LCL_10024.1","scaffold00001")
			desc = desc.replace("LCL_10024","scaffold00001")
			fastatup_list[i] = (desc,seq,)
		fasta.write_fasta_file_many(fastatup_list,outfilename)

if __name__ == "__main__":
	main()
