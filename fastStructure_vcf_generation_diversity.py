#Input parameters: length (bp between sites), vcf file, fai file
import argparse, statistics

parser = argparse.ArgumentParser(description="Generates one .vcf from another with sites no closer than a set value, to use as 'unlinked' data for fastStructure.")
parser.add_argument('-i', type=str,  metavar='input_file', required=True, help='Path to the input .vcf file')
parser.add_argument('-o', type=str, metavar='output_file', required=True, help='Path to the output .vcf file')
parser.add_argument('-r', type=str, metavar='reference_fai', required=True, help='Path to the index .fai for the reference.')
parser.add_argument('-g', type=int, metavar='minimum_gap', required=True, help='Minimum gap size between sites (and reference contig/chromosome ends)')
parser.add_argument('-f', type=float, metavar='allele_frequency', required=True, help='Lowest acceptable alternate allele frequency per site. Use 0 to accept all sites. Use 0.2 to accept only sites with 20 percent or more alternate alleles')
parser.add_argument('-m', type=float, metavar='missing_data', required=True, help='Highest acceptable ammount of missing data per site. Use 1 to accept all sites. Use 0.2 to exclude sites with 20 percent or more missing data.')
args = parser.parse_args()

#Extract useable data from the .fai file (lengths and chromosome/contig names)
contig_names=[]
contig_lengths=[]
contig_lengths_int=[]

reference=open(args.r, 'r').readlines()
for count1, line1 in enumerate(reference):
	line_list=line1.split("\t")
	contig_names.append(line_list[0])
	contig_lengths.append(line_list[1])
	contig_lengths_int.append(int(line_list[1]))

allele_frequency=[]
missing_data=[]
ploidy=[]
first_line=True
vcf=open(args.i, 'r')	
for count3, line3 in enumerate(vcf):
	vcf_line=line3.split("\t")
	if len(vcf_line) > 2 and vcf_line[0] != "#CHROM":
		if first_line == True:
			for count4, line4 in enumerate(vcf_line[9:]):
				ploidy.append(len(line4.split("/")))
			first_line=False
		info=vcf_line[7].replace("AF=", "").split(";")
		allele_frequency.append(float(info[1]))
		info2=vcf_line[7].replace("AN=", "").split(";")
		missing_data.append(1-(float(info2[2])/sum(ploidy)))
vcf.close()
acceptable= [i for i in allele_frequency if i >= args.f]
unacceptable= [i for i in allele_frequency if i < args.f]

#Give me some sumary statistics about the reference fragments
print("\nMean contig lenght = "+str(statistics.mean(contig_lengths_int))+" bp")
print("Median contig lenght = "+str(statistics.median(contig_lengths_int))+" bp")
print("Mode of contig lenghts = "+str(statistics.mode(contig_lengths_int))+" bp")
print("Longest contig lenght = "+str(max(contig_lengths_int))+" bp")
print("Shortest contig lenght = "+str(min(contig_lengths_int))+" bp")
print("\nMean allele frequency = "+str(statistics.mean(allele_frequency))+" bp")
print("Median allele frequency = "+str(statistics.median(allele_frequency))+" bp")
print("Mode of allele frequency = "+str(statistics.mode(allele_frequency))+" bp")
print("Highest allele frequency = "+str(max(allele_frequency))+" bp")
print("Lowest allele frequency = "+str(min(allele_frequency))+" bp")
print("\nTotal sites = "+str(len(allele_frequency)))
print("Total sites with acceptable allele frequency = "+str(len(acceptable)))
print("Total sites excluded because of low allele frequency = "+str(len(unacceptable))+"\n")
print("Mean missing data = "+str(statistics.mean(missing_data))+" ("+str(statistics.mean(missing_data)*100)+")%")
print("Median missing data = "+str(statistics.median(missing_data))+" ("+str(statistics.median(missing_data)*100)+")%")
print("Mode missing data = "+str(statistics.mode(missing_data))+" ("+str(statistics.mode(missing_data)*100)+")%")
print("Most missing data at a single site = "+str(max(missing_data))+" ("+str(max(missing_data)*100)+")%")
print("Least missing data at a single site = "+str(min(missing_data))+" ("+str(min(missing_data)*100)+")%\n")

#prepare the output file
output_file=open(args.o, "w+")

count=0
vcf=open(args.i, 'r')
#pick a site in the vcf, first site = NameError exception
old_contig="random_shite!"
for count2, line2 in enumerate(vcf):
	current_line=line2.split('\t')
	if len(current_line) <2:
		output_file.write(line2)
	elif current_line[0]=="#CHROM":
		output_file.write(line2)
	else: #for all lines with actual data
		frequency_info=current_line[7].replace("AF=", "").split(";")
		missing_info=current_line[7].replace("AN=", "").split(";")
		current_contig=current_line[0]
		current_contig_length=contig_lengths[contig_names.index(current_contig)]
		if current_contig==old_contig: #if we already have a 1st position check that it is more than the minimum gap from the last site and the contig end
			if int(current_line[1]) > int(args.g)+int(old_position) and int(current_contig_length)-int(current_line[1]) > int(args.g) and float(frequency_info[1]) >= args.f and 1-(float(missing_info[2])/sum(ploidy)) <= args.m:
				output_file.write(line2)
				old_contig=current_line[0]
				old_position=current_line[1]
				count+=1
		elif current_contig!=old_contig: #start of a new contig, check it is more than the minimum gap from the contig start
			if int(current_contig_length)-int(current_line[1]) < int(current_contig_length)-int(args.g) and float(frequency_info[1]) >= args.f and 1-(float(missing_info[2])/sum(ploidy)) <= args.m:
				output_file.write(line2)
				old_contig=current_line[0]
				old_position=current_line[1]
				count+=1

print("Total sites in new vcf = "+str(count))
output_file.close()