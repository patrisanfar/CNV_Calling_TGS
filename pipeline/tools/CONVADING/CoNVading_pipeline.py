from subprocess import call 
import sys
import os
import re   
import glob  
import numpy as np   
import pandas as pd   

#Arguments
bam_files_folder = sys.argv[1]
bedfile = sys.argv[2]
output_folder = sys.argv[3]
project_name = sys.argv[4]    
fai = sys.argv[5]   

#Creating CoNVaDING output folder if doesn't exist    # NEW #
if not os.path.exists(output_folder):                
	os.makedirs(output_folder)                       

#Creating StartWithMatchScore folder if doesn't exist,
StartWithMatchScore = output_folder + 'StartWithMatchScore'
if not os.path.exists(StartWithMatchScore):
    os.makedirs(StartWithMatchScore)

#Creating Controls folder if doesn't exist
Controls = output_folder + 'Controls'
if not os.path.exists(Controls):
	os.makedirs(Controls)

#Creating StartWithBestScore folder if doesn't exist
StartWithBestScore = output_folder + 'StartWithBestScore'
if not os.path.exists(StartWithBestScore):
	os.makedirs(StartWithBestScore)

#Creating CNVs_calling_results folder if doesn't exist
CNVs_calling_results = output_folder + 'CNVs_calling_results'
if not os.path.exists(CNVs_calling_results):
	os.makedirs(CNVs_calling_results)

#Creating TargetQcList folder if doesn't exist
TargetQcList = output_folder + 'TargetQcList'
if not os.path.exists(TargetQcList):
	os.makedirs(TargetQcList)

#Creating FinalList folder if doesn't exist
FinalList = output_folder + 'FinalList'
if not os.path.exists(FinalList):
	os.makedirs(FinalList)

#Calculating coverages for all samples
call('perl CoNVaDING.pl \
	-mode StartWithBam \
	-inputDir ' + bam_files_folder + ' \
	-outputDir ' + StartWithMatchScore + ' \
	-bed ' + bedfile + ' \
	-controlsDir ' + Controls + ' \
	-useSampleAsControl',shell = True)


#Remove 'chr' prefix from the .normalized.coverage.txt files before StartWithMatchScore.,
normalized_files = glob.glob(os.path.join(StartWithMatchScore, '*.normalized.coverage.txt'))

for filepath in normalized_files:
    with open(filepath, 'r') as file:
        lines = file.readlines()
    with open(filepath, 'w') as file:
        for line in lines:
            if line.startswith('chr'):
                line = line.replace('chr', '', 1)
            file.write(line)

#Copying coverages folder to Controls folder
call('cp ' + StartWithMatchScore+ '/*.normalized.coverage.txt' + Controls, shell = True)


#Calculating N Best Controls from StartWithMatchScore folder

tifCounter = len(glob.glob1(Controls,"*.normalized.coverage.txt"))

if(tifCounter < 30):
	SampleN=str(tifCounter-1)
	call('perl CoNVaDING.pl -mode StartWithMatchScore \
 		-inputDir ' + output_folder + 'StartWithMatchScore \
 		-controlsDir ' + Controls + ' \
 		-outputDir ' + StartWithBestScore + ' \
 		-controlSamples '+ SampleN, shell = True)

else:
	call('perl CoNVaDING.pl \
		-mode StartWithMatchScore \
 		-inputDir ' + output_folder + 'StartWithMatchScore \
 		-controlsDir ' + Controls + ' \
 		-outputDir ' + StartWithBestScore + ' \
 		-controlSamples 30', shell = True)


if len(os.listdir(StartWithBestScore)) == 0:
	sampleCounter = str(len(glob.glob(output_folder + 'StartWithMatchScore/*normalized.coverage.txt')) - 1)
	print("Available number of samples:")
	print(sampleCounter+"\n")
	print("WARNING: The number of Control Samples will be set to %s" %(sampleCounter))
	call('perl CoNVaDING.pl \
 		-mode StartWithMatchScore \
 		-inputDir ' + output_folder + 'StartWithMatchScore \
 		-controlsDir ' + Controls + ' \
 		-outputDir ' + StartWithBestScore + ' \
 		-controlSamples '+ sampleCounter, shell = True)


# Making the CNVs calling
call('perl CoNVaDING.pl \
	-mode StartWithBestScore \
 	-inputDir ' + StartWithBestScore + ' \
 	-controlsDir ' + Controls + ' \
 	-ratioCutOffHigh=1.3 \
 	-outputDir ' + CNVs_calling_results,shell = True)

#Combining all shortlist and all totallist


def fileCombiner(filesdir,filespattern):
		
	print(filesdir)
	cwd = os.getcwd() # Get current working directory
	os.chdir(filesdir) # Set working directory where the files are
	filenames = glob.glob(filespattern) # Get all shortlists file names
	if len(filenames)!=0:
		cont = 0
		for i in filenames:
			sampleinfo = pd.read_csv(i, sep = '\t') # Read files
			if len(sampleinfo) != 0: # Execute only if CNVs were detected
				samplename = re.sub(r"\_.*$|\..*$", "", i) # Extract the name of the sample from the name of the file
				sampleinfo['SAMPLE_NAME'] = samplename
				cont = cont + 1
				if cont == 1: 
					dfout = sampleinfo
				else:
					dfout = dfout.append(sampleinfo)
		
		cols = dfout.columns.tolist() # Change SAMPLE_NAME column the first one
		cols = cols[-1:] + cols[:-1]
		dfout = dfout[cols]
		
		dfout.index = range(len(dfout)) # Reset row index
		os.chdir(cwd) # Set back to the previous working directory
		return(dfout)
	else:
		#print("\nERROR: PROBLEMS WITH CoNVaDING CNV CALLING")
		sys.exit(1)



shortout = fileCombiner(CNVs_calling_results, '*.shortlist.txt')
shortout.to_csv(output_folder + 'CoNVaDING' + ".shortlist.txt", sep='\t', index = False)

totalout = fileCombiner(CNVs_calling_results, '*.totallist.txt')
totalout.to_csv(output_folder + 'CoNVaDING' + ".totallist.txt", sep='\t', index = False)


#Combining shortlist with scores of totallist
shortout["SAMPLE_NAME"]=shortout["SAMPLE_NAME"].astype(str)
shortout["CHR"]=shortout["CHR"].astype(str)
shortout["START"]=shortout["START"].astype(int)
shortout["STOP"]=shortout["STOP"].astype(int)

totalout["SAMPLE_NAME"]=totalout["SAMPLE_NAME"].astype(str)
totalout["CHR"]=totalout["CHR"].astype(str)
totalout["START"]=totalout["START"].astype(int)
totalout["STOP"]=totalout["STOP"].astype(int)

# shortout=shortout.astype({"SAMPLE_NAME":str,
#                           "CHR":str,
#                           "START":int,
#                           "STOP":int})
# totalout=totalout.astype({"SAMPLE_NAME":str,
#                           "CHR":str,
#                           "START":int,
#                           "STOP":int})

combiout = shortout.copy()
for i in range(0,len(combiout)):

	sampleNameRef = combiout.loc[i, "SAMPLE_NAME"]
	chrRef = str(combiout.loc[i, "CHR"])
	startRef = combiout.loc[i, "START"]
	endRef = combiout.loc[i, "STOP"]
	
	rows = totalout[(totalout["SAMPLE_NAME"] == sampleNameRef) & (totalout["CHR"] == chrRef) & (totalout["START"] >= startRef) & (totalout["STOP"] <= endRef)]

	combiout.loc[i, "AUTO_RATIO"] = np.mean(rows["AUTO_RATIO"])
	combiout.loc[i, "AUTO_ZSCORE"] = np.mean(rows["AUTO_ZSCORE"])
	combiout.loc[i, "AUTO_VC"] = np.mean(rows["AUTO_VC"])
	combiout.loc[i, "GENE_RATIO"] = np.mean(rows["GENE_RATIO"])
	combiout.loc[i, "GENE_ZSCORE"] = np.mean(rows["GENE_ZSCORE"])
	combiout.loc[i, "GENE_VC"] = np.mean(rows["GENE_VC"])
	combiout.loc[i, "SHAPIRO-WILK"] = np.mean(rows["SHAPIRO-WILK"])

combiout.to_csv(output_folder + 'CoNVaDING' + ".results.txt", sep='\t', index = False)