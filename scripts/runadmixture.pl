#!/usr/bin/perl

#timestamp to create random number
my $timestamp =  1000000 + int(rand(9999 - 1000));

# first input is the 'test individuals' file name the second is reference populations
# the test individuals are merged with plink
# the reference pops are prepended with name AA_Ref_
# $cores is number of cores to run
my ($fileName, $Refpop, $cores) = @ARGV;
# the test file is in the same folder as script. the reference file in REFERENCEFILE folder
#number of individuals in the reference population file
$fileFam=$fileName.".fam";
$count = `wc -l < $fileFam`;
my ($first_num) = $count =~ /(\d+)/;

# this writes out the markers in the reference panel
`./plink --bfile REFERENCEFILES/$Refpop  --write-snplist`;

# this opens and loads the test file individuals
open($fileH,$fileFam);
@listOfIndivs;
while($line=<$fileH>){
	print $line."\n";
	push(@listOfIndivs,$line);

}
close($fileH);

# deletes files from previous run used to generate csv
unlink("thisIter.txt");
unlink("uniq.list");
unlink("heer.list");

foreach (@listOfIndivs){

	@newAry=split(" ",$_);
	$strName=$newAry[0]."_".$newAry[1]."_".$timestamp;
	

	open($fileH2,">","$strName.keep") or die("keep not created");
	print $fileH2 $_;
	# this and the above files create a plink file to analyze where a single individual comes off the test file list and is run against reference
	 `./plink --bfile $fileName --keep $strName.keep --extract plink.snplist --make-bed --out $strName`;
	 	# this filters the SNPs of the reference population on the test individual
	 	 `./plink --bfile $fileName --keep $strName.keep --extract plink.snplist --make-bed --out ThisIsATest`;

		$strName2= $strName."Test";
		$strName2Miss=$strName2."-merge.missnp";
	 	if(-e "$strName2Miss") 
		{
			    unlink("$strName2Miss");
		}
	 	`./plink --bfile REFERENCEFILES/$Refpop --bmerge $strName.bed $strName.bim $strName.fam --make-bed --out $strName2`;
		# the above linke tries to merge. the below conditional checks to see if an error log was thrown by plink
		# if the error log exists it tries to flip the 
			if (-e "$strName2Miss") 
			{ 
				system("mv $strName2Miss flip");
				$flip1 = "./plink --bfile $strName --flip  flip --make-bed --out $strName";
				system($flip1);
				# now it tries to merge again with the SNPs flipped (strands of DNA can be read two ways)
				 `./plink --bfile REFERENCEFILES/$Refpop --bmerge $strName.bed $strName.bim $strName.fam --make-bed --out $strName2`;
			    #file exists - set variable and other things
			}

			# some times it still causes issues. so you can just remove all the problem SNPs from the files to be merged

			if (-e "$strName2Miss")
			{ 
				$exclude = "./plink --bfile $strName --exclude  $strName2Miss --make-bed --out $strName ";
				system($exclude);
				 `./plink --bfile REFERENCEFILES/$Refpop --bmerge $strName.bed $strName.bim $strName.fam --make-bed --out $strName2`;
			    #file exists - set variable and other things
			    # this should fix everything with the merging assumine that the SNP is coded the same
			}

			#unlink("$strName2Miss");
			#unlink("flip");	

		# assume the test file was created, rename it to the normal name
		`./plink --bfile $strName2 --make-bed --out $strName`;

	# this parses the file to generate a header with the population names
	`sed '/AA_Ref/! s/^.*-9/_/g'  $strName.fam > $strName.pop.txt `;
	`cut -d' '  -f1 $strName.pop.txt > $strName.pop `;
	`sed -i.bak 's/AA_Ref_//g'  $strName.pop `;
	`uniq  $strName.pop | grep '[a-zA-Z]' > $timestamp.uniq.list`;
	$count = `wc -l < $timestamp.uniq.list`;
	($K) = $count =~ /(\d+)/;
	`tr -s '\n'  ' '< $timestamp.uniq.list > $timestamp.heer.list `;
	`sed -i.bak '1s/^/Group ID /'  $timestamp.heer.list `;
	`echo "" >>  $timestamp.heer.list `;

	# run admixture with cores
	 `./admixture  -j$cores  -s time --supervised $strName.bed $K `;
	 `cut -d' ' -f1,2 $strName.fam > $timestamp.colLabels.txt`;
	`paste -d ' ' $timestamp.colLabels.txt $strName.$K.Q | grep -v "AA_Ref" >> $timestamp.thisIter.txt`;
	  `cat $timestamp.thisIter.txt >> $timestamp.heer.list `;


	  # ok, now create the .csv file
	  $outName=$fileName."_".$Refpop.".".$timestamp.".csv";
	`cp $timestamp.heer.list $outName`;
	# remove a bunch of files 
	unlink($strName.".keep");
	unlink($strName.".bed");
	unlink($strName.".fam");
	unlink($strName.".bim");
	unlink($strName.".log");
	unlink($strName.".nosex");
	unlink($strName.".pop");
	unlink($strName.".pop.bak");
	`rm This*`;
	unlink($strName.".pop.txt");
	unlink($strName.".".$K.".Q");
	unlink($strName.".".$K.".P");

	close($fileH2);

}
	unlink($timestamp.".thisIter.txt");
	unlink($timestamp.".heer.list");
	unlink($timestamp.".heer.list.bak");
	unlink($timestamp.".colLabels.list");
		unlink($timestamp.".colLabels.txt");

	unlink($timestamp.".uniq.list");
