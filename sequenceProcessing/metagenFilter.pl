#!/usr/bin/perl -w

#metagenomic filtering stats


use strict;
use File::Basename;
use Cwd;
use Getopt::Long;

###############################################################################################
# These are pointers to relevant directories. Modify before running script #
###############################################################################################

my $progDir = "path to program directory";
my $krakenDir = "path to kraken directory";
my $krakenDB = "path to kraken DB";
my $dataDir = "path to data directory";
my $sam = "path to samtools";


my $pwd = cwd();
my $reference = "";
my $help;
my $gtf;
my $remove;


GetOptions ( "reference=s" => \$reference,
			"gtf=s" => \$gtf,
			"remove=s" => \$remove,
			"help" => \$help);



my @files = @ARGV; #list of bam files. sample name should be first with "." separater. 

&help unless scalar @files > 0;
&help if $help;

#open (STATS, ">>", "kraken_stats_020919.txt") or die "couldn't open file for stats output: $?\n";
#print STATS "sample\ttotalSeqs\tnumberClassified\ttotalMTBC\tfinalMTBC\n";

foreach my $file (@files) {
	#start working on bam file, convert to fastq
	my $name = basename($file); 
	#print STDERR "$name\n";
	my ($sample, $temp1, $temp2) = split /\./, $name;
	my $fqOut = "${sample}.realn.fastq";
	my $filteredBam = "${sample}_classSeqs.bam";
	my $outVCF = "${sample}.vcf";
	my $mpileup1 = "${sample}_classSeqs.mpileup";
	my $indelGTF = "${sample}_indelRegs.gtf";
	my $mpileup2 = "${sample}_classSeqs_noIndel.mpileup";
	my $mpileup3 = "${sample}_classSeqs_noIndel_remRegs.mpileup";
	my $mpileup4 = "${sample}_classSeqs_noIndel_remRegs_noSB.mpileup";

        open (STATS, ">>", "${sample}_kraken_stats_020919.txt") or die "couldn't open file for stats output: $?\n";
        print STATS "sample\ttotalSeqs\tnumberClassified\ttotalMTBC\tfinalMTBC\n";

	mkdir "$sample", 0755 or warn "couldn't make directory $sample: $!";
	chdir "$sample" or die "can't change to directory $sample: $!";
	my $bwd = cwd();
	
	print STDERR "Processing $sample [bam2fq]...\n";
	system "$sam bam2fq ${pwd}/$file > ${bwd}/$fqOut";
	#system "gzip ${bwd}/$fqOut"; #can compress here, use other kraken option

#call kraken on fastqs
	
	#use below if inputting compressed fq
	#system "${progDir}/kraken/kraken -preload --threads 8 --fastq-input --gzip-compressed --db ${progDir}/kraken/DB/ ${bwd}/${fqOut}.gz > ${bwd}/${sample}.kraken";
	print STDERR "Processing $sample [kraken]...\n";
	system "${krakenDir}/kraken/kraken --threads 8 --fastq-input --db $krakenDB ${fqOut} > ${sample}.kraken";
#translate kraken file
	print STDERR "Processing $sample [kraken translate]...\n";
	system "${krakenDir}/kraken/kraken-translate --db $krakenDB ${bwd}/${sample}.kraken > ${bwd}/${sample}.kraken.labels";

#work on .kraken files
	open (IN, "<", "${bwd}/${sample}.kraken") or die "couldn't open $file: $?\n";
	print STDERR "Processing $sample [kraken output]...\n";
	my $numClass = 0;
	my $numUnclass = 0;
	# my $lengthsClass = 0;
	# my $lengthsUnclass = 0;
	# my $totalLength = 0;

	while (my $line = <IN>) {
		chomp $line;
		my ($class, $name, $tax, $length, $taxChain) = split "\t", $line;
		if ($class eq "C") {
			$numClass++;
			#$lengthsClass += $length;
		} elsif ($class eq "U") {
			$numUnclass++;
			#$lengthsUnclass += $length;
		} else {
			warn "found an unexpected value in 1st column: ${sample}.kraken";
			next;
		}
		#$totalLength += $length;
	}
	my $totalSeqs = $numClass + $numUnclass;
	# my $fracClass = $numClass / $totalSeqs;
	# my $avgLength = $totalLength / $totalSeqs;
	# my $avgLengthC = $lengthsClass / $numClass;
	# my $avgLengthU = $lengthsUnclass / $numUnclass || "no seqs";
	#print STDOUT "Total sequences: $totalSeqs\nFraction classified: $fracClass\nAverage Length total: $avgLength\nAverage Length C: $avgLengthC\nAverage Length U: $avgLengthU\n";
	close IN;
	
#now work on output from kraken-translate
	print STDERR "Processing $sample [kraken labels]...\n";
	open (LABELS, "<", "${bwd}/${sample}.kraken.labels") or die "couldn't open labels file: $?\n";
	#open (ALLNAMES, ">>", "${bwd}/${sample}_MTBC_allSeqs.txt") or die "couldn't open file for all the seqs: $?\n";
	my %seqNames;
	my $totalMTBC = 0;
	while (my $line = <LABELS>) {
		chomp $line;
		my ($seq, $taxList) = split /\t/, $line;
	#now filter out non MTBC hits
		my @taxa = split /;/, $taxList;
		if (scalar @taxa < 9) {
			next;
		} elsif ($taxa[7] ne "Mycobacterium") {
			next;
		} else {
		#	print ALLNAMES "$seqName\n";
		 	$totalMTBC++;
		# }
	#now count the number if times each sequence appears, removing the trailing /[12]
		my $seqName = (split /\//, $seq)[0];
		$seqNames{$seqName}++;
		}
	}
#now filter against sequences where only 1 of the pair was assigned
	open (NAMES, ">>", "${sample}_classSeq_names.txt") or die "couldn't open output for sequence names: $?\n";
	my $numMTBC = 0;
	foreach my $seqName (keys %seqNames) {
		if ($seqNames{$seqName} != 2) {
			next;
		} else {
			print NAMES "$seqName\n";
			$numMTBC++;
		}
	}

#now print some stats
	my $finalMTBC = 2 * $numMTBC;
	#print STDERR "Single end seqs: $diff\n";
	print STATS "$sample\t$totalSeqs\t$numClass\t$totalMTBC\t$finalMTBC\n" or die "couldn't print stats...\n";

#now filter the original bam file by sequence names
	print STDERR "Processing $sample [Picard FilterSamReads]...\n";
	system "java -Xmx8g -jar ${progDir}/picard-tools-1.138/picard.jar FilterSamReads INPUT=${pwd}/$file OUTPUT=${bwd}/$filteredBam READ_LIST_FILE=${bwd}/${sample}_classSeq_names.txt FILTER=includeReadList";

#now convert filtered bam to mpileup and do some further filtering on the mpileup
	print STDERR "Processing $sample [mpileup generation and filtering]...\n";
	#open (INBAM, "<", "${sample}_MTBCseqs.bam") or die "couldn't open BAM to generate mpileup: $?\n";
	system "$sam mpileup -B -Q 20 -f ${pwd}/$reference ${bwd}/$filteredBam > ${bwd}/$mpileup1";
	system "perl ${progDir}/popoolation_1.2.2/basic-pipeline/identify-genomic-indel-regions.pl --input ${bwd}/$mpileup1 --output ${bwd}/$indelGTF";
	system "perl ${progDir}/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input ${bwd}/$mpileup1 --gtf ${bwd}/$indelGTF --output ${bwd}/$mpileup2";
	system "perl ${progDir}/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input ${bwd}/$mpileup2 --gtf ${pwd}/$remove --output ${bwd}/$mpileup3";


#strand bias filter
	print STDERR "Processing $sample [generateing VCF]...\n";
	system "$sam mpileup -B -Q 20 -f ${pwd}/$reference -uv -t DP,DP4,SP ${bwd}/$filteredBam > ${bwd}/$outVCF";
	print STDERR "Processing $sample [strand bias filter]...\n";
	open (VCF, "<", "${bwd}/$outVCF") or die "couldn't open $outVCF: $?\n";
	#open (SBOUT, ">", "${bwd}/${sample}_SB_positions.txt") or die "couldn't open output file for strand bias positions: $?\n";
	my %SBpos;
	while (my $line = <VCF>) {
    	if ($line =~ /^#/) {
			next;
    	}
    	my ($pos, $field) = (split /\t/, $line)[1,9];
    	my $val = (split /:/, $field)[2];
    	if ($val >= 13) {
			$SBpos{$pos}=$val;
    	} else {
    		next;
    	}
    }
    open (MPILEUP, "<", "${bwd}/$mpileup3") or die "couldn't open mpileup for strand bias filtering: $?\n";
    open (MPOUT, ">>", "${bwd}/$mpileup4") or die "couldn't open file for strand bias output: $?\n";
    while (my $line = <MPILEUP>) {
    	my $position = (split /\t/, $line)[1];
    	if ($SBpos{$position}) {
    		next;
    	} else {
    		print MPOUT "$line"
    	}
    }
#now call automatePoolSeqDiversityStats.py script
	print STDERR "Processing $sample [popoolation]...\n";
	system "${pwd}/automatePoolSeqDiversityStats.py -p ${sample} -g ${pwd}/${gtf}";

#remove all but the final mpileup (.noIndel.remRegs.mpileup)
#	unlink "$mpileup1";
#	unlink "$mpileup2";
#	unlink "${mpileup2}.params";
#	unlink "$fqOut";
#	unlink "${sample}.kraken.label";
#	unlink "$filteredBam";
#	unlink "$mpileup3";
#	unlink "${mpileup3}.params";
#	unlink "${sample}_classSeqs.reads";
#	unlink "$indelGTF";
#	unlink "${indelGTF}.params";
#	unlink "${sample}.kraken";
#	unlink "${sample}.kraken.labels";
#	unlink "${sample}_classSeq_names.txt";
#	unlink "${sample}.realn.reads";
#	unlink "$outVCF";
#get averages of subsample
	print STDERR "Processing $sample [getting subsample averages]...\n";
	system "${pwd}/get_subsampleAverage.pl -output ${bwd}/${sample}_subAvg.pi ${bwd}/${sample}*rand*n10K.pi";
	system "${pwd}/get_subsampleAverage.pl -output ${bwd}/${sample}_subAvg.theta ${bwd}/${sample}*rand*n10K.theta";
	system "${pwd}/get_subsampleAverage.pl -output ${bwd}/${sample}_subAvg.td ${bwd}/${sample}*rand*n10K.td";
	system "${pwd}/get_subsampleAverage_genes.pl -output ${bwd}/${sample}_subAvg_genes.pi ${bwd}/${sample}*rand*gene.pi";
	system "${pwd}/get_subsampleAverage_genes.pl -output ${bwd}/${sample}_subAvg_genes.theta ${bwd}/${sample}*rand*gene.theta";
	chdir "$pwd";
}

sub help {
    print"
[usage]
./metagenFilter.pl [-options] <bam1> <bam2> ...

[description]
This script performs various filtering on aligned sequences in a bam file, namely running sequences through Kraken and kicking out sequences that aren't classified to a specific taxon. After filtering, it generates an mpileup and runs this through PoPoolation.

[output files]
kraken_stats.txt             contains some stats on the number of classified sequences at various steps.
prefix_classSeqs_noIndel_remRegs_noSB.mpileup      final filtered mpileup. This file is fed to PoPoolation.
prefix_rand[0..9]  various subsampled PoPoolation output files
prefix_subAvg	various genomewide averages calculated from all subsampled files

[options]
-reference      specify a reference in fasta format
-gtf 			provide a gtf file containing CDS information
-remove			gtf file of regions that should be removed
";
exit;
}

