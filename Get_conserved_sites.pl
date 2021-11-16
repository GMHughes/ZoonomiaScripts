$file=$ARGV[0];
if($file=~m/ENST[0-9A-Z]+\.([A-Za-z0-9]+)\./){#code assumes a specific naming format used in TOGA output files
    $query_gene=$1;#take gene name from file
    print "Gene:".$query_gene."\n";
}
open(IN, "$file");#open in file to memory
@seqs;
$len=0;
$seqlenraw=0;
while(<IN>){#loop through each line in input file
    if($_=~m/>/){#Assuming MSA is in fasta format, if the line is the fasta header for sequence, do nothing
    }
    else{
	if($_=~m/[\S]+/){#otherwise, if line contains data, it is assumed to be sequence alignment data
	    $_=~s/\s//g;#remove whitespaces
	    $len=length($_);#get length of genes in alignment, assumes sequence is one full line!!!
	    push @seqs, $_;#add sequences onto an array
	    $raw=$_;#variable for the raw unaligned seq
	    $raw=~s/\-//g;#remove indels to unalign
	    $seqlenraw=length($raw);#get the length of the raw sequence
	    $avgraw+=$seqlenraw;#add it for getting average length (not sure this is used downstream)
	}
    }
}
$total_conserved=0;
$total_conserved_gaps=0;
$taxanum=scalar(@seqs)-1;#get number of taxa, but subtract 1 to account for the reference not being there!
open(TMP, "$file");#open filr again
@tmparray=<TMP>;
@tmparray=split(/\>/,join('',@tmparray));#read each aligned sequence as an array entry
close TMP;
for($k=0;$k<scalar(@tmparray);$k++){#loop through all sequences
    if($tmparray[$k]=~m/REFERENCE\n([\S]+)/){#find one that is human (REFERENCE)
	$refseq=$1;#reference sequence taken as first sequence
    }
}
$refseq=~s/\s//g;#remove any white characters
$rawseq=$refseq;
$rawseq=~s/\-//g;
$rawlen=length($rawseq);
$ungapped_num=0;
@conserved_sites=();
for($i=0;$i<$len;$i++){#for each character in the human sequence alignment
    $sitematch=0;
    $sitegapmatch=0;
    $gaps=0;
    $refchar=substr($refseq,$i,1);#take each character position, take character into memory
    if($refchar eq "-"){
    }
    else{
	$ungapped_num++;
    }
    for($j=1;$j<scalar(@seqs);$j++){#loop through each subsequent sequence
	$targettaxa=$seqs[$j];#take sequence into memory
	$targettaxa=~s/\s//g;#remove whitespace
	$targetchar=substr($targettaxa,$i,1);#take character in the sequence at same position
	if($targetchar eq "-"){#if its a gap, add to gap count
	    $gaps++;
	}
	if($targetchar eq $refchar || $targetchar eq "-" || $targetchar eq "X"){#if the character equals the same character as in the reference sequence
	    if($targetchar=~m/[\S]+/ || $targetchar eq "-" || $targetchar eq "X"){#make sure the character is alphanumerical  ###note just added the gap allowance here
		$sitematch++;#add one to number of matching istes
	    }
	}
	elsif($targetchar eq "-"){#if other taxa has a gap
	    $sitegapmatch++;#increase site gap count
	}
    }
    if($sitematch eq $taxanum){#if number of matching characters equals taxa numbers, it is conserved
	if($refchar ne "X" && $refchar ne "-"){
	    $total_conserved++;
	    $matchsite=$i+1;
	    $data="Aln_Site_".$matchsite."|real_site_".$ungapped_num.":".$refchar;
	    push @conserved_sites,$data;
	}
    }
    elsif($sitegapmatch eq $taxanum){#if number of gaps equals taxa, then is a gap across all, probably not gonna be a thing!!
	$total_conserved_gaps++;
    }
    elsif(($taxanum-$sitematch) ne $gaps && ($taxanum-$gaps)>1){
    }
}
my $minumum=$total_conserved;#minimum conserved sites is sites where there is no difference!
my $maximum=$total_conserved + $total_conserved_gaps;# maximum is the total with no difference plus ones where there is a difference due to a gap in one sequence
#print "Minimum conserved sites acoss taxa: ".$minumum."\/".$rawlen."\n";#$total_conserved."\/".$len."\n";
$avglen=$avgraw/scalar(@seqs);
$minpercseq=sprintf("%.2f",(($minumum/$avglen)*100));

$maxpercseq=sprintf("%.2f",(($maximum/$len)*100));
###print "Minimum conserved sites acoss taxa: ".$minumum."\/".$avglen.", in unaligned sequencs, using avg length: ".$minpercseq."% \n";
###print "Maximum conserved sites across taxa: ".$maximum."\/".$len.": ".$maxpercseq."\n";#$total_conserved_gaps."\/".$len."\n";

$maxpercunaligned=sprintf("%.2f",(($maximum/$avglen)*100));
###print "Maximum conserved sites across taxa: ".$maximum."\/".$avglen.", in unaligned sequencs, using avg length: ".$maxpercunaligned."% \n";

foreach $y(@conserved_sites){
    if($y=~m/real_site_([0-9]+)\:([A-Z]+)/){
	$check=0;
	$loc=$1;$aa=$2;
	$aaname=&name($aa);
	print "\n\n=>".$y;
	$varfile="/home/people/ghughes/scratch/ZoonomiaGenes/codon_ali_portion_2/".$query_gene."\.cv";
	open(IN2, "$varfile");
	while(<IN2>){
	    if($_=~m/\(([\S]+)\)\:c/){
		$targethit=$1;
	    }
	    if($targethit eq $query_gene){
		if($_=~m/[\s]+Pathogenic[\s]+/ || $_=~m/Likely[\s]+pathogenic/){
		    if($_=~m/\(p\.$aaname$loc[A-Z]+/){
			print "\n";
			print "ClinVar_entry: ".$_;
			$check++;
		    }
		}
	    }
	}
    }
    if($check eq "0"){
	print " None\n";
    }
}

sub name{
    $query_aa=$_[0];
    my (%genetic_code) = ('A' => 'Ala', 'R' => 'Arg', 'N' =>'Asn', 'D' => 'Asp', 'B' => 'Asx', 'C' => 'Cys', 'E' => 'Glu', 'Q' => 'Gln', 'Z' => 'Glx', 'G' => 'Gly', 'H' => 'His', 'I' => 'Ile', 'L' => 'Leu', 'K' => 'Lys', 'M' => 'Met', 'F' => 'Phe', 'P' => 'Pro', 'S' => 'Ser', 'T' => 'Thr', 'W' => 'Trp', 'Y' => 'Tyr', 'V' => 'Val'); 
    if(exists($genetic_code{$query_aa})){
	return $genetic_code{$query_aa};
    }
}
