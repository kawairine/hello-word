#input is a list of PDB codes
#foreach code we download the transformed file from OPM, assign the Nterminal, keep the sequence and construct topology (if TM annotation is available)
#we keep only the 1st amino-acid for each chain and check its Z-coordinate
#if Z<0 -> Nterm: in, else Nterm: out

#perl grep_topo_OPM <LIST_OF_PDB_CODES_NO_CHAINS>

open NTERM_INFO, ">$ARGV[0].only_Nterminal";
open TOPO_INFO, ">$ARGV[0].full_topology";

use List::Util qw/min max/;
`rm *pdb`;
`rm protein*`;

my %three2one = (	ALA=>'A',TYR=>'Y',MET=>'M',LEU=>'L',CYS=>'C',GLY=>'G',
       			 	ARG=>'R',ASN=>'N',ASP=>'D',GLN=>'Q',GLU=>'E',HIS=>'H',TRP=>'W',
       				LYS=>'K',PHE=>'F',PRO=>'P',SER=>'S',THR=>'T',ILE=>'I',VAL=>'V'
       			);

my $opm_list_codes=$ARGV[0];
open LIST, $opm_list_codes;
while(<LIST>)
{
	my $opm_code=$_;
	chomp $opm_code;

	#download the respective PDB file from OPM (transformed)
	my $pdb_file_link="http://opm.phar.umich.edu/pdb/".$opm_code.".pdb";
  `wget $pdb_file_link >/dev/null 2>&1`;

	#download the respective HTML file to get the TM boundaries
	my $html_file_link="http://opm.phar.umich.edu/protein.php?search=".$opm_code;
	`wget –quiet $html_file_link >/dev/null 2>&1`;
	my $downloaded_filename="protein.php?search=".$opm_code;	#this is how the file is saved after download

	my %HoA_Nterminal=();	#to keep the information for the first amino-acid in OPM sequence (aka the Nterminal)
	my %HoA_sequences=();  	#to construct the sequence in the OPM entries
	my $Nterminal_OPM='';	#the assigned N-terminal based on OPM annotation
	my %sequences_all='';	#hash of hashes to keep all sequences with each position and AA (in the case of multi-chain examples)

	open OPM, "$opm_code.pdb";
	while(<OPM>)
	{
		if($_=~/^ATOM\s+\d+\s+N\s+(\w{3})\s+(\w{1}[\s\d]+)\s+[\d\.\-]+\s+[\d\.\-]+\s+([\d\-\.]+).*/)
		{
			my $amino_acid=$1;
				my $aa=$three2one{$amino_acid};
			my $chain_pos=$2;
				my $chain=substr($chain_pos,0,1);	#keep only the chain
				my $position=substr($chain_pos,1);	#keep the rest (position and possible gaps (pos<1000))

			my $Zcoordinate=$3;
			if($position=~/\s+/) {$position=~s/\s+//g;}	#this is for cases like A1024 or A 765 (they may or may not include gaps between chain and position)

			$sequences_all{$chain}{$position} = $aa;

			#this part is only to keep the smallest position (aka Nterminal)
			if (not exists $HoA_Nterminal{$chain} or $position < $HoA_Nterminal{$chain}[0])
			{
				$HoA_Nterminal{$chain} = [$position, $Zcoordinate];
			}
		}
	}
	close OPM;

	#this block of code loops through all chains for a given PDB code, checks each position in the respective sequence
	#from the OPM-PDB file, and, if there are some positions between the lowest and the maximum (start and end) as denoted
	#in the OPM file, it fills them with '-' (like in entry 1xio)
	for my $chain_separate1 ( keys %sequences_all )
	{
		if($chain_separate1 ne '')
		{
		    my ($min, $max) = (min(keys %{ $sequences_all{$chain_separate1} }), max(keys %{ $sequences_all{$chain_separate1} }));
		    for my $check_position ($min .. $max)
		    {
		    	unless (exists $sequences_all{$chain_separate1}{$check_position})
		    	{
		    		$sequences_all{$chain_separate1}{$check_position} = '-';
		    	}
		    }
		}
	}

	for my $chain_separate ( keys %sequences_all )
	{
		my $opm_seq='';			#the sequence in the OPM entry
		if($chain_separate ne '')
		{
	   		$fullid=$opm_code.'_'.$chain_separate;

	   		for my $seq_pos( sort {$a<=>$b} keys %{ $sequences_all{$chain_separate} } )
	   		{
	   			$opm_seq.= $sequences_all{$chain_separate}{$seq_pos};
	   		}
			my $seq_len=length($opm_seq);

			my $starting_position=$HoA_Nterminal{$chain_separate}[0]-1;			#find the position that counting starts for OPM (we need that for TM-parts) [-1 because counting starts at 0]
			if($HoA_Nterminal{$chain_separate}[1]>0) 	{$Nterminal_OPM='o';}	#assign the Nterminal for OPM topology based on the Z-coordinate
			elsif($HoA_Nterminal{$chain_separate}[1]<0) {$Nterminal_OPM='i';}

			#for each chain, get the annotation for the TM parts
			my $TM_all1='';			#all TM segments information initial
			my $TM_all='';			#all TM segments information for topology construction
			open HTML_OPM, $downloaded_filename;
			while(<HTML_OPM>)
			{
				if($_=~/\<td align\=\"left\"\>\<b\>$chain_separate\<\/b\>\s+\-\s+Tilt\:\s+\d+\&\#\d+\;\s+\-\s+Segments\:\s+(.*?)\<\/td\>/)
				{
					$TM_all1=$1;
				}
			}
			close HTML_OPM;

			if($TM_all1 ne '')
			{
				my @tms_opm=split(/\,/, $TM_all1);
				for (my $k=0; $k<=$#tms_opm; $k++)
				{
					if($tms_opm[$k]=~/\d+\((.*?)\)/)
					{
						my $actual_tm_opm=$1;
						$actual_tm_opm=~s/\s+//g;
						$TM_all.=$actual_tm_opm.',';
					}
				}
			}
			chop $TM_all;

			if($TM_all ne '')
			{
				print TOPO_INFO ">$fullid|".($starting_position+1)."|Z(".($starting_position+1).")=$HoA_Nterminal{$chain_separate}[1]|$Nterminal_OPM|$TM_all\n$opm_seq\n";

				my %opposite= (i=>'o', o =>'i');
				######################################### for results with just 1 TM helix predicted #######################################################################################
				if ($TM_all!~/\,/)
				{
					my $type_first=$Nterminal_OPM;
					my $type_end=$opposite{$Nterminal_OPM};
					my @single_tm_split=split("-", $TM_all);
					my $single_tm_start=$single_tm_split[0]-$starting_position-1;
					my $single_tm_end=$single_tm_split[1]-$starting_position-1;
					my $single_len_tm=$single_tm_end-$single_tm_start+1;

					my $length_initial=$single_tm_start;
					print TOPO_INFO $type_first x $length_initial;
					print TOPO_INFO 'M' x $single_len_tm;
					my $start_final=$single_tm_end+1;
					my $length_final=$seq_len-$start_final;
					print TOPO_INFO $type_end x $length_final;
				}

				######################################### for results with more than 1 TM helices predicted ################################################################################
				else
				{
					my @AoA_tm=();
					my @split_tm=split(",", $TM_all);
					for(my $j=0; $j<=$#split_tm; $j++)
					{
						my $tm=$split_tm[$j];
						my @split_tm_part=split("-", $tm);
						my $tm_start=$split_tm_part[0]-$starting_position-1;
						my $tm_end=$split_tm_part[1]-$starting_position-1;
						push @AoA_tm, [$j, $tm_start, $tm_end];
					}
					my $how_many_tm=@split_tm;
					my $real_parts=$how_many_tm-1;
					my $type_first='';
					my $type_for_previous='';

					for my $i ( 0 .. $#AoA_tm )
					{
						#case 1 -> initial part
						if($i == 0)
						{
							#construct previous label part, which will be same as N-terminal x length before start of first TM
							$type_first=$Nterminal_OPM;
							my $length_initial=$AoA_tm[$i][1];
							print TOPO_INFO $type_first x $length_initial;

							#now print the TM that follows the initial part of the label
							my $start_tm=$AoA_tm[$i][1];
							my $end_tm=$AoA_tm[$i][2];
							my $length_tm=$end_tm-$start_tm+1;
							print TOPO_INFO 'M' x $length_tm;
						}

						#case 2 -> intermediate part
						elsif($i>=1)
						{
							#$type_first remains, so it can be used for the first time (i=1) and then become ''
							if ($type_first ne '') 	{$type_for_previous=$opposite{$type_first};}
							else		       		{$type_for_previous=$opposite{$type_for_previous};}
							#print the label part before the TM
							my $start_between=$AoA_tm[$i-1][2]+1;
							my $end_between=$AoA_tm[$i][1]-1;
							my $length_between=$end_between-$start_between+1;
							print TOPO_INFO $type_for_previous x $length_between;

							# print the TM that is inbetween the two label parts
							my $start_tm=$AoA_tm[$i][1];
							my $end_tm=$AoA_tm[$i][2];
							my $length_tm=$end_tm-$start_tm+1;
							print TOPO_INFO 'M' x $length_tm;

							$type_first='';
							$type_for_previous=$type_for_previous;

							# now construct final label part, which will be the opposite of what was the part before the last TM x length (total_seq_length-end of final tm)
							if ($i==$real_parts)
							{
								my $type_last=$opposite{$type_for_previous};
								my $length_last=$seq_len-$AoA_tm[$i][2]-1;

								print TOPO_INFO $type_last x $length_last;
							}
						}
					}
				}
				print TOPO_INFO "\n";
			}
			else
			{
				print NTERM_INFO ">$fullid|".($starting_position+1)."|Z(".($starting_position+1).")=$HoA_Nterminal{$chain_separate}[1]|$Nterminal_OPM|$TM_all\n$opm_seq\n";
			}
		}
	}
}
close LIST;
close NTERM_INFO;
close TOPO_INFO;
