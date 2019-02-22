#!/usr/bin/env perl
#/gsc/btl/linuxbrew/bin/perl

#AUTHOR
#   Rene Warren, Hamid Mohamadi, Jessica Zhang, Lauren Coombe
#   rwarren at bcgsc.ca
#   October/November 2018

#NAME
#   ntEdit prototype / formerly haplotor

#SYNOPSIS

#DOCUMENTATION
#   README.md distributed with this software @ www.bcgsc.ca
#   We hope this code is useful to you -- Please send comments & suggestions to rwarren * bcgsc.ca
#   If you use ntedit, its code or ideas, please cite our work

#LICENSE
#ntedit Copyright (c) 2018-2019 British Columbia Cancer Agency Branch.  All rights reserved.
#ntedit is released under the GNU General Public License v3

use strict;
use POSIX;
use FindBin;
use lib "$FindBin::Bin/./lib";
use BloomFilter;
use Time::HiRes;
use Getopt::Std;
use vars qw($opt_f $opt_d $opt_k $opt_v $opt_b $opt_r $opt_z $opt_i $opt_x $opt_y);
getopts('f:r:k:d:i:z:b:v:x:y:');

#### SET DEFAULT VARIABLES ####
my ($bf_file,$base_name,$k,$min_size,$verbose)=("","",35,100,0);
my $i = 3;###factor1 low helps prevent FP.  factor2 high help make changes that are confirmed by subset of k kmers
my ($MAXBASEINS,$MAXBASEDEL,$MAXSTRINGLEN) = (0,0,5000000);
my ($factor1,$factor2) = (5,9);### was 5,9 for ecoli,celegans,human
my $version = "[v1.0.1]";

#-------------------------------------------------

if(! $opt_f || ! $opt_r){
   print "\nUsage: $0 $version\n";
   print "-f  draft genome assembly (Multi-FASTA format, required)\n"; 
   print "-r  Bloom filter of sequence reads (ntHits format, required)\n";
   print "-k  k-mer value (default -k $k, optional, same value of k to build -r Bloom filter)\n";
   print "-d  maximum number of base deletions (default -d $MAXBASEDEL, optional, range 0-5)\n";
   print "-i  maximum number of base insertions (default -i $MAXBASEINS, optional, range 0-5, values higher than 4 will impact run speed)\n";
   print "-x  leniency factor 1 : determines whether a missing kmer should be edited (default -x $factor1, optional, lower value=less permissive)\n";
   print "-y  leniency factor 2 : determines whether a change should be kept (default -y $factor2, optional, lower value=less permissive)\n";
   print "-z  minimum contig length to consider (default -z $min_size, optional)\n";
   print "-b  base name for your output files (optional)\n";
   print "-v  Runs in verbose mode (-v 1 = yes, default = no, optional)\n";
   print "\nNOTE: BLOOM FILTER MUST BE DERIVED FROM THE SAME SEQUENCE READS FILE(S) SUPPLIED to ntHits, WITH SAME -k VALUE SPECIFIED\n";
   die "\nError: Missing mandatory options -f, -r\n\n";
}

################## ASSIGN VARIABLES ###################

my $assemblyfile = $opt_f;
$k = $opt_k if($opt_k);
$verbose = $opt_v if($opt_v);
$min_size = $opt_z if($opt_z);
$base_name = $opt_b if($opt_b);
$bf_file = $opt_r if($opt_r);
$MAXBASEDEL = $opt_d if($opt_d);
$MAXBASEDEL = 5 if($MAXBASEDEL>5);
$MAXBASEINS = $opt_i if($opt_i);
$MAXBASEINS = 5 if($MAXBASEINS>5);
$factor1 = $opt_x if($opt_x);
$factor2 = $opt_y if($opt_y);

if($MAXBASEINS==0 || ($MAXBASEINS==1 && $MAXBASEDEL>1)){
   print "\n*WARNING: Because insertions are set at -i $MAXBASEINS, the deletions will also be set at -d $MAXBASEINS)\n";
   print "-" x 60, "\n";
   $MAXBASEDEL = $MAXBASEINS;
}

my ($checkkct,$subskct,$inskct,$delkct) = (($k/$factor1),($k/$factor2),($k/$factor2),($k/$factor2));
#eg.
#k=40. factor =10 ratio = 4
#k=40. factor =5  ratio = 8.  ctmissing would have to be higher to investigate or accept a change. Therefore, lowering factor may increase stingency

my $outofpossible = int($k/$i) + 1;

###create insert combinations
my $inserthash;
my @dnabases = ('A','C','G','T');
foreach my $initialbase(@dnabases){

   my @singlearr = ("$initialbase");
   my (@level1arr,@level2arr,@level3arr,@level4arr);
   foreach my $dnabase1 (@dnabases){
      push @level1arr, $initialbase . $dnabase1;  ### 4 possibilities for given initial base
      foreach my $dnabase2 (@dnabases){
         push @level2arr, $initialbase . $dnabase1 . $dnabase2; ###
         foreach my $dnabase3 (@dnabases){
            push @level3arr, $initialbase . $dnabase1 . $dnabase2 . $dnabase3;
            foreach my $dnabase4 (@dnabases){
               push @level4arr, $initialbase . $dnabase1 . $dnabase2 . $dnabase3;
            }
         }
      }
   }
   my @insertarr;                                                           ### 1 possibility for given base
   if($MAXBASEINS == 5){
      @insertarr = (@singlearr,@level1arr,@level2arr,@level3arr,@level4arr);### 1 + 4 + 16 + 64 + 256 = 341 possibilites for given base 
   }elsif($MAXBASEINS == 4){
      @insertarr = (@singlearr,@level1arr,@level2arr,@level3arr);           ### 1 + 4 + 16 + 64 = 85 possibilites for given base
   }elsif($MAXBASEINS == 3){
      @insertarr = (@singlearr,@level1arr,@level2arr);                      ### 1 + 4 + 16 = 19 possibilites for given base
   }elsif($MAXBASEINS == 2){
      @insertarr = (@singlearr,@level1arr);                                 ### 1 + 4  = 5 possibilities for given base
   }elsif($MAXBASEINS == 1){
      @insertarr = @singlearr; 
   }
   $inserthash->{$initialbase} = \@insertarr; ### hash or insert bases arrays, where all insertbase combos up to 
}

################## END ASSIGN VARIABLE ##################



my $assemblyruninfo="";


if(! -e $assemblyfile){
   die "Invalid file: $assemblyfile -- fatal\n";
}


### Naming output files
if ($base_name eq ""){

   $base_name = $assemblyfile . "_r-" . $bf_file . "_k" . $k . "_d". $MAXBASEDEL . "_i" . $MAXBASEINS . "_x" . $factor1 . "_y" . $factor2 . "_z" . $min_size;
   my $pid_num = getpgrp(0);
   $base_name .= "_pid" . $pid_num;
}

my $log = $base_name . ".log";
my $editeddraft_file = $base_name . "_edited.fa";
my $changes_file = $base_name . "_changes.tsv";

open (LOG, ">>$log") || die "Can't write to $log -- fatal\n";

#-------------------------------------------------
#
my $init_message = "\nRunning: $0 $version\n-f $assemblyfile\n";
$init_message .= "-k $k\n-d $MAXBASEDEL\n-i $MAXBASEINS\n-x $factor1\n-y $factor2\n-z $min_size\n-b $base_name\n-r $bf_file\n\n\n----------------- Verifying files -----------------\n\n";

print $init_message;
print LOG $init_message;
$assemblyruninfo=$init_message;

#-------------------------------------------------
my $file_message = "";

if(! -e $assemblyfile){
   $file_message = "\nInvalid file: $assemblyfile -- fatal\n";
   print $file_message;
   print LOG $file_message;
   exit;
}else{
   $file_message = "Checking sequence target file $assemblyfile...ok\n";
   print $file_message;
   print LOG $file_message;
}

#-------------------------------------------------
#PRGM STARTS
#-------------------------------------------------

my $date = `date`;
chomp($date);
my $bloom;

eval{

my $bffile_message="";
my $bfreusemessage = "A Bloom filter was supplied ($bf_file)\n";
print LOG $bfreusemessage;
print $bfreusemessage;

if(! -e $bf_file){
   $bffile_message = "\nInvalid file: $bf_file -- fatal\n";
   print $bffile_message;
   print LOG $bffile_message;
   exit;
}else{
   $bffile_message = "Checking Bloom filter file $bf_file...ok\n";
   print $bffile_message;
   print LOG $bffile_message;
}

$bffile_message="Loading bloom filter from $bf_file\n";
print $bffile_message;
print LOG $bffile_message;

$bloom = new BloomFilter::BloomFilter($bf_file);


#-------------------------------------------------
$date = `date`;
chomp($date);


my $reading_seqs_message = "\n\n=>Reading sequence contigs (to edit), tracking k-mer positions : $date\n";
print $reading_seqs_message;
print LOG $reading_seqs_message;
$assemblyruninfo.=$reading_seqs_message;
&readContigs($assemblyfile,$k,$min_size,$bloom,$editeddraft_file,$changes_file);

#-------------------------------------------------
$date = `date`;
chomp($date);

my $sc_start_message = "\n=>Process completed : $date\n\nThe new draft is: $editeddraft_file\nwith accompanying changes: $changes_file\n";
print $sc_start_message;
print LOG $sc_start_message;

};###end eval block

#-------------------------------------------------
$date = `date`;
chomp($date);

if($@){
   my $message = $@;
   my $failure = "\nSomething went wrong running $0 $date\n$message\n\n";
   print $failure;
   print LOG $failure;
   $assemblyruninfo.=$failure . "\n";
}else{
   my $success = "\n$0 executed normally $date\n\n";
   print $success;
   print LOG $success;
   $assemblyruninfo.=$success . "\n";
}

print "$0 $version terminated on $date\n\n";
print LOG "$0 $version terminated on $date\n\n";

close LOG;

exit;





#----------------
sub readContigs{
   my ($file,$k,$min_size,$bloom,$editeddraft_file,$changes_file) = @_;

   my $prevhead = "";
   my $seq = "";
   my $cttig=0;

   open(FA, ">$editeddraft_file") || die "Can't write to $editeddraft_file -- fatal.\n";
   open(CHG,">$changes_file") || die "Can't write to $changes_file -- fatal.\n";

   ###LIST ALL VARIATIONS AND CHANGES MADE
   print CHG "ID\tbpPosition+1\tOriginalBase\tNewBase\tSupport $k-mers (out of $outofpossible)\tAlternateNewBase\tAlt.Support $k-mers\n";

   ###Support for compressed files MAR2016
   if($file=~/zip$/i){
      open(IN,"unzip -p $file|") || die "Error reading $file -- fatal\n";
   }elsif($file=~/gz$/i || $file=~/gzip$/i){
      open(IN,"gunzip -c $file|") || die "Error reading $file -- fatal\n";
   }else{
      open(IN,$file) || die "Error reading $file -- fatal\n";
   }

   my $contigs_processed_message = "Contigs (>= $min_size bp) processed k=$k:\n";
   print $contigs_processed_message;
   print LOG $contigs_processed_message;
   ###
   while(<IN>){
      chomp;
      if(/^\>(.*)/){
         my $head=$1;

         if ($head ne $prevhead && $seq ne '' && $prevhead ne ''){
            $cttig++;
            print "\r$cttig";
            $|++;

            my $tiglen = length($seq);
            if($tiglen >= $min_size){
               my $track_all;
               ($track_all,$seq) = &kmerizeCheck(uc($seq),$track_all,$k,$prevhead);
               print FA ">$prevhead\n$seq\n";
               my $poslist = $track_all->{$prevhead};
               foreach my $coord(sort {$a<=>$b} keys %$poslist){
                  print CHG "$prevhead\t$coord\t$poslist->{$coord}{'draft'}\t$poslist->{$coord}{'change'}\t$poslist->{$coord}{'support'}\t$poslist->{$coord}{'altchgbase'}\t$poslist->{$coord}{'altchgsupport'}\n" if($poslist->{$coord}{'draft'} ne $poslist->{$coord}{'change'} || $poslist->{$coord}{'altchgbase'} ne "");
               }
            }else{
               print FA ">$prevhead\n$seq\n";
            }
         }
         $seq = '';
         $prevhead = $head;
      }else{
         $seq .= $_;
      }
   }
   $cttig++;
   print "\r$cttig";
   print LOG "$cttig\n\n";
   $|++;

   if(length($seq) >= $min_size){
      my $track_all;
      ($track_all,$seq) = &kmerizeCheck(uc($seq),$track_all,$k,$prevhead);
      print FA ">$prevhead\n$seq\n";
      my $poslist = $track_all->{$prevhead};
      foreach my $coord(sort {$a<=>$b} keys %$poslist){
         print CHG "$prevhead\t$coord\t$poslist->{$coord}{'draft'}\t$poslist->{$coord}{'change'}\t$poslist->{$coord}{'support'}\t$poslist->{$coord}{'altchgbase'}\t$poslist->{$coord}{'altchgsupport'}\n" if($poslist->{$coord}{'draft'} ne $poslist->{$coord}{'change'} || $poslist->{$coord}{'altchgbase'} ne "");
      }
   }else{
      print FA ">$prevhead\n$seq\n";
   }
   ###
   close IN;
   close FA;
   close CHG;

}

#----------------
sub kmerizeCheck{

   my ($masterseq,$track_all,$k,$head) = @_;

   ### sequence/sub-specific variables
   my $newmaster = "";
   my $last2k = "";
   my $adjustedpos = 0; 
   my $mastercoord = 0;
   my $chunknum = 0;
   #my $failsafe = 0;### not used at the moment, this led to abrupt terminations. Used for logs.

   print "\nprocessing chunks...\n";

   for(my $mpos=0; $mpos<= length($masterseq); $mpos+=$MAXSTRINGLEN){ 
      $chunknum++;
      my $seq = substr($masterseq,$mastercoord,$MAXSTRINGLEN);###This is needed for better memory management on large strings
      my $len = length($seq);
      print "Adjusted position: $adjustedpos Tile-mpos:$mpos, mastercoord: $mastercoord step: $MAXSTRINGLEN (length of seq = $len)\n";

      POSITION:
      for(my $pos=0;$pos<=(length($seq)-$k);$pos++){
         my $kmer = substr($seq,$pos,$k);
         print "\n$pos\t$kmer\n" if($verbose);

         if($kmer=~/^[ACGT]*$/ && ! $bloom->contains($kmer)){# the reads do not contain this kmer and only contain known bases
         #if($kmer=~/^[ACGTRYSWKMBDHV]*$/ && ! $bloom->contains($kmer)){# the reads do not contain this kmer and only contain known bases SUPPORT IUPAC
            ###check the next k kmers
            my $ctmissing=0;

            for(my $subset=$pos+1;$subset<=($pos+$k);$subset+=$i){###important to catch areas of no kmer coverage
               my $subkmer = substr($seq,$subset,$k);
               last POSITION if(length($subkmer)<$k);
               next POSITION if($subkmer=~/[EFIJLNOPQUXZ]/);  ###this prevents stalling when bases other than ATCG are encountered
               $ctmissing++ if(! $bloom->contains($subkmer)); ###expect all (most) missing if base is different in haploid set
            }

            if($ctmissing >= $checkkct){### this could be changed, depends on false pos rate I guess; bottom line is if base isn't in readset ctmissing==k, unless FP exists || the larger factor is, the smaller this fraction is and the least stringent is the requirement MOST ARE MISSING

               my $draftBase = substr($seq,($pos+$k-1),1);
               my $substitution = 0;
               my $coord = $pos + $k - 1;
               my $oneBasedCoord = $adjustedpos + $coord + 1; ### 1-base coordinate system

               my @basearr;### assignment of 3' bases
               if($draftBase eq "A"){ 
                 @basearr=('T','C','G');### seems strange to put this in conditional, but the base order in the array is important (the first base that is)
               }elsif($draftBase eq "T"){ 
                 @basearr=('A','C','G');    
               }elsif($draftBase eq "C"){ 
                 @basearr=('A','T','G');
               }elsif($draftBase eq "G"){ 
                 @basearr=('A','T','C');
               }else{
                 @basearr=('A','T','C','G');### could be used for IUPAC. Must uncomment line 356 and comment out 355. At moment, IUPAC ambiguity codes are ignored.
               }

               BASE:
               foreach my $changeBase(@basearr){

                  ###confirmation that the base is not in the haploid set
                  my $repkmer = substr($seq,$pos,($k-1)) . $changeBase;

                  print "\tmissing:$ctmissing chck:$checkkct  LOOKING AT BLOOM PRESENCE FOR $repkmer >$changeBase< because $kmer and k kmers missing\n" if($verbose);

                  #SUBSTITUTION
                  if($bloom->contains($repkmer)){
                     my $ctchg=0;
                     $_ = $seq;
                     ###This changes the base in the $seq string
                     substr($_, $pos+$k-1, 1) = $changeBase;#lc($changeBase);

                     ###check the next k kmers
                     for(my $subset=$pos+1;$subset<=($pos+$k);$subset+=$i){   #### this code prevents changing the seq string when base change caused by indels
                        my $subkmer = substr($_,$subset,$k);
                        $ctchg++ if($bloom->contains($subkmer)); #### expect very few to be missing if correction is valid
                     }

                     #### IF MET, CHECK SUBSTITUTIONS, ELSE GO TO INDELS.  SAME BASE WILL GO TO INDELS
                     if($ctchg >= $subskct){### this indicates that the change was likely valid. As factor increases, ratio decreases making it harder to meet condition
                        ###This changes the base in the original
                        if($track_all->{$head}{$oneBasedCoord}{'change'} ne ""){### a change was made before
                           if($ctchg > $track_all->{$head}{$oneBasedCoord}{'support'}){###current change better supported
                              $seq = $_; ### make change permanent
                              $track_all->{$head}{$oneBasedCoord}{'altchgbase'} = $track_all->{$head}{$oneBasedCoord}{'change'};
                              $track_all->{$head}{$oneBasedCoord}{'altchgsupport'} = $track_all->{$head}{$oneBasedCoord}{'support'};
                              $track_all->{$head}{$oneBasedCoord}{'draft'} = $draftBase;
                              $track_all->{$head}{$oneBasedCoord}{'change'} = $changeBase;
                              $track_all->{$head}{$oneBasedCoord}{'support'} = $ctchg;
                              print "\t\tSUBS: FIXED $draftBase\t$changeBase\t$head @ $oneBasedCoord" if($verbose);
                              $substitution++;
                              #$failsafe = 0;
                              next POSITION if($substitution==3);### if this is at 3 alternate base, means a change WAS MADE USING SUBS, NO POINT INSPECTING INDELS
                              next BASE;
                           }else{
                              $track_all->{$head}{$oneBasedCoord}{'altchgbase'} = $changeBase;
                              $track_all->{$head}{$oneBasedCoord}{'altchgsupport'} = $ctchg;
                              $substitution++;
                              next POSITION if($substitution==3);
                              next BASE;
                           }
                        }else{### no previous change
                              $seq = $_; ### make change permanent
                              $track_all->{$head}{$oneBasedCoord}{'draft'} = $draftBase;
                              $track_all->{$head}{$oneBasedCoord}{'change'} = $changeBase;
                              $track_all->{$head}{$oneBasedCoord}{'support'} = $ctchg;
                              print "\t\tSUBS: FIXED $draftBase\t$changeBase\t$head @ $oneBasedCoord" if($verbose);
                              $substitution++;
                              #$failsafe = 0;
                              next POSITION if($substitution==3);
                              next BASE;
                        }
                     }elsif(! $substitution){### NOT A SIMPLE SUBSTITUTION, CHECK INDEL(S)
                        print "\tSUBS not found, checking indels..\n" if($verbose);
                        ###INDEL BLOCK
                        my $ctelement = 0;
                        foreach my $insertnt (@{$inserthash->{$changeBase}}){
                           $ctelement++;
                           #### INSERTIONS
                           $_ = $seq;
                           substr($_, $pos+$k-1, 0) = $insertnt;### insert single

                           ###check the next k kmers
                           my $ctins=0;
                           for(my $subset=$pos+1;$subset<=($pos+$k-1);$subset+=$i){   #### this code prevents changing the seq string when base change caused by indels
                              my $subkmer = substr($_,$subset,$k);
                              next POSITION if($subkmer=~/[EFIJLNOPQUXZ]/);  ###this prevents stalling when bases other than ATCG are encountered
                              $ctins++ if($bloom->contains($subkmer));
                           }
                           if($ctins >= $inskct){###
                              ###This changes the base in the original
                              $seq = $_;
                              $track_all->{$head}{$oneBasedCoord}{'draft'} = $draftBase;
                              $track_all->{$head}{$oneBasedCoord}{'change'} = "+" . $insertnt;
                              $track_all->{$head}{$oneBasedCoord}{'support'} = $ctins;
                              print "\t\tINS: FIXED $draftBase\t$insertnt\t$head @ $oneBasedCoord" if($verbose);
                              #$failsafe = 0;
                              next POSITION;
                           }elsif($ctelement<=$MAXBASEDEL){### used to be else.  This insures the delete loop doesnt run 86x or howmany times there are inserts. Placing delete loop outside the indel block doesn't work as well.FIX FOR CASES WHERE INSERT<DEL
                              #### DELETIONS
                              $_ = $seq;
                              my $delbase = substr($seq, $pos+$k-1, $ctelement);
                              substr($_, $pos+$k-1, $ctelement) = '';###delete the last base

                              ###check the next k kmers
                              my $ctdel=0;
                              for(my $subset=$pos;$subset<=($pos+$k-2);$subset+=$i){   #### this code prevents changing the seq string when base change caused by indels
                                 my $subkmer = substr($_,$subset,$k);
                                 next POSITION if($subkmer=~/[EFIJLNOPQUXZ]/);  ###this prevents stalling when bases other than ATCG are encountered. This is more likely to occur with a deletion, even though k kmers were checked above for the presence of non-ATCG. Insertions would just pushed the Ns further 3'
                                 $ctdel++ if($bloom->contains($subkmer));
                              }
                              if($ctdel >= $delkct){###
                                 ###This changes the base in the original
                                 $seq = $_;
                                 $changeBase = "-" . $delbase;
                                 $track_all->{$head}{$oneBasedCoord}{'draft'} = $draftBase;
                                 $track_all->{$head}{$oneBasedCoord}{'change'} = $changeBase;
                                 $track_all->{$head}{$oneBasedCoord}{'support'} = $ctdel;
                                 print "\t\tDEL: FIXED $draftBase\t$changeBase\t$head @ $oneBasedCoord" if($verbose);
                                 #$failsafe = 0;
                                 next POSITION;
                              }else{
                                 print "\tDEL nothing to do for $ctelement\n" if($verbose);
                              }

                           }### end indel block
                        }###end insert for loop
                     }### END SUBSTITUTION
                  }###end contains if
                  print "ALTERNATE KMER MISSING --- possible SV or very low coverage.\n" if($verbose);
               }###Go through each alternate 3' end base

               print "ALTERNATE KMERs MISSING --- possible SV or very low coverage.\n" if($verbose);
               #$failsafe++; ### this will keep a count of number of bases/kmers in a row with no match in BF
               #if($failsafe >= 10000){### this is how many times missing kmers have been investigated in a row without resolution
               #   my $errormsg = "Many kmers have been flagged as needing editing on $head, but without resolution at $oneBasedCoord. Ensure you have sufficient read coverage >= 15X\n";
               #   print LOG "$errormsg\n";
               #}

            }### Though the kmer was missing, other kmers downstream that contain the 3' base were found (in sufficient amount in the BF)
         }### ^if kmer doesn't exist, other wise -- kmer in draft exists, move on
      }###for loop iterating through kmers

      my $seqminuslastk = substr($seq, 0, (length($seq) - (2*$k)));###The last kmer can never be resolved (in fact last <2*k can's be resolved)
      $last2k = substr($seq, (length($seq) - (2*$k)), (2*$k)); 
      $newmaster .= $seqminuslastk;
      $adjustedpos = length($newmaster);
      $mastercoord = ($mpos + $MAXSTRINGLEN) - ((2*$k)*$chunknum);

   }###master seq loop

   $newmaster .= $last2k;

   return $track_all,$newmaster;
}


## We hope this code is useful to you -- Please send comments & suggestions to rwarren at bcgsc.ca
