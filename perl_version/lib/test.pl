#!/gsc/btl/linuxbrew/bin/perl

use strict;
use POSIX;
use FindBin;
use lib "$FindBin::Bin/.";
use BloomFilter;
my $filter;

eval{
  $filter = new BloomFilter::BloomFilter("solidBF_k25.bf");
};

print "Reading existing Bloom filter....";

if($@){
   print "FAILED\n";
}else{
   print "PASSED.\n";
}


eval{
   $filter->insert("ACGTCAGATACGACTAATCAGCGAC");
};


print "Inserting kmer into existing Bloom filter....";

if($@){
   print "FAILED\n";
}else{
   print "PASSED.\n";
}

eval{
  if($filter->contains("ACGTCAGATACGACTAATCAGCGAC")) {
	print "Filter contained expected. \n";
  }else{
        print "Filter contained unexpected. \n";
  }
};


print "Checking kmer....";

if($@){
   print "FAILED\n";
}else{
   print "PASSED.\n";
}


print "Bloom filter tests complete.\n";

exit;
