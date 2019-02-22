### Instructions for building the nthits Bloom filter utility PERL module
-------

1. BUILD a PERL5 module

Make sure you have swig installed and included in your path.

http://www.swig.org/


TO BUILD a Perl5 module (run in swig/):
```
a) /home/<user>/<path to swig>/preinst-swig -Wall -c++ -perl5 BloomFilter.i

b) g++ -c BloomFilter_wrap.cxx -I/usr/lib64/perl5/CORE -fPIC -Dbool=char -O3

The location of CORE may change on your system.
If you use linuxbrew, this works:
-I/<yourpath>/linuxbrew/Cellar/perl/5.28.0/lib/perl5/5.28.0/x86_64-linux-thread-multi/CORE/

c) g++ -Wall -shared BloomFilter_wrap.o -o BloomFilter.so -O3
```

TO COMPILE, swig needs the following Perl5 headers:
```C++
#include "Extern.h"
#include "perl.h"
#include "XSUB.h"
```
If they are not located in /usr/lib64/perl5, you can run "perl -e 'use Config; print $Config{archlib};" to locate them.


2. VERIFY your install

in the lib folder, execute:
$ ./test.pl
All tests should pass


3. CHANGE the relative path to BloomFilter.pm in ntEdit.pl/test.pl 

You only need to change if you have re-built in a relative directory different
from:
<pre>
use lib "$FindBin::Bin/./lib/"; (for LINKS and test.pl)
</pre>


