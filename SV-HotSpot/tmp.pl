open my $file, '<', "test_data/sv.bedpe"; 
my $firstLine = <$file>; 
print "$firstLine\n"; 

if ($firstLine =~ /^chrom1|^stuart1$/) {
  print "header is good\n";
} else {
  print "header is wrong\n";
}

close $file;

