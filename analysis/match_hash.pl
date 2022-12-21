%hash1;
open(F1,"$ARGV[0]") or die;
while(<F1>){
chomp $_;
$hash1{$_}=1;
}

%hash2;
open(F2,"$ARGV[1]") or die;
while(<F2>){
chomp $_;
($gene,$trf)=split(/==/,$_);
$hash2{$gene}=$trf;
}

foreach $ky (keys %hash1){
if (exists $hash2{$ky}){
#print "$ky\n";
print "$ky\t$hash2{$ky}\n";
}
}
