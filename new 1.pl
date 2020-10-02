$a = "./builder_sample.exe ";
for ($i = 1.0; $i < 10.0; $i=$i+1.0){
	$j=1.1 + $i/2.0;
	my @b = $a;
	push(@b,"res".$i.".txt");
	print($j);
	push(@b,$j);
	print(@b);
	system(@b);
}
