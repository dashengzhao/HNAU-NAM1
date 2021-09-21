#!/usr/bin/perl -w
#按照SNP间距差值划分，如0~100kb  100~200kb  200~300kb   求此bin窗口内所有R^2的均值
#purpose: to calculate the average of corresponding R^2 in each window.
#usage: perl $0  chr1_LDdecay.out chr1_LDdecay_average.out
#input file: chr1_LDdecay.out （extracted from the raw output of TASSEL_version5.2.74 software）
#output file: chr1_LDdecay_average.out
use List::Util qw/sum/;
my($window,$n,@line,$m,$z);
$window=100000;# 间隔100kb统计R^2均值
$n=10000000/$window; # 一般认为距离超过10Mb就不存在连锁关系
#==================================================================
open IN0,$ARGV[0];
open OUT,">./$ARGV[1]";
while(<IN0>){chomp;
	if($.==1){
		print OUT "Window_No\tID\tstart\tend\tDist_bp_aver\tR^2_average\n";
		next;
	}
	@line=split;
	next if ($line[13] eq "NaN");  ### R^2
	next if ($line[12]!~/\d+/);    ### Dist_bp
	$m=sprintf("%.6f",$line[12]/$window);
	$z=$1 if $m=~/(\d+)\./; #整数位
	if($m<1){
		my$tem_r='r_R0';
		my$tem_d='d_R0';
		push(@{$tem_r},$line[13]);
		push(@{$tem_d},$line[12]);
	}else{
		my$t_r='r_R'.$z;
		my$t_d='d_R'.$z;
		push(@{$t_r},$line[13]);
		push(@{$t_d},$line[12]);
	}
}
close IN0;
my($id_r,$id_d,$he_r,$he_d,$aver_r,$aver_d,$r1,$r2);
for(0..$n-1){
	$id_r="r_R".$_;
	$id_d="d_R".$_;
	$he_r=sum @{$id_r}; ##窗口内R方值之和
	$he_d=sum @{$id_d}; ##窗口内dis值之和
	$he_r=0 if (!$he_r);
	$he_d=0 if (!$he_d);
	my $num= @{$id_r};
	$num=1 if ($num ==0);
	$aver_r=sprintf("%.4f",$he_r/$num);	
	$aver_d=sprintf("%.4f",$he_d/$num);
	$r1=$_*$window;
	$r2=$r1+$window;
	print OUT "$id_r\t$_\t$r1\t$r2\t$aver_d\t$aver_r\n";
}

close OUT;		
