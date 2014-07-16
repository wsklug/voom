#!/usr/bin/perl

printf "Enter value for l_B: ";
$l_B=<STDIN>;
chomp $l_B;
if(!($l_B=~m/^[\de\+-\.]+$/)) {
    $l_B = .012;
    printf "Value for l_B was bad, using default value of .012\n";
}

printf "Enter value for L: ";
$L=<STDIN>;
chomp $L;
if(!($L=~m/^[\de\+-\.]+$/)) {
    $L = 2.0;
    printf "Value for L was bad, using default value of 2.0\n";
}

printf "Enter starting value for L/l_C: ";
$Ll_Cstart=<STDIN>;
chomp $Ll_Cstart;
if(!($Ll_Cstart=~m/^[\de\+-\.]+$/)) {
    $Ll_Cstart = 8.0;
    printf "Value for L/l_C was bad, using default start value of 8.0\n";
}

printf "Enter ending value for L/l_C: ";
$Ll_Cend=<STDIN>;
chomp $Ll_Cend;
if(!($Ll_Cend=~m/^[\de\+-\.]+$/)) {
    $Ll_Cend = 40.0;
    printf "Value for L/l_C was bad, using default end value of 40.0\n";
}

printf "Enter number of steps to take: ";
$nSteps=<STDIN>;
chomp $nSteps;
if(!($nSteps=~m/^[\de\+-\.]+$/)) {
    $nSteps = 4;
    printf "Value for # of steps was bad, using default value of 4\n";
}

printf "Enter crosslink spring constant: ";
$kcl=<STDIN>;
chomp $kcl;
if(!($kcl=~m/^[\de\+-\.]+$/)) {
    $kcl = -1.0;
    printf "Value for crosslink spring constant was bad, using default value of -1.0 (pinned nodes)\n";
}

printf "Enter system size as Wx/L,Wy/L or as a single number: ";
$ssize=<STDIN>;
chomp $ssize;
if(!($ssize=~m/^[\de\+-\.]+[,\de\+-\.]*$/)) {
    $Wx = 8.0;
    $Wy = 8.0;
    printf "Value for system size was bad, using default value of 8.0,8.0 \n";
}
else {
    @sizes = split(/,/,$ssize);
    $Wx = $sizes[0];
    $Wy = $sizes[0];
    if($#sizes > 0) {
	$Wy = $sizes[1];
    }
}

printf "Enter starting shear: ";
$shr_start=<STDIN>;
chomp $shr_start;
if(!($shr_start=~m/^[\de\+-\.]+$/)) {
    $shr_start = 0.0;
    printf "Value for starting shear was bad, using default value of 0.0\n";
}

printf "Enter ending shear: ";
$shr_end=<STDIN>;
chomp $shr_end;
if(!($shr_end=~m/^[\de\+-\.]+$/)) {
    $shr_end = 1.0;
    printf "Value for ending shear was bad, using default value of 1.0\n";
}

printf "Enter # of shear steps: ";
$nShrSteps=<STDIN>;
chomp $nShrSteps;
if(!($nShrSteps=~m/^[\de\+-\.]+$/)) {
    $nShrSteps = 20;
    printf "Value for # of shear steps was bad, using default value of 20\n";
}

printf "Enter name of parameter file: ";
$parFileName=<STDIN>;
chomp $parFileName;

open(initfile,"<$parFileName") or die "Must give a valid parameter filename!\n";
@paramdata=<initfile>;
close(initfile);
foreach $paramline (@paramdata) {
    if($paramline=~m/^kcl[^\s]*\s+/) {
	$paramline=~s/(\s)[\de\+-\.]+(\s)/$1$kcl$2/;
    }
    elsif($paramline=~m/^l_B[^\s]*\s+/) {
	$paramline=~s/(\s)[\de\+-\.]+(\s)/$1$l_B$2/;
    }
    elsif($paramline=~m/^Wx[^\s]*\s+/) {
	$paramline=~s/(\s)[\de\+-\.]+(\s)/$1$Wx$2/;
    }
    elsif($paramline=~m/^Wy[^\s]*\s+/) {
	$paramline=~s/(\s)[\de\+-\.]+(\s)/$1$Wy$2/;
    }
    elsif($paramline=~m/Shear_start/) {
	$paramline=~s/(\s)[\de\+-\.]+(\s)/$1$shr_start$2/;
    }
    elsif($paramline=~m/Shear_final/) {
	$paramline=~s/(\s)[\de\+-\.]+(\s)/$1$shr_end$2/;
    }
    elsif($paramline=~m/Shear_steps/) {
	$paramline=~s/(\s)[\de\+-\.]+(\s)/$1$nShrSteps$2/;
    }
    elsif($paramline=~m/^L\s+/) {
	$paramline=~s/(\s)[\de\+-\.]+(\s)/$1$L$2/;
    }
    
}

open(initfile,">$parFileName");
print initfile @paramdata;
close(initfile);

$Ll_Cstep = ($Ll_Cend-$Ll_Cstart)/$nSteps;
for($i=0; $i<=$nSteps; $i++) {
    $Ll_C = $Ll_Cstart + $i*$Ll_Cstep;
    open(initfile,"<$parFileName");
    @paramdata=<initfile>;
    close(initfile);
    foreach $paramline (@paramdata) {
	if($paramline=~m/l_C/) {
	    $paramline=~s/(\s)[\de\+-\.]+(\s)/$1$Ll_C$2/;
	}
    }
    open(initfile,">$parFileName");
    print initfile @paramdata;
    close(initfile);
    @progandargs = ("./shear", "$parFileName");
    system(@progandargs);
}

exit(0);
