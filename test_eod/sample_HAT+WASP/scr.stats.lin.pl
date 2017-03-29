#!/usr/bin/perl

use POSIX;

$binsize = ($ARGV[0]) || 10;
$dispfield = ($ARGV[1]) || 1;

@value_list = ();
$count = 0.;
$sum = 0.;
$sum_of_square = 0.;
$sum_of_cube = 0.;
$sum_of_fourth = 0.;

while (<STDIN>) 
{
    next if (/^\#/); 
    $fieldvalue = (/(\S+)/g)[$dispfield-1];
    $output_hash{$binsize*(floor(($fieldvalue)/($binsize)))}++;

    $value_list[$count] = $fieldvalue;
    $sum += $fieldvalue;
    $sum_of_square += $fieldvalue * $fieldvalue;
    $sum_of_cube += $fieldvalue * $fieldvalue * $fieldvalue;
    $sum_of_fourth += $fieldvalue * $fieldvalue * $fieldvalue * $fieldvalue;
    ++$count;
}

if($count <= 0) { print STDOUT "# Count = $count\n"; exit; }

# Make histogrm and Calc Mode
print STDOUT "# Histogram \n";
$mode_key = "";
$mode_count = -99999999;
$cum = 0.;
$cumcount = 0;
$percentile05 = -1.;
$percentile50 = -1.;
$percentile95 = -1.;
foreach (sort {$a <=> $b} keys %output_hash) {
    $key = $_;
    $value = $output_hash{$key};
    if($value > $mode_count)
    {
	$mode_count = $value;
	$mode_key = $key; 
    }
    $cumcount +=$value;
    $cum += $value/$count;
    if($percentile05<0.)    {  if($cum>=0.159) { $percentile05 = $key; }    }
    if($percentile50<0.)    {  if($cum>=0.50) { $percentile50 = $key; }    }
    if($percentile95<0.)    {  if($cum>=0.841) { $percentile95 = $key; }    }
    print $key, "  ",  $key+$binsize, " ", $value/$count, " ",$cum," ", $value, " ", $cumcount, "\n";
}
$mode = $mode_key + ($binsize);
print STDOUT "# Mode: $mode  ($mode_count in $mode_key - ", $mode_key+$binsize , ") \n";


# Calc Median
@value_list_sorted = sort { $a <=> $b } @value_list;
if($count % 2 == 1 ) 
{
 $median_index = int(($count-1)/2); 
 $median = $value_list_sorted[$median_index];
}
else
{
 $median_index = int($count/2); 
 $median = ($value_list_sorted[$median_index] + $value_list_sorted[$median_index-1] ) /2.;
}

print STDOUT "# Median: $median \n";


# WARNING: I've never actually verified that this does what you expect
# The first few moments are probably fine, but I'd check before using 
# the higher momemts for anything important.

# Calc Moments
$average  = $sum / $count;
$average_of_square = $sum_of_square / $count;
$average_of_cube = $sum_of_cube / $count;
$average_of_fourth = $sum_of_fourth / $count;

$variance = $average_of_square - $average * $average;
$skewness = $average_of_cube - 3. * $average * $average_of_square + 2. * $average * $average * $average;
$kurtosis = $average_of_fourth - 4. * $average_of_cube * $average + 6. * $average_of_square * $average * $average - 3. * $average * $average * $average * $average;

$std_dev = sqrt($variance);
if($skewness<0.) { $skewness = -(-$skewness)**(1./3.); }
else             { $skewness =  ( $skewness)**(1./3.); }
if($kurtosis>=0.) { $kurtosis = sqrt(sqrt($kurtosis)); }
else { $kurtosis = -1.; }

printf(STDOUT "# Moments: %d  %.8f %.8f  %.8f %.8f \n",$count,$average,$std_dev,$skewness,$kurtosis );
printf(STDOUT "# Percentiles 5, 50, 95:  %.8f %.8f %.8f\n",$percentile50,$percentile95-$percentile50,$percentile50-$percentile05 );





