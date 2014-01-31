#!/usr/bin/perl

use strict; 
my ( $n, $nl, $il, $key, $nb, $line, $lineo, $atm, $atm1, 
     $seq, $must_write, $k, $lineok, $lc_seq, $nlo,
     @atm_base, @t, @linek, @xyz, @xb, 
     %atm_cen, %atm_cen1, %has_CA, %has_C4, %nside
     );

# Extract three atoms from each residue

@atm_base = ( " C2 ", " C4 ",  " C6 " ); 

# Read in the data on the central atom of each residue

$n = 0;
while ( <DATA> ) { 
#  19591 SER  CB OG
  chop; @t = split ( / +/, $_ );
  if ( !$t[2] ) { next }
  $n++;  
  $atm_cen{$t[2]} = $t[3]; $nside{$t[2]} = 1;
  if ( $t[4] ) { $atm_cen1{$t[2]} = $t[4]; $nside{$t[2]} = 2 };
  }
# print STDERR "Have read in central atoms for $n residues.\n";

# print STDERR "atm_base= @atm_base\n";

# Read in all lines

$nl = 0; while ( <STDIN> ) { if ( /^ENDMDL/ ) { last }; $linek[$nl] = $_; $nl++ }

# Find and eliminate residues with no CA

for ( $il = 0; $il < $nl; $il++ ) { 

  $_ = $linek[$il]; chop;

#ATOM    426  CA  ASP    56A     26.191  53.931  14.169  1.00 49.50           C  
#            1234 123456789
#012345678901234567890123456789012345678901234567890123456789012345

  if ( substr($_,12,4) eq " CA "  ) { 
    $key = substr($_,17,9);
    $has_CA{$key}++;
  }
  if ( substr($_,12,4) eq " C4\*" ) { 
    $key = substr($_,17,9);
    $has_C4{$key}++;
  }
}

$nb = 0; $nlo = 0;
for ( $il = 0; $il < $nl; $il++ ) { 

  $_ = $linek[$il]; chop;

  if ( !/^ATOM  / ) { 
    if ( /^HEADER/ ) { print "$_\n"}
    if ( /^TER/    ) { print "$_\n"}
    next;
    }

  $key = substr($_,17,9);
#  print "key1= '$key', has_CA= '$has_CA{$key}' line=$_\n";
  if ( !$key ) { next }
  if ( !$has_CA{$key} && !$has_C4{$key} ) { next }

  $line = $_; $lineo = $line;

#ATOM   2281  C4*   G A 111      -7.624  73.628-124.140  1.00 32.96           C  
#                                     1       2       3
#                               fffffff fffffff fffffff
#012345678901234567890123456789012345678901234567890123456789012345
#          1         2         3         4         5         6
#            aaaa^    ^    ^

  $xyz[0] = substr($_,30,8); $xyz[1] = substr($_,38,8); $xyz[2] = substr($_,46,8);

  $atm = substr($line, 12, 4 ); $atm1 = $atm; $atm1 =~ s/ //g;
  $seq = substr($line, 17, 3 ); 

  $must_write = 0;
  if ( $has_CA{$key} ) { 
    if ( $atm eq " CA "  || $atm eq " O  " ) { $must_write = 1; if ( $atm eq " CA " ) { $nb = 0 } }
    elsif ( $nside{$seq} == 1 && $atm1 eq $atm_cen{$seq} ) {
      substr($lineo, 12, 4 ) = " CMA";
      $must_write = 1; $nb = 0;
    }
    elsif (  $nside{$seq} == 2 && ( $atm1 eq $atm_cen{$seq} || $atm1 eq $atm_cen1{$seq} ) ) {  
      if ( $nb == 0 ) { 
        $nb++; 
        for ( $k = 0; $k < 3; $k++ ) { $xb[$k] = $xyz[$k] }
        print STDERR "A: nb=$nb, atm1='$atm1'\n";
      }
      elsif ( $nb > 0 ) { 
        $nb++; 
        for ( $k = 0; $k < 3; $k++ ) { $xb[$k] += $xyz[$k] }
        print STDERR "B: nb=$nb, atm1='$atm1'\n";
        if ( $nb == 2 ) { 
          for ( $k = 0; $k < 3; $k++ ) { $xb[$k] /= $nb }
          $nb = 0; 
          $must_write = 2;
          substr($lineo, 12, 4 ) = " CMA";
        }
        print STDER "C: nb=$nb, atm1='$atm1',must_write=$must_write\n";
      }
    } 

    if ( $seq eq "PCA" ) { $seq = "GLU" }; if ( $seq eq "CYH" ) { $seq = "CYS" }
    substr($lineo,17,3) = $seq;
  }

  if ( $has_C4{$key} ) {
    $lc_seq = $seq; $lc_seq = lc ( $lc_seq );
    substr($lineo, 17, 3 ) = $lc_seq;
    if ( $atm eq " C4\*" || $atm eq " P  " ) { $must_write = 1 } 
    elsif ( $atm eq $atm_base[0] || $atm eq $atm_base[1] || $atm eq $atm_base[2] ) {
      if ( $nb == 0 ) { 
        $nb++; 
        for ( $k = 0; $k < 3; $k++ ) { $xb[$k] = $xyz[$k] }
      }
      elsif ( $nb > 0 ) { 
        $nb++; 
        for ( $k = 0; $k < 3; $k++ ) { $xb[$k] += $xyz[$k] }
        if ( $nb == 3 ) { 
          for ( $k = 0; $k < 3; $k++ ) { $xb[$k] /= $nb }
          $nb = 0; 
          $must_write = 2;
          substr($lineo, 12, 4 ) = " CMN";
        }
      }
    } 
  }

  if ( $must_write == 1 ) { 
    print "$lineo\n"; $nlo++;
    if ( $seq eq "GLY" && $atm eq " CA " ) { $lineok = $lineo; substr($lineok, 12, 4 ) = " CMA" }
    if ( $seq eq "GLY" && $atm eq " O  " ) { print "$lineok\n"; $nlo++ }
  }
  elsif ( $must_write == 2 ) { $nlo++; printf ( "%28s  %8.3f%8.3f%8.3f%s\n", substr($lineo,0,28), $xb[0], $xb[1], $xb[2], substr($lineo,54) ) }

}

print STDERR "Have written $nlo lines to STDOUT.\n";


__END__
  50306 LEU  CG  
  45937 ALA  CB  
  43003 GLY  CA  
  39788 VAL  CB  
  34559 GLU  CD  
  32774 ASP  CG  
  31911 THR  CB  
  31000 LYS  CD  
  27992 ILE  CG1 
  26906 ARG  NE  
  25008 ASN  CG  
  24878 PRO  CG  
  21935 PHE  CG  
  20666 GLN  CD 
  19591 SER  CB OG
  12890 HIS  CG  
  11810 MET  CG  
   9043 TYR  CD1 CD2
   8267 TRP  CD2 
   4983 CYS  CB SG

