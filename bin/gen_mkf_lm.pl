#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

Generate SFBA Makefile

=head1 SYNOPSIS

./gen_mkf.pl [options] 

Options:
  -h      help messages
  -d      degub
  --mf    output make file

=head1 OPTIONS

=over 8

=item B<-help>

  -w         work directory : location for all output files
  -t         directory for the C++ executable file for Estep (running MCMC)
  --geno     specify genotype format: (gzipped) vcf, (gzipped) genotxt, or bed
  --gd       genotype file directory
  --pheno    phenotype file
  -f         file with a list of fileheads for genotype files
  -G         genotype format: GT(genotype data 0/0, 0/1, 1/1) or EC (dosage data)
  --maf      maf threshold: default 0.5% 
  --lm       specify lm mode: 1, 2, 3, 4
  --mem      specify maximum memory usage: default 3000MB
  --time     specify time of runing MCMC per block per MCMC iteration: default 24hr
  -l         specify how job will be submited (options: slurm, mosix, local): default slurm
  --mf       output make file
  --nice     SLURM option for scheduling priority
  --xnode    compuation nodes to be excluded
  -j         specify SLURM job names
  --wnode    specify compuation nodes to be used
  --part     specify SLURM partition

=back

=head1 DESCRIPTION

B<gen_mkf.pl> will generate a makefile for conducting Bayesian GWASs by SFBA. 

=cut

# define default option variables
my $help;
my $verbose;
my $debug;
my $man;
my $launchMethod = "slurm";
my $wkDir=getcwd();
my $makeFile = "lm_assoc.mk";
my $genofile = "vcf";

my $toolE="/bin/Estep_mcmc";
my $genoDir = "";
my $pheno="";
my $filelist = "/net/fantasia/home/yjingj/GIT/SFBA_example/ExData/fileheads_4region.txt";

my $GTfield="GT";
my $maf="0.005";

my $maxmem = "3000";
my $time = "24:00:00";
my $nice = "0";
my $jobid="";
my $xnode="";
my $wnode="";
my $part="nomosix";

my $lm="3";


#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help, 'v'=>\$verbose, 'd'=>\$debug, 'm'=>\$man,
                'w:s'=>\$wkDir, 'Estep:s' =>\$toolE, 
                'geno:s'=>\$genofile,
                'gd:s'=>\$genoDir, 'pheno:s'=>\$pheno,
                'G:s'=>\$GTfield, 'maf:s'=>\$maf,
                'mem:s'=>\$maxmem, 'time:s'=>\$time, 'f:s'=>\$filelist, 
                'l:s'=>\$launchMethod, 'mf:s'=>\$makeFile, 'nice:s'=>\$nice, 
                'j:s' =>\$jobid, 'xnode:s'=>\$xnode, 'wnode:s'=>\$wnode, 'part:s'=>\$part)
  || !defined($wkDir) || scalar(@ARGV)!=0)
{
    if ($help)
    {
        pod2usage(1);
        exit(0);
    }
    elsif($man) {
        pod2usage(-exitval=>0, -verbose => 2);
    }
    else
    {
        pod2usage(1);
        exit(0);
    }
}

if ($help)
    {
        pod2usage(1);
        exit(0);
    }
elsif($man) {
        pod2usage(-exitval=>0, -verbose => 2);
    }


if ($launchMethod ne "local" && $launchMethod ne "slurm" && $launchMethod ne "mosix" && $launchMethod ne "qsub")
{
    print STDERR "Launch method has to be local, slurm, mosix, or qsub\! \n";
    exit(1);
}

##############
#print options
##############
printf("Options\n");
printf("\n");
printf("launch method : %s\n", $launchMethod);
printf("work directory : %s\n", $wkDir);
print "Estep: ", $toolE, "\n";
print "genoDir: ", $genoDir, "\n", 
        "pheno: ", $pheno,  
        "GTfield ", $GTfield, "; maf ", $maf, "\n";
printf("\n");


#arrays for storing targets, dependencies and commands
my @tgts = ();
my @deps = ();
my @cmds = ();

#temporary variables
my $tgt;
my $dep;
my @cmd;

mkpath($wkDir);

########################################
# Initial Set up before EM iterations
########################################

$tgt = "$wkDir/pre_em.OK";
$dep = "";
@cmd = "rm -f -r $wkDir/output $wkDir/OUT";
push(@cmd, "mkdir -p $wkDir/output $wkDir/OUT");
makeJob("local", $tgt, $dep, $wkDir, @cmd);  


###### EM step 0 without dependencies ###########
my $premcmcOK="";
my @filehead;
my $line;

open(my $FILELIST, $filelist)
    or die "Can not open $filelist \!";
 while ($line = <$FILELIST>) {
    chop $line;
    push(@filehead, $line);
}
close $FILELIST;

if(@filehead == 0) {
    print STDERR "file list is empty\! \n";
    exit(1);
} 
else{ print "Total \# of fileheads: ", scalar(@filehead), "\n \n"; }


for(my $j=0; $j< @filehead; ++$j)
    {
        $line = $filehead[$j];
        $premcmcOK .= "$wkDir/OUT/$line.OK ";
        $tgt = "$wkDir/OUT/$line.OK";
        $dep = "$wkDir/pre_em.OK";
        if($genofile eq "vcf"){
           @cmd = "$toolE -vcf $genoDir/$line.vcf.gz -p $pheno -GTfield $GTfield -maf $maf -lm $lm -o $line > $wkDir/OUT/$line.output.txt";
        }elsif ($genofile eq "genotxt"){
          @cmd = "$toolE -g $genoDir/$line.geno.gz -p $pheno -GTfield $GTfield -maf $maf -lm $lm -o $line  > $wkDir/OUT/$line.output.txt";
        }elsif ($genofile eq "bed") {
          @cmd = "$toolE -bfile $genoDir/$line -GTfield $GTfield -maf $maf -lm $lm -o $line > $wkDir/OUT/$line.output.txt";
        }else{
          die "genoDir need to be one of the vcf, genotxt, or bed file types\!\n"
        }

        makeJob($launchMethod, $tgt, $dep, $wkDir, @cmd);
    }

my $assocfile="$wkDir/assoc_all.txt";
my $logfile="$wkDir/log_all.txt";

$tgt = "$wkDir/cp_assoc.OK";
$dep = "$premcmcOK";
@cmd = "cat \`ls -d -1 $wkDir/output/** | grep assoc | sort\` | sort -k1,1n -k3,3n > $assocfile";
push(@cmd, "cat \`ls -d -1 $wkDir/output/** | grep log | sort\` > $logfile");
makeJob("local", $tgt, $dep, $wkDir, @cmd);


#*******************
#Write out make file
#*******************
open(MAK,">$makeFile") || die "Cannot open $makeFile\n";
#this tells GNU Make to delete files that were generated during a command that did not end successfully.
print MAK ".DELETE_ON_ERROR:\n\n";
#this ensures that all the targets are executed; exclude clean
print MAK "all: @tgts\n\n";

######## Create clean jobs command #######
mkpath("$wkDir/slurm_err");
push(@tgts, "clean_err");
push(@deps, "");
push(@cmds, "\t-rm -rf $wkDir/slurm_err/*.err");


push(@tgts, "clean");
push(@deps, "");
push(@cmds, "\t-rm -rf $wkDir/*.OK $wkDir/OUT/*.OK $wkDir/slurm_err/*.err");


for(my $i=0; $i < @tgts; ++$i)
{
    print MAK "$tgts[$i]: $deps[$i]\n";
    print MAK "$cmds[$i]\n";
}
close MAK;

##########
#functions
##########

#run a job either locally or by slurm
sub makeJob
{
    my ($method, $tgt, $dep, $wkd, @cmd) = @_;

    if ($method eq "local")
    {
        makeLocalStep($tgt, $dep, @cmd);
    }
    elsif ($method eq "slurm")
    {
        makeSlurm($tgt, $dep, $wkd, @cmd);
    }
    elsif($method eq "mosix")
    {
      makeMosix($tgt, $dep, $wkd, @cmd);
    }
}

#run mosix jobs
sub makeMosix
{
    my ($tgt, $dep, $wkd, @cmd) = @_;
    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmdtemp = "";
    for my $c (@cmd)
    {
      if($wnode){
        $cmdtemp .= "\tmosbatch -E$wkd -m$maxmem -j$wnode " . $c . "\n";
      }else{
        $cmdtemp .= "\tmosbatch -E$wkd -m$maxmem -b " . $c . "\n";
      }

    }
    $cmdtemp .= "\ttouch $tgt\n";
    push(@cmds, $cmdtemp);
}

#run slurm jobs
sub makeSlurm
{
    my ($tgt, $dep, $wkd, @cmd) = @_;
    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmdtemp = "";
    for my $c (@cmd)
    {
        $cmdtemp .= "\tsrun --exclude=$xnode --partition=$part --mem-per-cpu\=$maxmem --time\=$time --nice\=$nice --error\=$wkDir/slurm_err/\%N.\%j.err -J $jobid -D $wkd $c \n";
    }
    $cmdtemp .= "\ttouch $tgt\n";
    push(@cmds, $cmdtemp);
}

#run a local job
sub makeLocalStep
{
    my ($tgt, $dep, @cmd) = @_;

    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = "";
    for my $c (@cmd)
    {
        $cmd .= "\t" . $c . "\n";
    }
    $cmd .= "\ttouch $tgt\n";
    push(@cmds, $cmd);
}

# Check empty directories
sub dir_is_empty
{
    my ($path) = @_;
    opendir DIR, $path;
    while(my $entry = readdir DIR) {
        next if($entry =~ /^\.\.?$/);
        closedir DIR;
        return 0;
    }
    closedir DIR;
    return 1;
}
