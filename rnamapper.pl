#!/usr/bin/perl
use strict;
use warnings;

# Pre: valid input file
# Post: All contact strings for the aligned sequences
my (@seq_ref, @seq_str, @aligned_seqs);
my $input;
my $show_gaps = 0;

if (scalar(@ARGV) == 1) {
    $input = $ARGV[0];
}
elsif (scalar(@ARGV) == 2) {
    if ($ARGV[0] eq "-g") {
        $show_gaps = 1;
    }
    else {
        &usage();
    }
    $input = $ARGV[1];
}
else {
    &usage();
}

# Parsing data from files to variables
open(IN, "<$input");

my @line;
@line = split(' ', <IN>);
@seq_ref = split('', $line[1]);
@line = split(' ', <IN>);
@seq_str = split('', $line[1]);

while (<IN>) {
    my @line = split(' ', $_);
    my @i = split('v', $line[0]);
    push @{ $aligned_seqs[$i[1]-1] }, split('', $line[1]);
};
close(IN);

print "SQref\t".join("",@seq_ref)."\n";
print "SQstr\t".join("",@seq_str)."\n";

# Compute contact strings
for (my $i = 0; $i < scalar(@aligned_seqs); $i++) {
    my @secondary_structure = @{$aligned_seqs[$i]};
    my @contact_string = &get_contact_string(\@seq_ref, \@seq_str, \@{ $aligned_seqs[$i] });
    
    if (!$show_gaps) {
        &remove_gaps(\@secondary_structure, \@contact_string);
    }

    print "\n";
    print "SQv".($i+1)."\t".join("", @secondary_structure)."\n";
    print "\t".join("", @contact_string)."\n";
}
print "\n";

# Methods
sub get_contact_string() {
    my ($seq_ref, $seq_str, $aligned_seq) = @_;

    my @result = @{$seq_str};
    my %hairpin_positions = &extract_hairpin_positions($seq_str);
    
    foreach my $position (keys %hairpin_positions) {
        if (defined $hairpin_positions{$position}) {
            my $complementary_position = &get_complementary_hairpin_position(\%hairpin_positions, $position);
            if (@{$seq_ref}[$position] ne @{$aligned_seq}[$position] && (
                @{$aligned_seq}[$position] eq '-' ||
                !&is_complementary_nucleotide(
                    @{$aligned_seq}[$position], 
                    @{$aligned_seq}[$complementary_position]
                )
            )) {
                $result[$position] = '.';
                $result[$complementary_position] = '.';
                delete($hairpin_positions{$position});
                delete($hairpin_positions{$complementary_position});
            }
        }
    }
    return @result;
}

sub extract_hairpin_positions() {
    my ($contact_string) = @_;

    my %result;
    for (my $i = 0; $i < scalar(@{$contact_string}); $i++) {
        if ($$contact_string[$i] ne '.') {
            $result{$i} = $$contact_string[$i];
        }
    }
    return %result;
}

sub is_complementary_nucleotide() {
    my ($nucleotide_a, $nucleotide_b) = @_;
    if (($nucleotide_a eq 'a' && $nucleotide_b eq 'u') ||
        ($nucleotide_a eq 'u' && $nucleotide_b eq 'a') ||
        ($nucleotide_a eq 'c' && $nucleotide_b eq 'g') ||
        ($nucleotide_a eq 'g' && $nucleotide_b eq 'c')
    ) {
        return 1;
    }
    return 0;
}

sub get_complementary_hairpin_position() {
    my ($hairpin_positions, $position) = @_;
    my $result;

    my @pos_stack;
    if ($$hairpin_positions{$position} eq '(') {
        foreach my $key (sort {$a <=> $b} keys %$hairpin_positions) {
            my $value = $$hairpin_positions{$key};
            if ($value eq '(') {
                push @pos_stack, $key;
            }
            elsif (pop @pos_stack == $position) {
                $result = $key;
            }
        }
    }
    else {
        foreach my $key (reverse sort { $a <=> $b} keys %$hairpin_positions) {
            my $value = $$hairpin_positions{$key};
            if ($value eq ')') {
                push @pos_stack, $key;
            }
            elsif (pop @pos_stack == $position) {
                $result = $key;
            }
        }
    }
    return $result;
}

sub remove_gaps () {
    my ($sequence, $contact_string) = @_;
    for (my $i = 0; $i < scalar(@{$sequence}); $i++) {
        if ($$sequence[$i] eq "-") {
            splice(@{$sequence}, $i, 1);
            splice(@{$contact_string}, $i, 1);
            $i--;
        }
    }
}

sub usage() {
    print "Usage: perl rnamapper.pl [-g] FILE\n";
    print "Mapping RNA Secondary Structure into an Alignment .\n";
    print "\n";
    print "Mandatory arguments to long options are mandatory for short options too.\n";
    print "-g                  print the gaps of the secondary structures .\n";

    exit();
}
