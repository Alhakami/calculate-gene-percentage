#!/usr/bin/env perl

#  gene_coverage.pl
#
#
#  Created by Hind Alhakami on 23/1/16.
#


use warnings;
use lib "PATH_TO/perl5";
use Bio::Perl;
use Bio::SearchIO;


#input path to output.blast.txt files
foreach my $file_name (@ARGV) {
    my $in = new Bio::SearchIO(-format => 'blast', -file   => $file_name);

    my $total_genes_length = 0;
    my $total_coverage = 0;

    while( my $result = $in->next_result ) {

        ##  get query length
        my $query_length = $result->query_length;
        $total_genes_length +=  $query_length;

        if ( $result->num_hits == 0 ){next;}
        if ( $result->num_hits == 1 ){
            ## @coverage array will contain one element
            my $hit = $result->next_hit;
            while( my $hsp = $hit->next_hsp ) {
                if ( $hsp->rank == 1  &&  $hsp->percent_identity >= 75 ){
                    $total_coverage += $hsp->length('query');
                    last;
                }
            }
        }else{ ## in case of multiple hits, calculate accumulated coverage
            ## the overall coverage equals the sum of hits' coverage minus overlaps
            my @query_ranges;
            while( my $hit = $result->next_hit ) {
                while( my $hsp = $hit->next_hsp ) {
                    ## for each hit, get best alignment with 75% minimum identity
                    if ( $hsp->rank == 1 &&  $hsp->percent_identity >= 75 ){
                        my @range = $hsp->range('query');
                        push @query_ranges, [ @range[0,1], $hsp->length('query') ];
                        last;
                    }
                } ## while $hsp
            } ## while $hit

            if ( @query_ranges ){
            ## calculate hits accumulated coverage and overlaps
            ## sort hits by end positions
            my @sorted = sort {$a->[1] <=> $b->[1]} @query_ranges;
            ## get first sequence
            my $current = shift(@sorted);
            ## end position
            my $prev_end = (@$current)[1];
            ## initialize
            my $accumelated_coverage = 0;
            my $overlaps = 0;
            
            foreach my $row (@sorted){
                $accumelated_coverage +=(@$row)[2];
                my $current_start = (@$row)[0];
                if ( $current_start <= $prev_end ){
                    $overlaps += $prev_end -  $current_start + 1;
                }
                $prev_end = (@$row)[1];
            }
            
            ## calculated total coverage
            $accumelated_coverage -= $overlaps;
            $total_coverage += $accumelated_coverage;
            } 
        } ## else
    } ## while $results
    if ( $total_genes_length == 0) {$total_genes_length = 1;}
        my $coverage_ratio = $total_coverage / $total_genes_length;

    print $file_name, "\t",  $coverage_ratio, "\n";
    my $dirs = substr($file_name, 0, -length("output.blast.txt "));
    my $ofile = $dirs . "gene_coverage.txt";

    open(OUTPUT,">", $ofile) or die "Can't open data";
    print OUTPUT $coverage_ratio;
    close(OUTPUT);

} ##foreach assembly $file_name
