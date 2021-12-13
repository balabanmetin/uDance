#!/usr/bin/perl -w
# Compares bootstrap and length for splits that match
use strict;
use lib "$ENV{HOME}/Genomics/Perl/modules";
use MOTree;
use Getopt::Long;

# Known bugs:
# If comparing rooted trees, the length of a split at the root should be the
# sum of the two lengths at the root, but it only uses the 1 branch length

{
    my $usage = "CompareTree.pl [-noboot] [ -leaf ] [ -simplify ] [ -wanders ] [-debug]\n"
	. "                     -tree treeFile1 -versus treeFile2\n"
	. "   Compares bootstrap and length for splits in tree1 that are found\n"
	. "   in tree2. If a split in tree1 is absent from tree2, has no\n"
	. "   values for tree2\n"
	. "   Warning: this is designed for unrooted trees in which every internal node\n"
	. "   has more than one child. For rooted trees (i.e. the root has two children),\n"
	. "   the splits around the root's children may be reported more than once.\n"
	. "   This could also lead to small errors in the reported Frac of splits that match\n"
	. "   If -simplify is set, then leaves that are absent from either tree are removed\n"
	. "   before doing the comparison\n"
	. "   If -leaf is set, includes a row for each leaf (to compare branch lengths)\n";
    
    my $file1 = undef;
    my $file2 = undef;
    my $debug = 0;
    my $simplify = 0;
    my $wanders = 0;
    my $leaf = 0;
    my $noboot = 0;
    if (@ARGV == 2) {
	($file1,$file2) = @ARGV;
    } else {
	die $usage unless GetOptions('tree=s' => \$file1,
				     'versus=s' => \$file2,
				     'debug' => \$debug,
				     'simplify' => \$simplify,
				     'wanders' => \$wanders,
				     'leaf' => \$leaf,
				     'noboot' => \$noboot)
	    && defined $file1 
	    && defined $file2
	    && @ARGV == 0;
    }
    open(FILE1,"<",$file1) || die "Cannot read $file1";
    open(FILE2,"<",$file2) || die "Cannot read $file2";

    my $nTree = 0;
    while(1) {
	my $tree1 = MOTree::new(fh => \*FILE1);
	if (!defined $tree1) {
	    die "Cannot read tree from $file1" if $nTree == 0;
	    last;
	    next;
	}
	my $tree2 = MOTree::new(fh => \*FILE2);
	$nTree++;
	die "Cannot read tree $nTree from $file2" unless defined $tree2;

	# remove leading and trailing underscores, created as artifacts of some tree conversions e.g. Rose
	foreach my $node (@{ $tree1->depthfirst() }) {
	    if ($tree1->is_Leaf($node)) {
		my $id = $tree1->id($node);
		$id =~ s/^_+//;
		$id =~ s/_+$//;
		$tree1->set_id($node,$id);
	    }
	}

	if ($simplify) {
	    my %id2 = map {$tree2->id($_) => $_} $tree2->get_leaf_nodes();
	    my @remove1 = grep {!exists $id2{$tree1->id($_)}} $tree1->get_leaf_nodes;
	    $tree1->removeNodes(\@remove1);
	    my %id1 = map {$tree1->id($_) => $_} $tree1->get_leaf_nodes();
	    my @remove2 = grep {!exists $id1{$tree2->id($_)}} $tree2->get_leaf_nodes;
	    $tree2->removeNodes(\@remove2);
	}

	foreach my $node (@{ $tree2->depthfirst() }) {
	    if ($tree2->is_Leaf($node)) {
		my $id = $tree2->id($node);
		$id =~ s/^_+//;
		$id =~ s/_+$//;
		$tree2->set_id($node,$id);
	    }
	}

	my $match = $tree1->Compare($tree2);
	if ($nTree==1) {
	    print join("\t", qw{Len1 Boot1 Len2 Boot2 NDesc LeafA LeafB Node1 Node2})."\n";
	}
	my $maxlenDiff = undef;
	my $maxbootDiff = undef;
	my $nDesc1 = $tree1->countAllLeaves(); # node => number of leaves beneath it
	my $totLeaf1 = $nDesc1->{$tree1->get_root_node()};
	my $nFound = 0;
	my $nInternal = 0;
	my %wander = map { $_ => 0 } $tree1->get_leaf_nodes();
	my $sumlen1 = 0;
	my $sumlen2 = 0;
	foreach my $node1 (@{ $tree1->depthfirst() }) {
	    next if $nDesc1->{$node1} <= 1 || $nDesc1->{$node1} >= $totLeaf1 - 1;
	    $nInternal++;

	    my $len1 = $tree1->branch_length($node1);
	    my $boot1 = $tree1->id($node1);
	    my $len2 = "";
	    my $boot2 = "";
	    my $node2 = undef;

	    my $match1 = $match->{$node1};
	    if (defined $match1) {
		$nFound++;
		$node2 = $match1->[0];
		$len2 = $tree2->branch_length($node2);
		$boot2 = $tree2->id($node2);

		if (defined $len1 && defined $len2 && $len1 ne "" && $len2 ne "") {
		    my $diff = abs($len1-$len2);
		    $maxlenDiff = $diff if !defined $maxlenDiff || $diff > $maxlenDiff;
		}
		if (!$noboot && (defined $boot1 && defined $boot2 && $boot1 ne "" && $boot2 ne "")) {
		    my $diff = abs($boot1-$boot2);
		    $maxbootDiff = $diff if !defined $maxbootDiff || $diff > $maxbootDiff;
		}
		if (defined $len1 && defined $len2) {
		    $sumlen1 += $len1;
		    $sumlen2 += $len2;
		}

		$len2 = "NA" if !defined $len2;
		$boot2 = "NA" if !defined $boot2;
	    } elsif ($wanders) {
		# Increment wander count on both sides of this not-found split
		# Note this increases time to O(N**2) but could probably be reduced to O(N) by using nDesc
		# and summing down from the root to get the value for each leaf
		# (but would need careful bookkeeping)
		my $beneath = $tree1->all_descendents($node1);
		my %beneath = map {$_=>1} grep{exists $wander{$_}} @$beneath;
		my $nLeafBeneath = scalar(keys %beneath);
		foreach my $leaf (keys %wander) {
		    if (exists $beneath{$leaf}) {
			$wander{$leaf} += 1/$nLeafBeneath;
		    } else {
			$wander{$leaf} += 1.0/($totLeaf1 - $nLeafBeneath);
		    }
		}
	    }
	    my @children1 = $tree1->children($node1);
	    my $leafA = $children1[0];
	    while (! $tree1->is_Leaf($leafA)) { my @l = $tree1->children($leafA); $leafA = $l[0]; };
	    my $leafB = $children1[-1];
	    while (! $tree1->is_Leaf($leafB)) { my @l = $tree1->children($leafB); $leafB = $l[-1]; };
	    my @fields = map { defined $_ ? $_ : "" } ($len1, $boot1,
						       $len2, $boot2,
						       $nDesc1->{$node1},
						       $tree1->id($leafA), $tree1->id($leafB),
						       $node1, $node2);

	    print join("\t", @fields)."\n";
	}
	if ($leaf) {
	    my %map2 = map { $tree2->id($_) => $_ } $tree2->get_leaf_nodes();
	    foreach my $leaf1 ($tree1->get_leaf_nodes()) {
		my $id = $tree1->id($leaf1);
		my $leaf2 = $map2{$id};
		die "Cannot match leaf $leaf1 id $id" unless defined $leaf2;
		my $len1 = $tree1->branch_length($leaf1);
		my $len2 = $tree2->branch_length($leaf2);
		print join("\t", defined $len1 ? $len1 : "",
			   "",
			   defined $len2 ? $len2 : "",
			   "",
			   1, $id, $id, $leaf1, $leaf2)."\n";
		if (defined $len1 && defined $len2) {
		    $sumlen1 += $len1;
		    $sumlen2 += $len2;
		}
	    }
	}
	print STDERR join("\t","Splits","Found",$nFound,"Total",$nInternal,
			  "Frac",$nInternal > 0 ? sprintf("%.4g", $nFound/$nInternal) : "",
			  "MaxLnDf",defined $maxlenDiff ? sprintf("%.3g",$maxlenDiff) : "",
			  "Ratio", $sumlen1 > 0 && $sumlen2 > 0 ? sprintf("%.3g",$sumlen1/$sumlen2) : "",
			  "MaxBtDf",defined $maxbootDiff? sprintf("%.3g",$maxbootDiff) : ""
			  )."\n";
	if ($wanders) {
	    # also print the top few wanderers
	    my @leaves = sort {$wander{$b} <=> $wander{$a}} (keys %wander);
	    if ($wanders==1) {
		printf STDERR "Top Wanderers: ";
		for (my $i = 0; $i < 10 && @leaves > 0; $i++) {
		    print STDERR " " . $tree1->id(shift @leaves);
		}
		print STDERR "\n";
	    } else {
		foreach my $leaf (@leaves) {
		    print STDERR join("\t", $tree1->id($leaf), $wander{$leaf})."\n";
		}
	    }
	}
    }
}

