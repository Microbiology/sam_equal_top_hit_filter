# SAM File Equal Top Hit Filter

SAM file presents a single alignment for a given sequence as a "primary alignment", but this can be arbitrarily chosen during ties. The purpose of this code is to provide all primary alignments as well as other alignments from the same read which have identical quality to the primary alignment (defined as matching CIGAR string and mapping quality score).

## Test Input

The input file is a toy sam file with five references and two read pairs. This is captures in `test1.sam`.

- `r001`=`ref`   This is not a primary alignment but is identical to the primary alignment.
- `r001`=`ref`   This is another equal alignment to ref, although at a different position.
- `r001`=`ref2`   This is the primary alignment for r001
- `r001`=`ref3`   This is an alignment in which both pairs are lower quality than the primary alignemnt
- `r001`=`ref4`   This is an alignment in which one of the paired reads is lower quality than the primary alignment.
- `r001`=`ref5`   This is an alignment in which one of the paired reads is lower quality than the primary alignment.
- `r002`=`ref`   This is a primary alignment for the second read pair. This is not a "multi-mapper".

To run a test: `julia TopBestHits.jl -i test1.sam -o testout.txt`

## Expected Output

The output should include four sets of paired reads. `r001` should present the primary alignment against `ref2` as well as both matching secondary alignments to `ref`. `r002` should also be present. The `tsv` output should be the following, which includes the alignment sequence, reference genome, and the positions of the pairs as reference.

```
r002	ref	9	27
r001	ref	1	28
r001	ref	7	37
r001	ref2	5	31
```

## How It Works

This script runs in three steps. First a loop will identify "multimapper" reads by counting the frequency of the read pairs.

```julia
julia> counts
Dict{String,Int64} with 4 entries:
  "r001-Second" => 6
  "r002-Second" => 1
  "r001-First"  => 6
  "r002-First"  => 1
```

Then the script will create a dictionary with the primary alignments and their associated information.

```julia
julia> primaryd
Dict{String,Array} with 4 entries:
  "r001-Second" => Any["ref2", 5, 31, 30, "8M4I4M1D3M"]
  "r002-Second" => Any["ref", 9, 27, 40, "8M4I4M1D3M"]
  "r001-First"  => Any["ref2", 5, 31, 30, "8M4I4M1D3M"]
  "r002-First"  => Any["ref", 9, 27, 40, "8M4I4M1D3M"]
```

Finally the script puts it all together to create a dictionary with all primary alignments and the secondary alignments that match the primary alignments. A count is included so that we can pick only those with both matching pairs (denoted as a count of `2`).

```julia
julia> primarypairs
Dict{String,Int64} with 7 entries:
  "r001-ref3-6-32"  => 0
  "r002-ref-9-27"   => 2
  "r001-ref-1-28"   => 2
  "r001-ref4-15-36" => 1
  "r001-ref5-2-34"  => 1
  "r001-ref-7-37"   => 2
  "r001-ref2-5-31"  => 2
```

Then that's it. We just parse and print the results to the output.
