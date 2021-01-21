#!/usr/bin/env julia
# Geof Hannigan
# Purpose of this script is to find all SAM file seconary alignments
# whose quality match their associated primary alignment.

# using Pkg
# Pkg.add("Bio")
# Pkg.add("ArgParse")

# Use the package BioAlignments
# https://biojulia.net/BioAlignments.jl/latest/hts-files.html#SAM-and-BAM-Records-1
using BioAlignments
using DelimitedFiles
using ArgParse

# Setup the command line arguments
s = ArgParseSettings(
	version = "0.1",
	add_version = true
	)
@add_arg_table s begin
	"-i", "--input"
		help = "Input file name."
		required = true
	"-o", "--output"
		help = "Output file name."
		required = true
end

parsed_args = parse_args(ARGS, s)

println("Parsed args:")
for (arg,val) in parsed_args
    println("  $arg  =>  $val")
end


# Open the SAM file
reader = open(SAM.Reader, parsed_args["input"])

# Testing can look like this:
# reader = open(SAM.Reader, "test1.sam")

# Make a count hash to capture multi-mappers
# Here we do it only at the read level
counts = Dict{String,Int64}()
for record in reader
	if SAM.ismapped(record)
		# Get the name
		tname = SAM.tempname(record)
		tflag = SAM.flag(record)
		flagresult = collect(last(bitstring(tflag), 12))
		if Int64(flagresult[6]) == 49
			rpair = "First"
			# Make sure to include position for better ID
			tpos = SAM.position(record)
			trpos = SAM.nextposition(record)
		elseif Int64(flagresult[5]) == 49
			rpair = "Second"
			# Reverse for the other paired read
			trpos = SAM.position(record)
			tpos = SAM.nextposition(record)
		else
			# println("Looks like its not paired")
			rpair = "Fail"
		end

		# For the multimappers, we only want to look at the sequence read name
		# and whether it if first or second. Location doesn't matter.

		tcombo = string(tname, "-", rpair)

		# First or second in pair
		# println(tcombo)
		if tcombo in keys(counts)
			counts[tcombo] += 1
		else
			counts[tcombo] = 1
		end
	end
end

println("Printing Counts:")
for (key,val) in counts
    println("  $key  =>  $val")
end

# Get a dictionary of the primary alignments
# Should only be one primary alignment per read
# Key is the read, value is the reference

reader = open(SAM.Reader, parsed_args["input"])
# reader = open(SAM.Reader, "test1.sam")

primaryd = Dict{String,Array}()
for record in reader
	if SAM.ismapped(record)
		# This needs to be defined so that the scope carries
		# through to the rest of the loop
		mapq = 255

		tname = SAM.tempname(record)
		# println(tname)
		refn = SAM.refname(record)
		# println(refn)
		# Exception handling in case mapping quality is 255,
		# which means it does not exist
		try
			mapq = Int64(SAM.mappingquality(record))
		catch
			# 255 means it does not exist
			mapq = 255
		end
		# println(mapq)
		cig = SAM.cigar(record)

		tflag = SAM.flag(record)
		flagresult = collect(last(bitstring(tflag), 12))
		if Int64(flagresult[6]) == 49
			rpair = "First"
			# Make sure to include position for better ID
			tpos = SAM.position(record)
			trpos = SAM.nextposition(record)
		elseif Int64(flagresult[5]) == 49
			rpair = "Second"
			# Reverse for the other paired read
			trpos = SAM.position(record)
			tpos = SAM.nextposition(record)
		else
			# println("Looks like its not paired")
			rpair = "Fail"
		end

		tcombo = string(tname, "-", rpair)

		# println(tcombo)
		if tcombo in keys(primaryd)
			# println("Warning: Potentially multiple primary alignments!")
		else
			if SAM.isprimary(record)
				# println("Adding")
				# println(mapq)
				primaryd[tcombo] = [refn, tpos, trpos, mapq, cig]
			else
				# println("Not a primary alignment.")
			end
		end
	end
end

println("Printing Primary Aligment Stats:")
for (key,val) in primaryd
    println("  $key  =>  $val")
end

# Find matching best hits among each read
# Report the read and alignment pair
# reader = open(SAM.Reader, "test1.sam")
reader = open(SAM.Reader, parsed_args["input"])

# Count primary pairs
primarypairs = Dict{String,Int64}()

for record in reader
	# Only look at what mapped, always
	if SAM.ismapped(record)
		mapq = 255

		#Get the ID information
		tname = SAM.tempname(record)
		refn = SAM.refname(record)
		try
			mapq = Int64(SAM.mappingquality(record))
		catch
			# 255 means it does not exist
			mapq = 255
		end
		cig = SAM.cigar(record)

		tflag = SAM.flag(record)
		flagresult = collect(last(bitstring(tflag), 12))
		if Int64(flagresult[6]) == 49
			rpair = "First"
			# Make sure to include position for better ID
			tpos = SAM.position(record)
			trpos = SAM.nextposition(record)
		elseif Int64(flagresult[5]) == 49
			rpair = "Second"
			# Reverse for the other paired read
			trpos = SAM.position(record)
			tpos = SAM.nextposition(record)
		else
			# println("Looks like its not paired")
			rpair = "Fail"
		end

		tcombo = string(tname, "-", rpair, "-", tpos, "-", trpos)

		# We want to count the read-pair combos because we only want
		# to keep those with two (a full paired match)

		combotr = string(tname, "-", refn, "-", tpos, "-", trpos)
		multimapperid = string(tname, "-", rpair)

		# And establish the counts starting at zero
		if combotr in keys(primarypairs)
			#nothing
		else
			primarypairs[combotr] = 0
		end

		# println("Searching ", tcombo, "-", refn)
		# println("Multimapper ID is ", multimapperid)
		# Check if the read is a multimapper
		if counts[multimapperid] > 1
			# Check if it is a secondary hit
			if SAM.isprimary(record)
				# println(tcombo, "---", refn)
				# Keep track of primary alignment counts
				primarypairs[combotr] += 1
			else
				# Determine if this hit is the same as the
				# primary alignment of the read
				aref = primaryd[multimapperid][1]
				refstart = primaryd[multimapperid][2]
				refend = primaryd[multimapperid][3]
				aqual = primaryd[multimapperid][4]
				acigar = primaryd[multimapperid][5]
				if isequal(tpos, refstart) && isequal(trpos, refend)
						# nothing
				else
					if isequal(mapq, aqual)
						if isequal(cig, acigar)
							# println(tcombo, "---", refn)
							primarypairs[combotr] += 1
						end
					end
				end
			end
		elseif counts[multimapperid] == 1
			# println(tcombo, "---", refn)
			primarypairs[combotr] += 1
		end
	end
end

println("Printing Primary Pairs:")
for (key,val) in primarypairs
    println("  $key  =>  $val")
end

# Print output to a 2 dim array for now
outarray = Array{Char}(undef, 0, 4)
for (key, value) in primarypairs
	if isequal(value, 2)
		karray = split(key, "-")
		tread = karray[1]
		tref = karray[2]
		spos = karray[3]
		srpos = karray[4]
		# println(karray)
		global outarray = [outarray; [tread tref spos srpos]]
		# push!(outarray, karray)
	end
end

# Write the final array to tsv file for reference
# writedlm( "testout.txt",  outarray, '\t')
writedlm( parsed_args["output"],  outarray, '\t')


