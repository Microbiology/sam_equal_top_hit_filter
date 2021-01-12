#!/usr/bin/env julia
# Geof Hannigan
# Purpose of this script is to find all SAM file seconary alignments
# whose quality match their associated primary alignment.

# TODO: I should also add position in case the same read maps multiple
# portions of the same genome

using Pkg
Pkg.add("Bio")

# Use the package BioAlignments
# https://biojulia.net/BioAlignments.jl/latest/hts-files.html#SAM-and-BAM-Records-1
using BioAlignments
using DelimitedFiles

# Open the SAM file
reader = open(SAM.Reader, "test.sam")
# Make a count hash to capture multi-mappers
counts = Dict{String,Int64}()
for record in reader
	if SAM.ismapped(record)
		# Get the name
		tname = SAM.tempname(record)
		tflag = SAM.flag(record)
		flagresult = collect(last(bitstring(tflag), 12))
		if Int64(flagresult[6]) == 49
			rpair = "First"
		elseif Int64(flagresult[5]) == 49
			rpair = "Second"
		else
			println("Looks like its not paired")
			rpair = "Fail"
		end
		# We should also be including the position in case
		# the read is mapped multiple times on same reference
		tpos = SAM.position(record)
		tcombo = string(tname, "-", rpair, "-", tpos)
		# First or second in pair
		println(tcombo)
		if tcombo in keys(counts)
			counts[tcombo] += 1
		else
			counts[tcombo] = 1
		end
	end
end

# Get a dictionary of the primary alignments
# Key is the read, value is the reference
reader = open(SAM.Reader, "test.sam")
primaryd = Dict{String,Array}()
for record in reader
	if SAM.ismapped(record)
		tname = SAM.tempname(record)
		refn = SAM.refname(record)
		mapq = Int64(SAM.mappingquality(record))
		cig = SAM.cigar(record)

		tflag = SAM.flag(record)
		flagresult = collect(last(bitstring(tflag), 12))
		if Int64(flagresult[6]) == 49
			rpair = "First"
		elseif Int64(flagresult[5]) == 49
			rpair = "Second"
		else
			println("Looks like its not paired")
			rpair = "Fail"
		end
		tpos = SAM.position(record)
		tcombo = string(tname, "-", rpair, "-", tpos)

		println(tcombo)
		if tcombo in keys(primaryd)
			println("Warning: Potentially multiple primary alignments!")
		else
			if SAM.isprimary(record)
				println("Adding")
				println(mapq)
				primaryd[tcombo] = [refn, mapq, cig]
			else
				println("Not a primary alignment.")
			end
		end
	end
end

# Find matching best hits among each read
# Report the read and alignment pair
reader = open(SAM.Reader, "test.sam")
# Count primary pairs
primarypairs = Dict{String,Int64}()

for record in reader
	# Only look at what mapped, always
	if SAM.ismapped(record)
		#Get the ID information
		tname = SAM.tempname(record)
		refn = SAM.refname(record)
		mapq = Int64(SAM.mappingquality(record))
		cig = SAM.cigar(record)

		tflag = SAM.flag(record)
		flagresult = collect(last(bitstring(tflag), 12))
		if Int64(flagresult[6]) == 49
			rpair = "First"
		elseif Int64(flagresult[5]) == 49
			rpair = "Second"
		else
			println("Looks like its not paired")
			rpair = "Fail"
		end
		tpos = SAM.position(record)
		tcombo = string(tname, "-", rpair, "-", tpos)

		# We want to count the read-pair combos because we only want
		# to keep those with two (a full paired match)
		combotr = string(tname, "-", refn, "-", tpos)
		# And establish the counts starting at zero
		if combotr in keys(primarypairs)
			#nothing
		else
			primarypairs[combotr] = 0
		end

		println("Searching ", tcombo, "-", refn)
		# Check if the read is a multimapper
		if counts[tcombo] > 1
			# Check if it is a secondary hit
			if SAM.isprimary(record)
				println(tcombo, "---", refn)
				# Keep track of primary alignment counts
				primarypairs[combotr] += 1
			else
				# Determine if this hit is the same as the
				# primary alignment of the read
				aref = primaryd[tcombo][1]
				aqual = primaryd[tcombo][2]
				acigar = primaryd[tcombo][3]
				if isequal(mapq, aqual)
					if isequal(cig, acigar)
						println(tcombo, "---", refn)
						primarypairs[combotr] += 1
					end
				end
			end
		elseif counts[tcombo] == 1
			println(tcombo, "---", refn)
		end
	end
end

# Print output to a 2 dim array for now

outarray = Array{Char}(undef, 0, 3)
for (key, value) in primarypairs
	if isequal(value, 2)
		karray = split(key, "-")
		tread = karray[1]
		tref = karray[2]
		spos = karray[3]
		println(karray)
		global outarray = [outarray; [tread tref spos]]
		# push!(outarray, karray)
	end
end

# Write the final array to tsv file for reference
writedlm( "testout.txt",  outarray, '\t')





