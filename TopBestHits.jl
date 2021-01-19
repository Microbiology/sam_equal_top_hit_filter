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
reader = open(SAM.Reader, "test1.sam")
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
			println("Looks like its not paired")
			rpair = "Fail"
		end

		# For the multimappers, we only want to look at the sequence read name
		# and whether it if first or second. Location doesn't matter.

		tcombo = string(tname, "-", rpair)

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
# Should only be one primary alignment per read
# Key is the read, value is the reference
reader = open(SAM.Reader, "test1.sam")
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
			# Make sure to include position for better ID
			tpos = SAM.position(record)
			trpos = SAM.nextposition(record)
		elseif Int64(flagresult[5]) == 49
			rpair = "Second"
			# Reverse for the other paired read
			trpos = SAM.position(record)
			tpos = SAM.nextposition(record)
		else
			println("Looks like its not paired")
			rpair = "Fail"
		end

		tcombo = string(tname, "-", rpair)

		println(tcombo)
		if tcombo in keys(primaryd)
			println("Warning: Potentially multiple primary alignments!")
		else
			if SAM.isprimary(record)
				println("Adding")
				println(mapq)
				primaryd[tcombo] = [refn, tpos, trpos, mapq, cig]
			else
				println("Not a primary alignment.")
			end
		end
	end
end

# Find matching best hits among each read
# Report the read and alignment pair
reader = open(SAM.Reader, "test1.sam")
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
			# Make sure to include position for better ID
			tpos = SAM.position(record)
			trpos = SAM.nextposition(record)
		elseif Int64(flagresult[5]) == 49
			rpair = "Second"
			# Reverse for the other paired read
			trpos = SAM.position(record)
			tpos = SAM.nextposition(record)
		else
			println("Looks like its not paired")
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

		println("Searching ", tcombo, "-", refn)
		println("Multimapper ID is ", multimapperid)
		# Check if the read is a multimapper
		if counts[multimapperid] > 1
			# Check if it is a secondary hit
			if SAM.isprimary(record)
				println(tcombo, "---", refn)
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
							println(tcombo, "---", refn)
							primarypairs[combotr] += 1
						end
					end
				end
			end
		elseif counts[multimapperid] == 1
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




