# Author Chao Zhang

global R = [Dict("k"=>5, "p"=>0.25, "q"=>0.1, "L"=>30),
     Dict("k"=>9, "p"=>0.25, "q"=>0.25, "L"=>54),
     Dict("k"=>17, "p"=>0.1, "q"=>0.5, "L"=>Inf)]

PROGRAM_VERSION = v"1.0.0"
try
	using ArgParse
catch
	import Pkg
	Pkg.add("ArgParse")
	using ArgParse
end
import Statistics.median

function correction(c, output, k, X, MASK, pvalue, qvalue, threshold)
	upperx = ('a' <= X && X <= 'z') ? X - 'a' + 'A' : X
	upperc = [uppercase(str) for str in c]
	arrc = [Array{Char, 1}(str) for str in c]
	n = length(c[1])
	m = length(c)
	wo = zeros(n, m)
	# read sequences and compute per column profiles.
	for i in 1:n
		cnt = zeros(128)
		# read a column, internally represent in upper case, and count letters
		for j in 1:m
			cnt[UInt8(upperc[j][i])] += 1
		end
		cnt[UInt8(upperx)] = 0
		cnt[UInt8('-')] = 0
		unq = length([1 for t in cnt if t > 0])
		total = sum(cnt)
		for j in 1:m
			wo[i, j] = total == 0 ? 0 : total / (unq * cnt[UInt8(upperc[j][i])])
		end
	end
	w1 = [wo[:,j][(arrc[j] .!= '-') .& (arrc[j] .!= X)] for j in 1:m]
	w = [[median(arr[i:i+k-1]) for i in 1:length(arr)-k+1] for arr in w1]
	ws = [[sum(arr[i:i+k-1]) for i in 1:length(arr)-k+1] for arr in w1]
	wsorted = [sort(arr) for arr in w]
	wsum = [accumulate(+, arr) for arr in wsorted]
	f(x, y, n, m) = x^2 / n + (n == m ? 0 : (y - x)^2 / (m - n))
	var = [length(arr) > 0 ? f.(arr, arr[end], 1:length(arr), length(arr)) : [] for arr in wsum]
	wCutoff = [length(var[j]) > 0 ? wsorted[j][findmax(var[j])[2]] : 0 for j in 1:m]
	cutoffSorted = sort([wCutoff[j] for j in 1:m if length(var[j]) > 0])
	cutoffFloor = (pvalue >= 1) ? 0 : cutoffSorted[end - floor(Int, length(cutoffSorted) * pvalue)]
	s = zeros(n - k + 1, m, 2)
	tiebreaker = zeros(n - k + 1, m, 2)
	bt = zeros(Int64, n - k + 1, m, 2)
	for j in 1:m
		wj = w[j]
		wsj = ws[j]
		L = length(wj)
		if L <= k
			push!(output, c[j])
			continue
		end
		s = zeros(L, 2)
		tiebreaker = zeros(L, 2)
		bt = zeros(Int64, L, 2)
		cutoff = max(wCutoff[j], cutoffFloor, (qvalue >= 1) ? 0 : wsorted[j][end - floor(Int, length(wsorted[j]) * qvalue)], threshold)
		for i in 1:L
			v = (wj[i] > cutoff ? 0 : 1)
			if i == 1
				s[i, 1] = v
				s[i, 2] = 1 - v
				tiebreaker[i, 1] = 0
				tiebreaker[i, 2] = wsj[i]
			else
				s[i, 1] = s[i - 1, 1] + v
				s[i, 2] = s[i - 1, 2] + 1 - v
				tiebreaker[i, 1] = tiebreaker[i - 1, 1]
				tiebreaker[i, 2] = tiebreaker[i - 1, 2] + wsj[i]
				bt[i, 1] = 1
				bt[i, 2] = 2
			end
			if i > k && (s[i, 1], tiebreaker[i, 1]) < (s[i - k, 2] + v, tiebreaker[i - k, 2])
				s[i, 1] = s[i - k, 2] + v
				tiebreaker[i, 1] = tiebreaker[i - k, 2]
				bt[i, 1] = 2
			end
			if i > k && (s[i, 2], tiebreaker[i, 2]) < (s[i - k, 1] + 1 - v, tiebreaker[i - k, 1])
				s[i, 2] = s[i - k, 1] + 1 - v
				tiebreaker[i, 2] = tiebreaker[i - k, 1] + wsj[i]
				bt[i, 2] = 1
			end
		end
		str = arrc[j][(arrc[j] .!= '-') .& (arrc[j] .!= X)]
		icur = L
		if s[L, 1] < s[L, 2]
			str[L:L+k-1] .= MASK
			bcur = 2
		else
			bcur = 1
		end
		while true
			if bcur == 1 && bt[icur, bcur] == 1
				icur -= 1
				bcur = 1
			elseif bcur == 1 && bt[icur, bcur] == 2
				icur -= k
				bcur = 2
				str[icur : icur + k - 1] .= MASK
			elseif bcur == 2 && bt[icur, bcur] == 1
				icur -= k
				bcur = 1
			elseif bcur == 2 && bt[icur, bcur] == 2
				icur -= 1
				bcur = 2
				str[icur] = MASK
			elseif bcur == 1
				break
			else
				str[1:icur - 1] .= MASK
				break
			end
		end
		strout = ""
		i = 1
		for t in 1:length(c[j])
			if c[j][t] == X || c[j][t] == '-'
				strout *= c[j][t]
			else
				strout *= str[i]
				i += 1
			end
		end
		push!(output, strout)
	end
end

function ArgParse.parse_item(::Type{Char}, x::AbstractString)
	return x[1]
end

function parse_commandline()
	s = ArgParseSettings()

	@add_arg_table! s begin
		"--list", "-l"
			help = "running on a list of inputs; for every two lines of the list file, the first one should be the path to the input and the second should be the path to its output"
			action = :store_true
		"--mask", "-m"
			help = "the character to mask erroneous regions"
			arg_type = Char
			default = 'X'
		"--any", "-a"
			help = "the character to denote ambiguous positions or character to denote ANY in the input files"
			arg_type = Char
			default = 'X'
		"--cutoff", "-c"
			help = "set score cutoff to control the minimum aggressiveness of masking (should be > 1)"
			arg_type = Float64
			default = 3.0
		"--parameter", "-p"
			help = "(Advanced) load the list of k, p, q, and L from the input parameter file"
			arg_type = String
		"input"
			help = "a fasta file as input"
			required = true
	end

	return parse_args(s)
end

function correction_multi(args, fin, fout)
	inputText = read(fin, String)
	temp = split(inputText, ">")
	temp = temp[length.(temp) .> 0]
	temp = [split(arr, "\n") for arr in temp]
	header = [arr[1] for arr in temp]
	c = [string(arr[2:end]...) for arr in temp]
	arrc = [Array{Char, 1}(str) for str in c]
	n = length(c[1])
	m = length(c)
	MASK = args["mask"]
	for r in R
		output = []
		correction(c, output, get(r, "k", 9), args["any"], '*', get(r, "p", 0), get(r, "q", 0), args["cutoff"])
		L = get(r, "L", Inf)
		for j in 1:m
			start = 0
			cnt = 0
			for i in 1:n
				if start == 0 && output[j][i] == '*'
					start = i
					cnt = 1
				end
				if start != 0 && output[j][i] == '*'
					cnt += 1
				end
				if start != 0 && output[j][i] != '*' && output[j][i] != '-'
					if cnt < L
						for k in start:i-1
							arrc[j][k] = (output[j][k] == '-') ? '-' : MASK
						end
					end
					start = 0
					cnt = 0
				end
			end
			if start != 0 && cnt < L
				for k in start:n
						arrc[j][k] = (output[j][k] == '-') ? '-' : MASK
				end
			end
		end
	end
	for j in 1:m
		println(fout, ">" * header[j])
		println(fout, String(arrc[j]))
	end
end

function main()
	println(stderr, "Version " * string(PROGRAM_VERSION))
	args = parse_commandline()
	if isnothing(args["parameter"]) == false
		global R = eval(Meta.parse(read(open(args["parameter"], "r"), String)))
	end
	if args["list"] == false
		correction_multi(args, open(args["input"], "r"), stdout)
	else
		temp = split(read(open(args["input"], "r"), String), "\n")
		temp = temp[length.(temp) .> 0]
		for i = 2:2:length(temp)
			try
				println(stderr, "Processing " * temp[i - 1] * "...")
				correction_multi(args, open(temp[i - 1], "r"), open(temp[i], "w"))
			catch
				println(stderr, "Error happened when processing " * temp[i - 1] * ".")
				println(stderr)
			end
		end
	end
end

main()