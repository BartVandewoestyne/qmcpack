#!/bin/bash

#### Arguments ####
ignorenodes="$1"

#### Templates ####
node_template="\$nodename [label=\\\"\$nodename\\\",color=\\\"\$nodecolor\\\",style=\"filled\" ];"
link_template="\$nodename -> \$toname;"

#### Title of graph ####
date=`date`
title="Dependencies for QMCPack ($date)"

#### .DOT-file header ####
echo "digraph G {"
echo "  compound=true;"
echo "          fontsize=\"12\";"
echo "          label=\"$title\";"
echo "          node [shape=box, fontsize=16, fontname=\"sans-serif\",edgecolor=black,decorate=true ];"
echo "          edge [arrowsize=1, color=black];"
echo "          concentrate=true;"
echo "          rankdir=\"LR\";"
echo "          size=\"210,297\";orientation=\"portrait\";"
echo "          ratio=\"compact\";scale=\"auto\";"
echo "          bgcolor=white;center=\"true\""
echo "  "
echo "	      overlap=\"voronoy\";"
echo "  "
echo "          /* Nodes */"

#### Parse all dependency lines ####
# expected line format:
# 	<target.o> : <dependency1> <dependency2> ... <dependencyn>

cluster=""
while read line
do
	
	# Extract the node (=target)
	nodename=${line%:*}
	nodename=${nodename%.*}
	if [ -n "$nodename" ]; then # not empty
		if [ "$ignorenodes" = "${ignorenodes%$nodename*}" ]; then # not ignored

			# Asign color
			nodecolor="lightgray"; 
			new_cluster="$cluster";
			if [ ${nodename:0:2} = "f_" ]; then
				nodecolor="greenyellow"
				new_cluster="Testfuncties"
			fi
			if [ ${nodename:0:4} = "mod_" ]; then
				nodecolor="lightblue"
				new_cluster="Modules"
			fi
			if [ "$nodename" = "qmcpack" ]; then
				nodecolor="lightsalmon"
			fi

			# check cluster
# 			if [ "$cluster" != "$new_cluster" ]; then
# 				if [ -n "$cluster" ]; then
# 					echo "} "
# 				fi
# 				cluster="$new_cluster"
# 				echo "subgraph cluster$cluster{"
# 				echo "label=\"$cluster\";"
# 				echo "color=\"$nodecolor\";"
# 			fi

			# Write node statement base on template
			eval "echo \"$node_template\""

			# Iterate over all links
			n_entries=`echo "$line" | awk '{ print NF }'`
			for((i=3;i<=n_entries;i++)); do
				# Get dependency <i>
				eval "toname=\`echo \"$line\" | awk '{ print \$$i }'\`"
				if [ "o" = "${toname#*.}" ]; then # only object files
					toname=${toname%.*}
					if [ "$ignorenodes" = "${ignorenodes%$toname*}" ]; then # not ignored
						# Write link statement base on template	
						eval "echo \"$link_template\""
					fi
				fi
			done
			echo " "
		fi
	fi
done

#### Close last cluster ####
# if [ -n "$cluster" ]; then
# 	echo "} "
# fi

#### .DOT-file footer
echo "  "
echo "}" 
