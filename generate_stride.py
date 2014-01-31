#! /usr/bin/python

import sys

def main(argv):

	cmdargs = str(sys.argv)

	#Converting start and end inputs to flexible residue position for ease of use when
	#calculating centres (see below)

	start = int(sys.argv[1])-1
	end = int(sys.argv[2])+1
	
	# Get flex residues and convert strings in list to integers
	flex = sys.argv[3:]
	flex = map(int, flex)

	stride = "STRIDE >"
	
	# A function to generate the STRIDE, given start, end and flexible residues
	def get_stride(start, end, flex, stride):
		
		for i in range(start+1,end):

			if i in flex:
				stride += "C"
			else:
				stride += "R"

		print stride
		return

	# A function to calculate centres between flexible residues i.e. segments
	def calc_centres(start,end,flex):

		if len(flex) == 0:

			len_seg = (end-1) - (start+1)
		
                        if len_seg % 2 == 0:

                                temp_centre = (len_seg/2)
	                        centre = temp_centre + start + 1
				
				print "The centre is %d" % centre
                       
			 # If length is odd then divide by 2, round up (ceiling) and add to start residue
                        else:

        	                temp_centre = int(round(len_seg/2))
                                centre = temp_centre + start + 1 

				print "The centre is %d" % centre


		if len(flex) == 1:

			tmp_flex = flex[0]
			len_seg1 = tmp_flex - start + 2
			len_seg2 = end - 1 - tmp_flex + 2
			temp_centre1 = len_seg1/2
			temp_centre2 = len_seg2/2
			centre1 = int(round(temp_centre1)) + start - 1 
			centre2 = int(round(temp_centre2)) + tmp_flex - 1
			
			print "Centres are %d and %d" % (centre1,centre2)

		if len(flex) > 1:
			
			flex_res = flex

			# Add start and end residues to flexible residues list for calculation of segment lengths
			flex_res.insert(0,start)
			flex_res.append(end)

			len_seg = int()
			temp_centre = int()
			centre = []

			# loop through flexible residues list to calculate segment lengths and centres
			for i in range (0,len(flex)-1):

				seg_start = flex_res[i] + 1
				seg_end = flex_res[i + 1] - 1

				#print "seg_start :%d" % seg_start
				#print "seg_end :%d" % seg_end

				# Calculate length of segments
				len_seg = (seg_end - seg_start + 1)		
				
				print "length of segment %d = %d" % (i+1,len_seg)

				# If length is even then divide by 2 and add to start residue
				if len_seg % 2 == 0:

					temp_centre = (len_seg/2)
					centre.append(temp_centre-1 + seg_start)
				
				# If length is odd then divide by 2, round up (ceiling) and add to start residue
				else:

					temp_centre = int(round(len_seg/2))
					centre.append(temp_centre + seg_start)

			print "Centres are: %s" % ', '.join(map(str, centre))

		return

	# Call functions
	get_stride(start,end,flex,stride)
	calc_centres(start,end,flex)	


	pass

if __name__ == "__main__":
	main(sys.argv)


