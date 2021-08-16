for f in *.csv; do
	mv -- "$f" "${f%.csv}.txt"
done

# To run the script
	# 1. copy the batch_change_file_extension.sh file to the directory
	# 2. run the following command:
		# bash batch_change_file_extension.sh