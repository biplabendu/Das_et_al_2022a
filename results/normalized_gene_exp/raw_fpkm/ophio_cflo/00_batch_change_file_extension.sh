for f in *.csv; do
	mv -- "$f" "${f%.csv}.txt"
done
