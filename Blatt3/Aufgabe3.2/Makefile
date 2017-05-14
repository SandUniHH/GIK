
test: vienna2pairlist.py
	./vienna2pairlist.py Vienna-examples.txt > pairlist.txt
	./pairlist2vienna.py pairlist.txt | diff - Vienna-examples.txt

pairlist: vienna2pairlist.py
	./vienna2pairlist.py Vienna-examples.txt > pairlist.txt

vienna:
	./pairlist2vienna.py pairlist.txt