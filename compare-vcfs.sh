# Scaffold Dirs
mkdir -p output/{set1,set2,set3}/{bfc,bloocoo,lighter,nomodel,norealign}

# Run Python Scripts
python3 compare-vcf.py -s set1 -c bfc -o output/set1/bfc
python3 compare-vcf.py -s set2 -c bfc -o output/set2/bfc
python3 compare-vcf.py -s set3 -c bfc -o output/set3/bfc
python3 compare-vcf.py -s set1 -c nomodel -o output/set1/nomodel
python3 compare-vcf.py -s set2 -c nomodel -o output/set2/nomodel
python3 compare-vcf.py -s set3 -c nomodel -o output/set3/nomodel
python3 compare-vcf.py -s set1 -c norealign -o output/set1/norealign
python3 compare-vcf.py -s set2 -c norealign -o output/set2/norealign
python3 compare-vcf.py -s set3 -c norealign -o output/set3/norealign
python3 compare-vcf.py -s set1 -c bloocoo -o output/set1/bloocoo
python3 compare-vcf.py -s set2 -c bloocoo -o output/set2/bloocoo
python3 compare-vcf.py -s set3 -c bloocoo -o output/set3/bloocoo
python3 compare-vcf.py -s set1 -c lighter -o output/set1/lighter
python3 compare-vcf.py -s set2 -c lighter -o output/set2/lighter
python3 compare-vcf.py -s set3 -c lighter -o output/set3/lighter
