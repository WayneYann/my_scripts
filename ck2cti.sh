cp ~/Code/RMG-Py/examples/rmg/silane/chemkin/chem_annotated.inp .
cp ~/Code/RMG-Py/examples/rmg/silane/chemkin/tran.dat .
python /usr/local/lib/python2.7/site-packages/cantera/ck2cti.py --input=chem_annotated.inp --transport=tran.dat --output=silane.cti
