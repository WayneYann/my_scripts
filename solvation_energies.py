"""
Perform single point energy calculations in solvent on a set of gas phase
reactant and TS geometries in a database.
"""
import sqlite3

conn = sqlite3.connect('/scratch/westgroup/ts_data.db')
print "Opened database successfully"

geometries = conn.execute("SELECT ")

for g in geometries:
    # Do single point energy calc
