import sqlite3
import csv

def dump(version, con):
    # determine filename based on version/treeinform run
    if version == 0:
        filename='trinity.csv'
    else:
        filename='treeinform%s.csv' % version
    outfile = open(filename, 'wb')
    outcsv = csv.writer(outfile)

    cursor = con.execute('SELECT model_id,gene,isoform FROM agalma_genes WHERE version=%s AND run_id=6;' % version)

    # get and dump column names
    colnames = []
    for x in cursor.description:
        colnames.append(x[0])
    outcsv.writerow(colnames)

    outcsv.writerows(cursor.fetchall())
    outfile.close()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("database", help="agalma sqlite database")
    parser.add_argument("versions", nargs='*', type=int, help="version numbers of original trinity run (0) and additional treeinform runs")

    args = parser.parse_args()

    con = sqlite3.connect(args.database)
    print(args.versions)

    for i in args.versions:
        dump(i, con)
