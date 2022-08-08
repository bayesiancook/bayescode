#!python3
import argparse
import pandas as pd
import numpy as np


def transform(s):
    try:
        x = float(s)
        if x < 0:
            return np.nan
        else:
            return np.log(x)
    except ValueError:
        return s


def clean_line(line):
    while "  " in line:
        line = line.replace("  ", " ")
    while "\t\t" in line:
        line = line.replace("\t\t", "\t")
    return line.replace(" ", "\t").strip()


def main(args):
    r = [clean_line(line) for line in open(args.input, 'r').readlines()]
    assert r[0].strip() == '#TRAITS'
    n_rows = int(r[1].split('\t')[0])
    n_cols = int(r[1].split('\t')[1])
    print(f"Found {n_rows} taxa and {n_cols} traits.")
    columns = r[1].split('\t')[2:]
    assert len(columns) == n_cols
    table = [[transform(s) for s in line.split('\t')] for line in r[2:]]
    assert len(table) == n_rows

    df = pd.DataFrame(table, columns=["TaxonName"] + columns)
    assert len(df) == n_rows
    assert len(df.columns) == n_cols + 1
    df.to_csv(args.output, index=False, na_rep="NaN", sep="\t")
    print(f"Converted file '{args.input}' to mutsel format (in log-space).")
    print(f"Output file is '{args.output}'.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input', required=True, type=str, dest="input", help="Input tsv file")
    parser.add_argument('--output', required=True, type=str, dest="output", help="Output tsv file")
    main(parser.parse_args())
