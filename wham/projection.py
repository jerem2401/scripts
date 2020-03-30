#!/usr/bin/env python


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='give 2dpmf_clean.txt file', action='store', dest='f', type=str)
    parser.add_argument('-o', help='give path for plot output', action='store', dest='o')
    args = parser.parse_args()

    df = pandas.read_table(args.f, delim_whitespace=True, dtype={'#X': np.float64, 'Y': np.float64, 'Free': np.float64, 'Pro': np.float64})
    df['Pro2'] = np.exp((df['Free'])/-2.494)
    aggregation_functions = {'Pro2': 'sum'}
    df_new = df.groupby(df['#X']).aggregate(aggregation_functions)
    df_new=df_new.reset_index()

