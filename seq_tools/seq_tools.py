import click
import pandas as pd
import numpy as np
from seq_tools import data_frame


def is_seq_or_csv(input):
    input = input.upper()
    for e in input:
        if e == 'A' or e == 'C' or e == 'G' or e == 'U' or e == 'T' or e == 'N':
            continue
        else:
            return "CSV"
    return "SEQ"


def get_df_from_seq_and_ss(seq, ss):
    df = pd.DataFrame(columns="name,sequence,structure".split(","))
    df.loc[1] = ["seq", seq, ss]
    return df


def process_args(p):
    input_type = is_seq_or_csv(p["input"])
    if input_type == "CSV":
        df = pd.read_csv(p["input"])
    else:
        df = get_df_from_seq_and_ss(p["input"], p["ss"])

    calcs = {}
    if p["calc"] is not None:
        spl = p["calc"].split(",")
        for e in spl:
            calcs[e] = True

    if "l" in calcs or "len" in calcs or "length" in calcs:
        data_frame.get_length(df)

    return df


@click.command('a simple command for manipulating rna and dna sequences in dataframes')
@click.argument("input")
@click.option("-ss", required=False, default=None)
@click.option("-c", "--calc", required=False, default=None)
@click.option("-t", "--type", required=False, default=None)
@click.option("-trim5", required=False, default=None)
@click.option("-trim3", required=False, default=None)
@click.option("-add5", required=False, default=None)
@click.option("-add3", required=False, default=None)
@click.option("-o", "--output", required=False, default="output.csv")
@click.option("-to_fasta", required=False, is_flag=True, default=False)
@click.option("-remove_t7", required=False, is_flag=True, default=False)
@click.option("-add_t7", required=False, is_flag=True, default=False)
@click.option("-ds", required=False, is_flag=True, default=False)
@click.option("-fold", required=False, is_flag=True, default=False)
@click.option("-avg", required=False, is_flag=True, default=False)
def main(**p):
    if p["to_rna"] and p["to_dna"]:
        print("cannot use flag -to_rna and -to_dna")
        exit()
    if p["remove_t7"] and p["add_t7"]:
        print("cannot use flag -remove_t7 and -add_t7")
        exit()
    if p["add5"] and p["trim5"]:
        print("cannot use flag -trim5 and -add5")
        exit()
    if p["add3"] and p["trim3"]:
        print("cannot use flag -trim3 and -add3")
        exit()

    df = process_args(p)
    if p["to_fasta"]:
        f = open("test.fasta", "w")
        for i, row in df.iterrows():
            f.write(f">{row['name']}\n")
            f.write(row["sequence"] + "\n")
        f.close()
    if len(df) == 1:
        df = df.drop(["name"], axis=1)
        d = df.loc[1].to_dict()
        for k, v in d.items():
            if v == False or v == None:
                continue
            print(f"{k} : {v}")
    else:
        if "structure" is df.columns:
            if df.loc[1]["structure"] is None:
                df = df.drop("structure", axis=1)
        print("outputing results to csv: output.csv")
        df.to_csv("output.csv", index=False)
        if p["avg"]:
            excess = "name,sequence,structure,rc".split(",")
            for c in excess:
                if c in df.columns:
                    df = df.drop(c, axis=1)
            cols = df.columns
            print("averages: ")
            for c in cols:
                try:
                    avg = np.mean(df[c])
                    print(c + " : " + str(round(avg, 2)))
                except:
                    pass


if __name__ == "__main__":
    main()
