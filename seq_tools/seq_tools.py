import click
import pandas as pd
import numpy as np
from seq_tools import dataframe


def is_seq_or_csv(input):
    input = input.upper()
    for e in input:
        if e == "A" or e == "C" or e == "G" or e == "U" or e == "T" or e == "N":
            continue
        else:
            return "CSV"
    return "SEQ"


def get_df_from_seq_and_ss(seq, ss, name):
    if name is None:
        name = "seq"
    df = pd.DataFrame(columns="name,sequence,structure".split(","))
    df.loc[0] = [name, seq, ss]
    return df


def process_args(p):
    if p["type"] is not None:
        p["type"] = p["type"].upper()

    input_type = is_seq_or_csv(p["input"])
    if input_type == "CSV":
        df = pd.read_csv(p["input"])
    else:
        df = get_df_from_seq_and_ss(p["input"], p["ss"], p["name"])

    # ensure sequence is all upper case
    df['sequence'] = [s.upper() for s in df['sequence']]

    if "name" not in df:
        names = []
        for i, row in df.iterrows():
            names.append(f"seq_{i}")
        df.insert(0, "name", names)
    type = data_frame.get_nucleic_acid_type(df)
    if p["remove_t7"]:
        data_frame.remove_t7(df)
    if p["type"] == "RNA":
        data_frame.convert_to_rna(df)
        type = "RNA"
    if p["type"] == "DNA":
        data_frame.convert_to_dna(df)
        type = "DNA"
    if p["add_t7"]:
        if type != "DNA":
            raise ValueError("cannot add t7 to RNA only to DNA")
        data_frame.add_5p(df, "TTCTAATACGACTCACTATA")
    if "structure" not in df.columns and type == "RNA":
        df["structure"] = None
    if p["trim5"] is not None:
        data_frame.trim_5p(df, p["trim5"])
    if p["trim3"] is not None:
        data_frame.trim_3p(df, p["trim3"])
    if p["fold"]:
        if type == "DNA":
            raise ValueError("cannot fold DNA, use -t/--type RNA to convert")
        data_frame.get_folded_structure(df)

    ds = p["ds"]
    if ds is None:
        if type == "DNA":
            ds = True
        else:
            ds = False
    calcs = {}
    if p["calc"] is not None:
        spl = p["calc"].split(",")
        for e in spl:
            calcs[e] = True

    if "l" in calcs or "len" in calcs or "length" in calcs:
        data_frame.get_length(df)
    if "mw" in calcs:
        data_frame.get_molecular_weight(df, type, ds)
    if "ec" in calcs:
        data_frame.get_extinction_coeff(df, type, ds)
    if "rc" in calcs:
        data_frame.get_reverse_complement(df, type)

    return df


def run_main(p):
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
    print(f"outputing results to csv: {p['output']}")
    df.to_csv(p["output"], index=False)
    if len(df) == 1:
        df = df.drop(["name"], axis=1)
        d = df.loc[0].to_dict()
        for k, v in d.items():
            if v == False or v == None:
                continue
            print(f"{k} : {v}")
    else:
        if "structure" in df.columns:
            if df.loc[0]["structure"] is None:
                df = df.drop("structure", axis=1)
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


@click.command("a simple command for manipulating rna and dna sequences in dataframes")
@click.argument("input")
@click.option("-n", "--name", required=False, default=None)
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
@click.option("-ds", required=False, is_flag=True, default=None)
@click.option("-fold", required=False, is_flag=True, default=False)
@click.option("-avg", required=False, is_flag=True, default=False)
def main(**p):
    run_main(p)


if __name__ == "__main__":
    main()
