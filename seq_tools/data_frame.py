from seq_tools import sequence


def get_nucleic_acid_type(df):
    types = {"DNA": 0, "RNA": 0}
    for i, row in df.iterrows():
        types[sequence.get_nucleic_acid_type(row["sequence"])] += 1
    if types["DNA"] > types["RNA"]:
        return "DNA"
    else:
        return "RNA"


def convert_to_dna(df):
    seqs = []
    for i, row in df.iterrows():
        seqs.append(sequence.convert_to_dna(row["sequence"]))
    df["sequence"] = seqs


def convert_to_rna(df):
    seqs = []
    for i, row in df.iterrows():
        seqs.append(sequence.convert_to_rna(row["sequence"]))
    df["sequence"] = seqs


def get_reverse_complement(df, type):
    rcs = []
    for i, row in df.iterrows():
        rcs.append(sequence.get_reverse_complement(row["sequence"], type))
    df["rev_comp"] = rcs


def get_length(df):
    lens = []
    for i, row in df.iterrows():
        lens.append(len(row["sequence"]))
    df["len"] = lens


def get_molecular_weight(df, type, ds):
    mws = []
    for i, row in df.iterrows():
        mws.append(round(sequence.get_molecular_weight(row["sequence"], type, ds)))
    df["molecular weight"] = mws


def trim_5p(df, length):
    seqs = []
    ss = []
    include_ss = False
    if df.loc[1]["structure"] is not None:
        include_ss = True
    for i, row in df.iterrows():
        seqs.append(row["sequence"][length:])
        if include_ss:
            ss.append(row["structure"][length:])
    df["sequence"] = seqs
    if include_ss:
        df["structure"] = ss


def trim_3p(df, length):
    seqs = []
    ss = []
    include_ss = False
    if df.loc[1]["structure"] is not None:
        include_ss = True
    for i, row in df.iterrows():
        seqs.append(row["sequence"][:-length])
        if include_ss:
            ss.append(row["structure"][:-length])
    df["sequence"] = seqs
    if include_ss:
        df["structure"] = ss


def add_5p(df, p5):
    spl = p5.split(",")
    seq_p5 = p5
    ss_p5 = ""
    include_ss = False
    if len(spl) == 2:
        include_ss = True
        seq_p5 = p5[0]
        ss_p5 = p5[1]
    seqs = []
    ss = []
    if include_ss:
        if "structure" not in df.columns:
            print("cannot add to sequence and structure, structure is not defined!")
            exit()
        if df.loc[1]["structure"] is None:
            print("cannot add to sequence and structure, structure is not defined!")
            exit()
    for i, row in df.iterrows():
        seqs.append(seq_p5 + row["sequence"])
        if include_ss:
            ss.append(ss_p5 + row["structure"])
    df["sequence"] = seqs
    if include_ss:
        df["structure"] = ss


def add_3p(df, p3):
    spl = p3.split(",")
    seq_p3 = p3
    ss_p3 = ""
    include_ss = False
    if len(spl) == 2:
        include_ss = True
        seq_p3 = p3[0]
        ss_p3 = p3[1]
    seqs = []
    ss = []
    if include_ss:
        if "structure" not in df.columns:
            print("cannot add to sequence and structure, structure is not defined!")
            exit()
        if df.loc[1]["structure"] is None:
            print("cannot add to sequence and structure, structure is not defined!")
            exit()
    for i, row in df.iterrows():
        seqs.append(row["sequence"] + seq_p3)
        if include_ss:
            ss.append(row["structure"] + ss_p3)
    df["sequence"] = seqs
    if include_ss:
        df["structure"] = ss


def remove_t7(df):
    t7 = "TTCTAATACGACTCACTATA"
    seqs = []
    for i, row in df.iterrows():
        if row["sequence"].find(t7) == -1:
            print("sequence does not contain the t7 promoter cannot remove it!")
            seqs.append(row["sequence"])
            continue
        seqs.append(row["sequence"][20:])
    df["sequence"] = seqs

