from seq_tools import sequence, extinction_coeff
import vienna


def add(df, p5_seq, p3_seq) -> None:
    """
    adds a 5' and 3' sequence to the sequences in the dataframe
    :param df: dataframe
    :param p5_seq: 5' sequence
    :param p3_seq: 3' sequence
    :return: None
    """
    df["sequence"] = df["sequence"].apply(lambda x: p5_seq + x + p3_seq)
    if "structure" in df.columns:
        pass


def to_dna(df) -> None:
    """
    converts each sequence in dataframe to DNA
    :return: None
    """
    df["sequence"] = df["sequence"].apply(sequence.to_dna)



def get_molecular_weight(df, ntype, double_stranded) -> None:
    """
    :param df: pandas data frame
    :param ntype: nucleotide type, RNA or DNA
    :param double_stranded: is double stranded?
    :return: None
    """
    df["mw"] = df["sequence"].apply(
        lambda x: sequence.get_molecular_weight(x, ntype, double_stranded)
    )


def trim(df, p5_length, p3_length) -> None:
    """
    takes a data frame and trims the sequences. If there is a structure
    it will also trim the structure
    :param df: dataframe
    :param p5_length: length to trim from 5'
    :param p3_length: length to trim from 3'
    :return: None
    """
    # trim `sequence` column and `structure` column
    df["sequence"] = df["sequence"].str.slice(p5_length, p3_length)
    if "structure" in df.columns:
        df["structure"] = df["structure"].str.slice(p5_length, p3_length)



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


def get_extinction_coeff(df, type, ds):
    ecs = []
    for i, row in df.iterrows():
        if type == "RNA":
            ec = extinction_coeff.get_coefficient_rna(row["sequence"], row["structure"])
        else:
            ec = extinction_coeff.get_coefficient_dna(row["sequence"], ds)
        ecs.append(ec)
    df["extinction coeff"] = ecs


def get_folded_structure(df):
    structures = []
    mfes = []
    ensemble_diversities = []
    for i, row in df.iterrows():
        vr = vienna.fold(row["sequence"])
        structures.append(vr.dot_bracket)
        mfes.append(vr.mfe)
        ensemble_diversities.append(vr.ens_defect)
    df["structure"] = structures
    df["mfe"] = mfes
    df["ens_defect"] = ensemble_diversities
