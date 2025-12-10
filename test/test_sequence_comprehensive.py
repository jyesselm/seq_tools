"""
Comprehensive tests for sequence module to improve coverage
"""

from seq_tools.sequence import (
    get_max_stretch,
    get_molecular_weight,
    get_reverse_complement,
    to_dna,
    to_dna_template,
    to_rna,
)


class TestGetMaxStretch:
    """Tests for get_max_stretch function."""

    def test_get_max_stretch_all_same(self):
        """Test max stretch with all same nucleotides."""
        assert get_max_stretch("AAAAAAA") == 7

    def test_get_max_stretch_all_different(self):
        """Test max stretch with all different nucleotides."""
        assert get_max_stretch("ATCGATCG") == 1

    def test_get_max_stretch_in_middle(self):
        """Test max stretch in middle of sequence."""
        assert get_max_stretch("ATCCCCGG") == 4

    def test_get_max_stretch_at_start(self):
        """Test max stretch at start of sequence."""
        assert get_max_stretch("AAAATCG") == 4

    def test_get_max_stretch_at_end(self):
        """Test max stretch at end of sequence."""
        assert get_max_stretch("ATCGGGG") == 4

    def test_get_max_stretch_single_nucleotide(self):
        """Test max stretch with single nucleotide."""
        assert get_max_stretch("A") == 1

    def test_get_max_stretch_multiple_stretches(self):
        """Test max stretch with multiple stretches."""
        assert get_max_stretch("AAATTCCCCGG") == 4


class TestGetMolecularWeight:
    """Tests for get_molecular_weight function."""

    def test_get_molecular_weight_dna_single_stranded(self):
        """Test DNA molecular weight single stranded."""
        mw = get_molecular_weight("ATCG", "DNA", False)
        assert isinstance(mw, float)
        assert mw > 0

    def test_get_molecular_weight_dna_double_stranded(self):
        """Test DNA molecular weight double stranded."""
        mw_single = get_molecular_weight("ATCG", "DNA", False)
        mw_double = get_molecular_weight("ATCG", "DNA", True)
        # Double stranded should be about twice single stranded
        assert mw_double > mw_single

    def test_get_molecular_weight_rna_single_stranded(self):
        """Test RNA molecular weight single stranded."""
        mw = get_molecular_weight("AUCG", "RNA", False)
        assert isinstance(mw, float)
        assert mw > 0

    def test_get_molecular_weight_rna_double_stranded(self):
        """Test RNA molecular weight double stranded."""
        mw_single = get_molecular_weight("AUCG", "RNA", False)
        mw_double = get_molecular_weight("AUCG", "RNA", True)
        assert mw_double > mw_single

    def test_get_molecular_weight_auto_convert(self):
        """Test that function auto-converts to correct type."""
        # Pass RNA sequence with DNA type - should convert
        mw = get_molecular_weight("AUCG", "DNA", False)
        assert isinstance(mw, float)


class TestGetReverseComplement:
    """Tests for get_reverse_complement function."""

    def test_get_reverse_complement_dna(self):
        """Test DNA reverse complement."""
        assert get_reverse_complement("ATCG", "DNA") == "CGAT"

    def test_get_reverse_complement_rna(self):
        """Test RNA reverse complement."""
        assert get_reverse_complement("AUCG", "RNA") == "CGAU"

    def test_get_reverse_complement_palindrome(self):
        """Test reverse complement of palindrome."""
        # ATCGAT reverse complement is ATCGAT
        rc = get_reverse_complement("ATCGAT", "DNA")
        assert rc == "ATCGAT"

    def test_get_reverse_complement_auto_convert(self):
        """Test auto-conversion in reverse complement."""
        # Pass DNA sequence with RNA type - should convert
        rc = get_reverse_complement("ATCG", "RNA")
        assert "U" in rc
        assert "T" not in rc


class TestToDna:
    """Tests for to_dna function."""

    def test_to_dna_basic(self):
        """Test basic RNA to DNA conversion."""
        assert to_dna("AUCG") == "ATCG"

    def test_to_dna_already_dna(self):
        """Test converting DNA to DNA."""
        assert to_dna("ATCG") == "ATCG"

    def test_to_dna_long_sequence(self):
        """Test conversion of long sequence."""
        rna = "AUCGAUCGAUCG"
        dna = to_dna(rna)
        assert "U" not in dna
        assert dna.count("T") == rna.count("U")


class TestToDnaTemplate:
    """Tests for to_dna_template function."""

    def test_to_dna_template_basic(self):
        """Test basic DNA template generation."""
        result = to_dna_template("AUCG")
        assert result.startswith("TTCTAATACGACTCACTATA")
        assert result.endswith("ATCG")

    def test_to_dna_template_long_sequence(self):
        """Test DNA template with long sequence."""
        rna = "AUCGAUCGAUCGAUCG"
        template = to_dna_template(rna)
        assert template.startswith("TTCTAATACGACTCACTATA")


class TestToRna:
    """Tests for to_rna function."""

    def test_to_rna_basic(self):
        """Test basic DNA to RNA conversion."""
        assert to_rna("ATCG") == "AUCG"

    def test_to_rna_already_rna(self):
        """Test converting RNA to RNA."""
        assert to_rna("AUCG") == "AUCG"

    def test_to_rna_long_sequence(self):
        """Test conversion of long sequence."""
        dna = "ATCGATCGATCG"
        rna = to_rna(dna)
        assert "T" not in rna
        assert rna.count("U") == dna.count("T")
