class Alignment:

    def __init__(self, protein_seq, peptides_metadata):
        self.protein_seq = protein_seq
        self.peptides_metadata = peptides_metadata
        self.amount_of_the_same_peptides = {}

    def create_lps(self, peptide_seq: str, peptide_seq_len: int) -> list[int]:
        """
        Create lps table.

        :param peptide_seq: peptide sequence
        :param peptide_seq_len: length of peptide sequence
        :return: list containing integers which store information
        about length of prefix
        """
        lps = [0] * peptide_seq_len
        prefix_len = 0
        i = 1

        while i < peptide_seq_len:
            if peptide_seq[prefix_len] == peptide_seq[i]:
                prefix_len += 1
                lps[i] = prefix_len
                i += 1
            elif prefix_len != 0:
                prefix_len = lps[prefix_len - 1]
            else:
                lps[i] = 0
                i += 1

        return lps

    def search_peptide_in_protein_seq(self) -> None:
        """
        KMP (Knuth–Morris–Pratt algorithm) for pattern searching.

        """
        # index is not used by neccesacry to iter thorugh pepitde
        protein_seq_len = len(self.protein_seq)
        for index, peptide in self.peptides_metadata.iterrows():
            peptide_seq_len = len(peptide.Sequence)

            lps = self.create_lps(peptide.Sequence, peptide_seq_len)
            i = 0
            j = 0

            while (protein_seq_len - i) >= (peptide_seq_len - j):
                if j == peptide_seq_len:
                    coords = (i - j, i)
                    #TODO: think about name, maybe start and end position or maybe only start and in program add lenght to it?
                    peptide["Coords"] = coords
                    if peptide["Sequence"] in self.amount_of_the_same_peptides:
                        self.amount_of_the_same_peptides[peptide["Sequence"]] += 1
                    else:
                        self.amount_of_the_same_peptides[peptide["Sequence"]] = 1
                    j = lps[j - 1]
                elif self.protein_seq[i] == peptide.Sequence[j]:
                    i += 1
                    j += 1
                elif j > 0:
                    j = lps[j - 1]
                else:
                    i += 1
        print(self.peptides_metadata)
        #self.peptides_metadata.sort(key=lambda x: x["Coords"])

