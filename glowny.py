from brute import brute_force_motif_search
from mój_brute import m_brute_force_motif_search
from greedy import greedy_motif_search
from randomized import randomized_motif_search


dna = ['tagtggtcttttgagtgtagatctgaagggaaagtatttccaccagttcggggtcacccagcagggcagggtgacttaat',
'cgcgactcggcgctcacagttatcgcacgtttagaccaaaacggagttggatccgaaactggagtttaatcggagtcctt',
'gttacttgtgagcctggttagacccgaaatataattgttggctgcatagcggagctgacatacgagtaggggaaatgcgt',
'aacatcaggctttgattaaacaatttaagcacgtaaatccgaattgacctgatgacaatacggaacatgccggctccggg',
'accaccggataggctgcttattaggtccaaaaggtagtatcgtaataatggctcagccatgtcaatgtgcggcattccac',
'tagattcgaatcgatcgtgtttctccctctgtgggttaacgaggggtccgaccttgctcgcatgtgccgaacttgtaccc',
'gaaatggttcggtgcgatatcaggccgttctcttaacttggcggtgcagatccgaacgtctctggaggggtcgtgcgcta',
'atgtatactagacattctaacgctcgcttattggcggagaccatttgctccactacaagaggctactgtgtagatccgta',
'ttcttacacccttctttagatccaaacctgttggcgccatcttcttttcgagtccttgtacctccatttgctctgatgac',
'ctacctatgtaaaacaacatctactaacgtagtccggtctttcctgatctgccctaacctacaggtcgatccgaaattcg']


if __name__ == "__main__":
    with open("wyniki_E.txt", 'w') as plik:
        plik.write("Brute Force(1): "+str(brute_force_motif_search(dna, len(dna), len(dna[1]), 78)))
        plik.write("\nBrute Force(2): " + str(m_brute_force_motif_search(dna, len(dna), len(dna[1]), 78)))
        plik.write("\nGreedy: "+str(greedy_motif_search(dna, 78, len(dna))))
        plik.write("\nRandomized: "+str(randomized_motif_search(dna, 78, len(dna))))


