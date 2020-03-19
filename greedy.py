from przydatne_funkcje import Score, profil, podobienstwo

def greedy_motif_search(DNA, k, t):
    """"Funkcja zwracająca motyw otrzymany poprzez zastosowanie alorytmu zachlannego
       DNA - tabelka sekwencji nukleotydowych
       k - dlugosc meru
       t - ilosc sekwencji w dna"""

    bestMotifs = []
    for x in range(0, t):                  #Dla każdej sekwencji w dna
        bestMotifs.append(0)               #nadaj indeks 0
    for n in range(0, (len(DNA[0])-k)):    #Od 0 do długości sekwencji w dna pomniejszonej o długość meru
        Motif = [n]                        #Indeks meru w pierwszej sekwecji
        for i in range(1, t):              #Pętla od 1 do ilosci sekwencji
            macierz = profil(Motif, DNA, k)                     #wykorzystanie funkcji profil zwracajacej macierz wystapien nukleotydow dla merow
            Motif.append(podobienstwo(macierz, DNA, i, k))      #Dodajemy kolejne indeksy korzystajac z funkcji zwracajacej indeks najbardziej podobnego meru
        if(Score(Motif, DNA, k)>Score(bestMotifs, DNA, k)):     #Jesli Score od motywu wiekszy niz Score od BestMotywu
            bestMotifs = Motif                                  #BestMotyw mprzyjmuje wartosc motywu
    return bestMotifs

if __name__ == "__main__":
    dna = ['acgtatcttc', 'tttatcggg', 'gactataaaa']
    print("\nDziałanie greedy_motif_search dla dna: ", dna)
    best = greedy_motif_search(dna, 3, len(dna))
    print("Best Motif: ",best)
    print("Score: ", Score(best, dna, 3))

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
    print("Działanie greedy_motif_search dla dna: ", dna)
    best = greedy_motif_search(dna, 30, len(dna))
    print("Best Motif: ", best)
    print("Score: ", Score(best, dna, 30))